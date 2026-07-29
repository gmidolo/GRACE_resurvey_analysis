################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
################################################################################
# Description:
#  Observation-error sensitivity analysis for Jaccard dissimilarity.
#  Q1: How much pseudoturnover (from imperfect species detection) inflate 
#      Jaccard?
#  Q1: How much does additional imperfect species detection alter the regression 
#      coefficients (protected areas, elevation, etc)?
#
# To add/inflate error (imperfect detection) - see add_obs_error() function:
# We apply two independent error sources applied to each survey visit:
# (1) Pure false negatives (p_fn): present species entirely missed; overlooked 
#     species in the plot
# (2) Species incorrect identification (misidentification) (p_misid): this is 
#     the percentage of present species recorded as a different absent taxon 
#     (though could also be misidentified as an already present taxon); we do a 
#     1-for-1 swap driven by species misidentified to another taxa. Each swap 
#     removes a real species and adds a wrong one (taken randomly from the full 
#     pool of species). e.g. p_misid = 0.1 means 1 species is misidentified in a 
#     plot with 10 species.
################################################################################

# 0. Setup ####

# set seed
set.seed(123)

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(glmmTMB)
  library(broom.mixed)
})

pth2save <- 'data/output_pseudoturnover_sim'
pth2fig <- 'fig'
if(!dir.exists(pth2save)){dir.create(pth2save)}
if(!dir.exists(pth2fig)){dir.create(pth2fig)}

# 1. Load data ####

# species data (vascular plants only, no genus-level records)
spd <- read_csv('data/raw_veg_data/EVA&GRACE_species_20260202.csv', show_col_types = F) %>%
  filter(taxongroup == 'Vascular plant') %>%
  filter(!str_detect(species, ' species'))

# modeling dataset: 411 plots with all predictors already computed
mod_dat <- read_csv('data/model_data/model_input_data.csv', show_col_types = F) %>%
  mutate(
    delta_CA = grassland_ca_2023 - grassland_ca_1870,
    delta_NP = grassland_np_2023 - grassland_np_1870,
    log_area_ratio = log10(rsrv_plot_size / historic_plot_size),
    historic_plot_size_log = log10(historic_plot_size),
    protected = factor(pa_site, levels = c(0, 1), labels = c('unprotected', 'protected')),
    protected_restricted = factor(pa_site_restrict, levels = c(0, 1), labels = c('unprotected', 'protected')),
    country = factor(country, levels = c('AT', 'CZSK'))
  )

plot_ids_final <- mod_dat$plot_id
spd <- spd %>% filter(plot_id %in% plot_ids_final)

# useful summary
cat('Modeling plots:', length(plot_ids_final), '\n')
cat('Total species records:', nrow(spd), '\n')
cat('Unique species:', n_distinct(spd$species), '\n')


# 2. Build paired P/A matrices ####

# separate historical from resurvey data
make_pa_wide <- function(resurvey_flag) {
  spd %>%
    filter(resurvey == resurvey_flag, plot_id %in% plot_ids_final) %>%
    select(plot_id, species) %>%
    distinct() %>%
    mutate(pa = 1L) %>%
    pivot_wider(names_from = species, values_from = pa, values_fill = 0L) %>%
    arrange(match(plot_id, plot_ids_final)) %>%
    column_to_rownames('plot_id')
}

pa_hist <- make_pa_wide(0)
pa_rsrv <- make_pa_wide(1)

# align to full species union
all_species <- union(colnames(pa_hist), colnames(pa_rsrv))
pa_hist[, setdiff(all_species, colnames(pa_hist))] <- 0L
pa_rsrv[, setdiff(all_species, colnames(pa_rsrv))] <- 0L
pa_hist <- pa_hist[, all_species]
pa_rsrv <- pa_rsrv[, all_species]

# check/describe data
cat('PA matrix:', nrow(pa_hist), 'plots x', ncol(pa_hist), 'species\n')


# 3. Observation error function ####

# Description:
# p_fn = proportion of present species entirely missed (not detected at all).
#        e.g. a p_fn=0.10 with 20 species means 2 species are missed.
# p_misid = share of species present MISIDENTIFIED per survey visit.
add_obs_error <- function(pa_vec, p_fn, p_misid) {
  # pa_vec = r
  # p_fn = fn_null
  # p_misid = misid_null
  y <- as.integer(pa_vec)
  present_idx <- which(y == 1L)
  absent_idx <- which(y == 0L)
  
  s <- length(present_idx) # calc species richness 
  n_misid <- floor(s*p_misid) # calc no. species misidentified

  # 1) Add false negatives — present species entirely missed
  if (length(present_idx) > 0L) {
    fn <- rbinom(length(present_idx), 1L, p_fn)
    y[present_idx[fn == 1L]] <- 0L
  } # (a) removes the real species (linked false negative)

  # 2) Add misidentifications — swap a no. species (=n_misid) present species to random absent ones
  # re-derive after false negatives so we do nott swap already-missed species
  # adds randomly chosen absent species (linked false positive) yet this is a 1-for-1 species swap, 
  # this is NOT a proper 'false positive'
  # (N.B.  this would mean adding fully new species that are actually absent; this is rarer).
    
  still_present <- which(y == 1L)
  still_absent  <- which(y == 0L)
  if (n_misid > 0L && length(still_present) > 0L && length(still_absent) > 0L) {
    n_swap <- min(n_misid, length(still_present), length(still_absent))
    swap_from <- sample(still_present, n_swap)   # real species removd
    swap_to <- sample(still_absent,  n_swap)   # wrong taxon added
    y[swap_from] <- 0L
    y[swap_to] <- 1L
  }

  return(y)
}


# 4. Null model: pseudoturnover floor ####

# Null model: here compare error on same resurvey plots list vs itself (zero true change).
n_null <- 999 # no. repeats
fn_null <- 0.10 # 10% of present species missed
misid_null <- 0.10 # 10% of species present misidentified per survey

null_jaccard <- map_dbl(seq_len(n_null), function(i) {
  idx <- sample(nrow(pa_rsrv), 1)
  r <- as.integer(pa_rsrv[idx, ])
  r1 <- add_obs_error(r, fn_null, misid_null)
  r2 <- add_obs_error(r, fn_null, misid_null)
  mat <- rbind(r1, r2) # compare same resurvey plots each with their own sources of introduced errors
  # average pseudoturnover should be lower if comparison is made with original releve (`mat <- rbind(r, r1)`)
  if (sum(mat) == 0L) return(NA_real_)
  as.numeric(vegdist(mat, method = 'jaccard', binary = TRUE))
})

# N.B., we interpret this as the average share of Jaccard attributable to pseudoturnover(Q1)
cat('Mean pseudoturnover Jaccard among resurvey plots (w repeat):', round(mean(null_jaccard), 3), '\n') 
cat('95th percentile pseudoturnover:', round(quantile(null_jaccard, 0.95), 3), '\n')
# This is actual (observed) Jaccard
cat('Observed mean Jaccard (GRACE data):', round(mean(mod_dat$jaccard), 3), '\n')

# 5. Jaccard inflation across p_fn x p_misid grid scenarios ####
# Grid 'search' to see how introduced uncertainty responds to various falsenegatives / misidentification: 
# We compute Jaccard inflation across p_fn × p_misid scenarios (both error sources applied to historic and resurvey)
error_grid <- expand.grid(
  p_fn    = c(0.05, 0.10, 0.15, 0.20),
  p_misid = c(0.05, 0.10, 0.20)
)

n_boot <- 999

grid_error <- error_grid %>%
  mutate(mean_inflation = NA_real_,
         intercept = NA_real_,
         slope  = NA_real_,
         rmse = NA_real_)

# Export figure
# graphics.off()
st = Sys.time()
png(
  filename = file.path(pth2fig, 'SI_jaccard.simulation_error_grid.png'), 
  width = 9, height = 6.75, res = 500, units = 'in'
)
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
for (g in seq_len(nrow(grid_error))) {
  p_fn    <- grid_error$p_fn[g]
  p_misid <- grid_error$p_misid[g]

  jac_true <- rep(0, n_boot)
  jac_err  <- rep(0, n_boot)

  for (i in seq_len(n_boot)) {
    idx   <- sample(nrow(pa_hist), 1)
    h <- as.integer(pa_hist[idx, ]) # historical plot
    r <- as.integer(pa_rsrv[idx, ]) # resurvey plot
    h_err <- add_obs_error(h, p_fn, p_misid)
    r_err <- add_obs_error(r, p_fn, p_misid)

    mat_true <- rbind(h, r)
    mat_err  <- rbind(h_err, r_err)
    
    # compare historical vs. resurvey-like plots (ture; and each with their own added error)
    jac_true[i] <- if (sum(mat_true) == 0L) NA_real_ else
      as.numeric(vegdist(mat_true, method = 'jaccard', binary = TRUE))
    jac_err[i]  <- if (sum(mat_err)  == 0L) NA_real_ else
      as.numeric(vegdist(mat_err,  method = 'jaccard', binary = TRUE))
  }

  ok  <- !is.na(jac_true) & !is.na(jac_err) # check jaccard is computed in both
  fit <- lm(jac_true[ok] ~ jac_err[ok])
  
  # plot
  plot(jac_err[ok], jac_true[ok], xlab = 'Jaccard with Error', ylab = 'True Jaccard', 
       main = paste0('FN: ', grid_error$p_fn[g]*100,'%; ', 'MISID: ', grid_error$p_misid[g]*100, '%'), xlim=c(0,1), ylim=c(0,1))
  abline(a=0, b=1, lty=2, col = 'grey') 
  abline(lm(jac_true[ok] ~ jac_err[ok]), lty = 1, col = 'brown1')
  text(0.25, 0.85, paste0('RMSE = ', round(sqrt(mean((jac_true[ok] - jac_err[ok])^2)), 3)))
  
  # store res
  grid_error$mean_inflation[g] <- mean(jac_err[ok] - jac_true[ok])
  grid_error$intercept[g] <- coef(fit)[1]
  grid_error$slope[g] <- coef(fit)[2]
  grid_error$rmse[g] <- sqrt(mean((jac_true[ok] - jac_err[ok])^2))
}
dev.off()
Sys.time() - st

# results
print(grid_error, digits = 3)

# pretty heatmap
p_inflation <- grid_error %>%
  mutate(p_fn = factor(p_fn),
         p_misid = factor(p_misid)) %>%
  ggplot(aes(x = p_fn, y = p_misid, fill = mean_inflation)) +
  geom_tile(color = 'white', linewidth = 0.8) +
  geom_text(aes(label = round(mean_inflation, 3)), size = 3.5) +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1,
                       name = 'Mean Jaccard\ninflation') +
  labs(x = '% missed species',
       y = '% misidentified species',
       title = 'Jaccard inflation due to observation error',
       subtitle = 'Based on 999 random draws of GRACE paired plots (n = 411)') +
  theme_bw(base_size = 12)
p_inflation


# 6. Coefficient sensitivity: refit Jaccard models under each error scenario ####

# Here we refit M2 (using glmmTMB ordbeta; to save time) per scenario and compare slope stability.
f_m2_sim <- jaccard_sim ~ 
  scale(delta_CA) +
  scale(delta_NP) +
  protected_restricted +
  scale(gs_temp_change_decade) +
  scale(elevation) +
  scale(log_area_ratio) +
  scale(dist_meters) +
  country +
  scale(timespan)

key_terms <- c(
  'scale(delta_CA)',
  'scale(delta_NP)',
  'protected_restrictedprotected',
  'scale(gs_temp_change_decade)',
  'scale(elevation)',
  'scale(timespan)',
  'scale(log_area_ratio)',
  'scale(dist_meters)',
  'countryCZSK'
)

# define function to fit model (use glmmTMB to save time vs. brms, using ordered beta reg)
fit_m2_sim <- function(dat) {
  tryCatch(
    glmmTMB(f_m2_sim, family = ordbeta(), data = dat, REML = FALSE),
    error   = function(e) { message('  glmmTMB error: ', e$message); NULL },
    warning = function(w) { message('  glmmTMB warning: ', w$message)
      glmmTMB(f_m2_sim, family = ordbeta(), data = dat, REML = FALSE) }
  )
}

# baseline: observed Jaccard (no added error)
coef_list <- list()
m2_base <- fit_m2_sim(mod_dat %>% mutate(jaccard_sim = jaccard)) # model fitted on true Jaccard values
if (!is.null(m2_base)) {
  coef_list[['baseline']] <- tidy(m2_base, conf.int = TRUE, effects = 'fixed') %>%
    filter(term %in% key_terms) %>%
    mutate(scenario = 'Baseline\n(observed)', p_fn = NA_real_, p_misid = NA_real_)
}

# fit error scenarios
for (g in seq_len(nrow(error_grid))) {
  p_fn    <- error_grid$p_fn[g]
  p_misid <- error_grid$p_misid[g]
  scenario_id <- sprintf('FalseNegative=%.2f\nMisid=%.2f', p_fn, p_misid)
  cat('g: ', g, '/', nrow(error_grid), ':', sprintf('FalseNegative=%.2f  p_misid=%.2f', p_fn, p_misid), '\n')
  
  # compute jaccard on plots with introduced uncertainty
  jaccard_sim <- map_dbl(seq_len(nrow(pa_hist)), function(i) {
    h_err <- add_obs_error(as.integer(pa_hist[i, ]), p_fn, p_misid)
    r_err <- add_obs_error(as.integer(pa_rsrv[i, ]), p_fn, p_misid)
    mat   <- rbind(h_err, r_err)
    if (sum(mat) == 0L) return(NA_real_)
    as.numeric(vegdist(mat, method = 'jaccard', binary = TRUE))
  })

  sim_dat <- mod_dat %>%
    mutate(jaccard_sim = jaccard_sim) %>%
    filter(!is.na(jaccard_sim))

  # fit glmmTMB
  m2_sim <- fit_m2_sim(sim_dat)
  
  if (!is.null(m2_sim)) {
    coef_list[[scenario_id]] <- tidy(m2_sim, conf.int = TRUE, effects = 'fixed') %>%
      filter(term %in% key_terms) %>%
      mutate(scenario = scenario_id, p_fn = p_fn, p_misid = p_misid)
  }
}

# tidy results
coef_df <- bind_rows(coef_list) %>%
  mutate(
    term_clean = term %>%
      str_remove_all('scale\\(|\\)') %>%
      str_replace('protectedprotected', 'Protected') %>%
      str_replace('countryCZSK', 'Country: CZSK') %>%
      str_replace_all('_', ' ') %>%
      str_remove('historic ') %>%
      str_to_lower(),
    is_baseline = (scenario == 'Baseline\n(observed)'),
    scenario = scenario %>% str_replace('Negative', 'Neg.') %>% str_replace('Misid', 'Misidentif.')
  ) 
unique(coef_df$term_clean)
new_term_d <- data.frame(
  term_clean = unique(coef_df$term_clean),
  new_term = c(
    'ΔCore Area',
    'ΔNo. Patches',
    'Protected site',
    'Warming',
    'Elevation',
    'ΔPlot size',
    'Relocation distance',
    'Country (CZE-SVK)',
    'Timespan'
  )
)
coef_df <- coef_df %>%
  left_join(new_term_d)

# coefficient stability plot
p_coefs <- ggplot(coef_df,
                  aes(x = estimate, y = scenario,
                      colour = is_baseline, alpha = is_baseline)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey60') +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                orientation = 'y', width = 0.35, linewidth = 0.6) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c('FALSE' = 'darkorange', 'TRUE' = 'royalblue'),
                      guide = 'none') +
  scale_alpha_manual(values = c('FALSE' = 0.65, 'TRUE' = 1.0), guide = 'none') +
  facet_wrap(~ new_term, ncol = 3) +
  labs(x = 'Coefficient estimate (ordbeta link scale)',
       y = NULL,
       title = 'Regression coefficients under simulated observation error',
       subtitle = 'blue = baseline model (observed Jaccard); orange = error-added models') +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        axis.text.y = element_text(size = 7))

p_coefs

ggsave(filename = file.path(pth2fig, 'SI_jaccard.simulation_regression.png'), p_coefs, width = 6, height = 9.75, dpi = 500)

#7. Summary & export ####
# `baseline` estimates (model on observed Jaccard)
baseline_lookup <- coef_df %>%
  filter(is_baseline) %>%
  select(term_clean, baseline_est = estimate)

# check how which predictors have their sign always consistent to the baseline
coef_df %>%
  filter(!is_baseline) %>%
  left_join(baseline_lookup, by = 'term_clean') %>%
  group_by(term_clean) %>%
  summarise(
    baseline = round(first(baseline_est), 4),
    min_est = round(min(estimate), 4),
    max_est = round(max(estimate), 4),
    range = round(max(estimate) - min(estimate), 4),
    sign_stable = all(sign(estimate) == sign(first(baseline_est))),
    .groups = 'drop'
  ) %>%
  print()

# export
write_csv(grid_error, file.path(pth2save,'jaccard_error_inflation.csv'))
write_csv(coef_df, file.path(pth2save, 'jaccard_coef_sensitivity.csv'))