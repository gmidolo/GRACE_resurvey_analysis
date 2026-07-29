################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
################################################################################
# Description: Temporal analysis of biotic homogenization dynamics across 
#              time slices and country x decade subgroups.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
})

# 1. Data preparation ####

# load data & slice decades
dat <- read_csv('data/model_data/model_input_data.csv', show_col_types = FALSE) %>%
  mutate(
    delta_CA = grassland_ca_2023 - grassland_ca_1870,
    delta_NP = grassland_np_2023 - grassland_np_1870,
    log_area_ratio = log10(rsrv_plot_size / historic_plot_size),
    historic_plot_size_log = log10(historic_plot_size),
    protected = factor(pa_site, levels = c(0, 1), labels = c('unprotected', 'protected')),
    protected_restricted = factor(pa_site_restrict, levels = c(0, 1), labels = c('unprotected', 'protected')),
    country = factor(country, levels = c('AT', 'CZSK'))
  ) %>%
  mutate(
    sampling_year_slice = case_when(
      historic_sampling_year < 1960               ~ '<1960',
      between(historic_sampling_year, 1960, 1969) ~ '1960s',
      between(historic_sampling_year, 1970, 1979) ~ '1970s',
      between(historic_sampling_year, 1980, 1989) ~ '1980s',
      between(historic_sampling_year, 1990, 1999) ~ '1990s',
      between(historic_sampling_year, 2000, 2009) ~ '2000s',
      between(historic_sampling_year, 2010, 2020) ~ '2010s',
      TRUE ~ NA_character_
    )
  )

# load species data
spd <- read_csv('data/raw_veg_data/EVA&GRACE_species_20260202.csv', show_col_types = FALSE) %>%
  filter(
    taxongroup == 'Vascular plant',
    !str_detect(species, ' species'),
    plot_id %in% dat$plot_id
  ) %>%
  select(plot_id, resurvey, species) %>%
  distinct()

# species richness data
plot_S <- spd %>%
  group_by(plot_id, resurvey) %>%
  summarise(S = n_distinct(species), .groups = 'drop') %>%
  pivot_wider(names_from = resurvey, values_from = S, names_prefix = 'S_') %>%
  rename(hist_S = S_0, rsrv_S = S_1) %>%
  left_join(dat %>% select(plot_id, timespan), by = 'plot_id')

# 2. Alpha and Gamma calculation functions ####

compute_metrics <- function(plot_ids, sp_long, pl_rich, method = 'per_plot') {
  method  <- match.arg(method, choices = c('aggregate', 'per_plot'))
  pr_boot <- tibble(plot_id = plot_ids) %>% left_join(pl_rich, by = 'plot_id')
  
  alpha_hist <- mean(pr_boot$hist_S, na.rm = TRUE)
  alpha_rsrv <- mean(pr_boot$rsrv_S, na.rm = TRUE)
  duration   <- mean(pr_boot$timespan, na.rm = TRUE)
  
  sp_sub     <- sp_long %>% filter(plot_id %in% unique(plot_ids))
  gamma_hist <- sp_sub %>% filter(resurvey == 0) %>% pull(species) %>% n_distinct()
  gamma_rsrv <- sp_sub %>% filter(resurvey == 1) %>% pull(species) %>% n_distinct()
  
  delta_alpha <- if (method == 'aggregate') {
    log10(alpha_rsrv / alpha_hist) / duration
  } else {
    pr_boot %>% 
      mutate(rate = log10(rsrv_S / hist_S) / timespan) %>%
      pull(rate) %>% 
      mean(na.rm = TRUE)
  }
  delta_gamma <- log10(gamma_rsrv / gamma_hist) / duration
  
  tibble(
    method, n_plots = n_distinct(plot_ids), duration = round(duration, 2),
    alpha_hist = round(alpha_hist, 3), alpha_rsrv = round(alpha_rsrv, 3),
    gamma_hist, gamma_rsrv,
    delta_alpha = round(delta_alpha, 6),
    delta_gamma = round(delta_gamma, 6),
    delta_beta  = round(delta_gamma - delta_alpha, 6)
  )
}

bootstrap_metrics <- function(plot_ids, sp_long, pl_rich, n_boot = 999, method = 'per_plot') {
  map_dfr(seq_len(n_boot), function(i) {
    compute_metrics(sample(plot_ids, length(plot_ids), replace = TRUE), sp_long, pl_rich, method) %>%
      select(delta_alpha, delta_gamma, delta_beta)
  })
}

run_hom <- function(subgroups, method = 'per_plot', n_boot = 999, seed = 123) {
  set.seed(seed)
  map_dfr(names(subgroups), function(g) {
    ids  <- subgroups[[g]]
    pt   <- compute_metrics(ids, spd, plot_S, method) %>% mutate(group = g)
    boot <- bootstrap_metrics(ids, spd, plot_S, n_boot, method)
    ci   <- boot %>% summarise(
      delta_alpha_lo = quantile(delta_alpha, 0.025), delta_alpha_hi = quantile(delta_alpha, 0.975),
      delta_gamma_lo = quantile(delta_gamma, 0.025), delta_gamma_hi = quantile(delta_gamma, 0.975),
      delta_beta_lo  = quantile(delta_beta, 0.025),  delta_beta_hi  = quantile(delta_beta, 0.975)
    )
    bind_cols(pt, ci)
  })
}

delta_method <- 'per_plot'
n_boot_all   <- 999

# 3. Part 1: Historical Time Slices Analysis ####

ts_levels <- c('<1960', '1960s', '1970s', '1980s', '1990s', '2000s', '2010s')

subgroups_ts <- dat %>%
  filter(!is.na(sampling_year_slice)) %>%
  split(.$sampling_year_slice) %>%
  map(~ .x$plot_id) %>%
  .[ts_levels[ts_levels %in% names(.)]]

# bootstrap homogenization/differentiation analyses
st = Sys.time()
hom_ts <- run_hom(subgroups_ts, delta_method, n_boot_all, seed = 123)
print(Sys.time()-st) # ~ 1.5 minutes

ts_n   <- map_int(subgroups_ts, length)
ts_lbl <- paste0(names(ts_n), ' (', ts_n, ')')
names(ts_lbl) <- names(ts_n)

plot_dat_ts <- hom_ts %>%
  mutate(label = factor(ts_lbl[group], levels = ts_lbl))

clr_ts <- setNames(
  rev(viridis(nlevels(plot_dat_ts$label) + 1, option = 'plasma', direction = 1))[-1],
  levels(plot_dat_ts$label)
)

ax_ts <- plot_dat_ts %>%
  select(delta_alpha, delta_gamma, delta_alpha_lo, delta_alpha_hi, delta_gamma_lo, delta_gamma_hi) %>%
  unlist() %>% range(na.rm = TRUE) %>%
  (function(r) c(floor(r[1] * 1e4) / 1e4, ceiling(r[2] * 1e4) / 1e4))

p_typology_ts <- plot_dat_ts %>%
  ggplot(aes(x = delta_alpha, y = delta_gamma, colour = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_hline(yintercept = 0, colour = 'grey70', linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = 'grey70', linewidth = 0.4) +
  geom_errorbar(aes(ymin = delta_gamma_lo, ymax = delta_gamma_hi), width = 0, alpha = 0.75) +
  geom_errorbar(aes(xmin = delta_alpha_lo, xmax = delta_alpha_hi), width = 0, alpha = 0.75, orientation = 'y') +
  geom_point(size = 4, shape = 16) +
  scale_colour_manual(values = clr_ts, name = 'Historical sampling\nperiod (no. plots)') +
  labs(x = expression(Delta*alpha ~ '(annual log10-ratio of mean local richness)'),
       y = expression(Delta*gamma ~ '(annual log10-ratio of regional richness)')) +
  theme_bw(base_size = 12)

p_beta_ts <- plot_dat_ts %>%
  ggplot(aes(x = delta_beta, y = label, colour = label)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(xmin = delta_beta_lo, xmax = delta_beta_hi), width = 0.3, alpha = 0.7, orientation = 'y') +
  geom_point(size = 3, shape = 16) +
  scale_colour_manual(values = clr_ts, name = NULL) +
  labs(x = expression(Delta*beta ~ '(annual log10-ratio)'), y = 'Historical sampling period') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none')

# 4. Part 2: Country × Decade (1970s–2010s) ####

ts_keep <- c('1970s', '1980s', '1990s', '2000s', '2010s')

subgroups_tsc <- dat %>%
  filter(sampling_year_slice %in% ts_keep) %>%
  mutate(grp = paste0(sampling_year_slice, '_', country)) %>%
  split(.$grp) %>%
  map(~ .x$plot_id)

grp_order <- paste0(rep(ts_keep, each = 2), '_', c('AT', 'CZSK'))
grp_order <- grp_order[grp_order %in% names(subgroups_tsc)]
subgroups_tsc <- subgroups_tsc[grp_order]

# bootstrap homogenization/differentiation analyses
st = Sys.time()
hom_tsc <- run_hom(subgroups_tsc, delta_method, n_boot_all, seed = 456) %>%
  mutate(
    decade = factor(str_extract(group, '^[^_]+'), levels = ts_keep),
    country = factor(str_extract(group, '[^_]+$'), levels = c('AT', 'CZSK')),
    n_grp = map_int(group, ~ length(subgroups_tsc[[.x]])),
    label = paste0(decade, ' (', n_grp, ')')
  ) %>%
  mutate(label = factor(label, levels = unique(label[order(decade)])))
print(Sys.time()-st) # ~ 2 minutes

clr_decade <- setNames(clr_ts[3:7], ts_keep)

p_beta_tsc <- hom_tsc %>%
  ggplot(aes(x = delta_beta, y = label, colour = decade)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(xmin = delta_beta_lo, xmax = delta_beta_hi), width = 0.3, alpha = 0.7, orientation = 'y') +
  geom_point(size = 3, shape = 16) +
  scale_colour_manual(values = clr_decade, name = 'Decade') +
  facet_wrap(~ country, scales = 'free_y', ncol = 1,
             labeller = as_labeller(c('AT' = 'Austria', 'CZSK' = 'Czechia-Slovakia'))) +
  labs(x = expression(Delta*beta ~ '(annual log10-ratio)'), y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position  = 'none', strip.text = element_text(face = 2, hjust = 0), strip.background = element_blank())

# 5. Export figures ####

pth2save <- './fig/'
if(!dir.exists(pth2save)){dir.create(pth2save)}

ggsave(paste0(pth2save, 'fig5A_alpha.vs_gamma.pdf'), p_typology_ts, width = 5.45, height = 4.15, device = cairo_pdf)
ggsave(paste0(pth2save, 'fig5A_beta.pdf'), p_beta_ts + theme(axis.title.y = element_blank()), width = 3, height = 4.15, device = cairo_pdf)
ggsave(paste0(pth2save, 'fig5B_beta_country.pdf'), p_beta_tsc, width = 3.8, height = 5, device = cairo_pdf)