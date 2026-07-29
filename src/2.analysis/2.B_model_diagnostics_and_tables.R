################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
################################################################################
# Description: Model adequacy checks, posterior metrics, spatial autocorrelation 
#              (Moran's I), VIF collinearity tests, and summary tables export.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(bayestestR)
  library(spdep)
  library(car)
  library(patchwork)
})
source('src/utils.R')

# load data
dat <- read_csv('data/model_data/model_input_data.csv', show_col_types = FALSE) %>%
  mutate(
    delta_CA = grassland_ca_2023 - grassland_ca_1870,
    delta_NP = grassland_np_2023 - grassland_np_1870,
    log_area_ratio = log10(rsrv_plot_size / historic_plot_size),
    historic_plot_size_log = log10(historic_plot_size),
    protected_restricted = factor(pa_site_restrict, levels = c(0, 1), labels = c('unprotected', 'protected')),
    country = factor(country, levels = c('AT', 'CZSK'))
  )

# load main models
pth2models <- 'data/output_brms_models'
m1 <- read_rds(brms_name('m1_lnRR_S.rds'))
m2 <- read_rds(brms_name('m2_jaccard.rds'))
m3 <- read_rds(brms_name('m3_CM_N.rds'))
m4 <- read_rds(brms_name('m4_CM_M.rds'))
model_list <- list(
  'Species richness change (logRR_S)'       = m1,
  'Compositional dissimilarity (Jaccard)'   = m2,
  'Mean Nutrient Indicator change (ΔCM N)'  = m3,
  'Mean Moisture Indicator change (ΔCM M)'  = m4
)

# 1. Posterior predictive checks ####

plots <- imap(model_list, ~ pp_check(.x, ndraws = 100) + ggtitle(.y))
wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = 'Posterior predictive checks',
    theme = theme(plot.title = element_text(size = 16, face = 'bold'))
  )


# 2. Multicollinearity Check (VIF) ####

vif_check <- lm(
  logRR_S ~ scale(delta_CA) + scale(delta_NP) + protected_restricted +
    scale(gs_temp_change_decade) + scale(elevation) +
    scale(historic_plot_size_log) + scale(timespan) + scale(log_area_ratio) + 
    scale(dist_meters) + country,
  data = dat
)
print(vif(vif_check))


# 3. Residual Spatial Autocorrelation (Moran's I) ####

coords <- dat %>% select(historic_lon, historic_lat) %>% as.matrix()
nb <- knn2nb(knearneigh(coords, k = 10))
lw <- nb2listw(nb, style = 'W')

table_s1 <- map_dfr(names(model_list), function(mod_label) {
  mod <- model_list[[mod_label]]
  res <- residuals(mod, method = 'posterior_predict')[, 'Estimate']
  mt  <- moran.test(res, lw)
  tibble(
    `Response variable` = mod_label,
    `Moran's I` = round(mt$estimate[['Moran I statistic']], 4),
    `Expected I` = round(mt$estimate[['Expectation']], 4),
    `Variance` = round(mt$estimate[['Variance']], 6),
    `Z score` = round(mt$statistic, 3),
    `P value` = signif(mt$p.value, 3)
  )
})

table_s1

write_csv(table_s1, file.path(pth2models, 'table_s1_morans_I.csv'))


# 4. Main Results Table Export ####

extract_results <- function(model, model_name) {
  describe_posterior(model, centrality = 'median', ci = 0.95, ci_method = 'hdi', test = 'pd') %>%
    filter(str_detect(Parameter, '^b_')) %>%
    mutate(Model = model_name, .before = 1) %>%
    select(Model, Parameter, Median, CI_low, CI_high, pd) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

results_table <- bind_rows(
  extract_results(m1, names(model_list)[1]),
  extract_results(m2, names(model_list)[2]),
  extract_results(m3, names(model_list)[3]),
  extract_results(m4, names(model_list)[4])
)

results_table

write_csv(results_table, file.path(pth2models, 'main_results_table.csv'))