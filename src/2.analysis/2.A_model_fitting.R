################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
################################################################################
# Description: Fits Bayesian GLMMs using brms for plant biodiversity change
#              metrics (logRR_S, Jaccard dissimilarity, ΔCM_N, ΔCM_M, ΔCM_T).
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
})
source('src/utils.R')

pth2models <- 'data/output_brms_models'
if(!dir.exists(pth2models)){dir.create(pth2models, recursive = T)}

# 1. Data Preparation ####

dat <- read_csv('data/model_data/model_input_data.csv', show_col_types = FALSE) %>%
  mutate(
    delta_CA = grassland_ca_2023 - grassland_ca_1870,
    delta_NP = grassland_np_2023 - grassland_np_1870,
    log_area_ratio = log10(rsrv_plot_size / historic_plot_size),
    historic_plot_size_log = log10(historic_plot_size),
    protected = factor(pa_site, levels = c(0, 1), labels = c('unprotected', 'protected')),
    protected_restricted = factor(pa_site_restrict, levels = c(0, 1), labels = c('unprotected', 'protected')),
    country = factor(country, levels = c('AT', 'CZSK'))
  )

# 2. Model Formulas & Priors ####

# priors
priors_gaussian <- c(prior(normal(0, 1), class = 'b'))
priors_zoib <- c(prior(normal(0, 1), class = 'b'))

# main predictors formula
f_rhs <- ~ scale(delta_CA) +
  scale(delta_NP) +
  protected_restricted +
  scale(gs_temp_change_decade) +
  scale(elevation) +
  country +
  scale(historic_plot_size_log) +
  scale(log_area_ratio) +
  scale(dist_meters) +
  scale(timespan)

# interaction model formula
f_rhs_i <- ~ scale(delta_CA) * scale(timespan) +
  scale(delta_NP) * scale(timespan) +
  protected_restricted * scale(timespan) +
  scale(gs_temp_change_decade) * scale(timespan) +
  scale(elevation) +
  country +
  scale(historic_plot_size_log) +
  scale(log_area_ratio) +
  scale(dist_meters)

# formula for broader protected area definition
f_pabd <- ~ scale(delta_CA) +
  scale(delta_NP) +
  protected +
  scale(gs_temp_change_decade) +
  scale(elevation) +
  country +
  scale(historic_plot_size_log) +
  scale(log_area_ratio) +
  scale(dist_meters) +
  scale(timespan)

# settings for brms
brm_settings <- list(
  chains = 4, iter = 4000, warmup = 1000, cores = 4, 
  seed = 123, file_refit = 'on_change'
)

# function for brms
fit_brm_model <- function(formula, family, prior, file_path, data_subset = dat) {
  rlang::exec(
    brm,
    formula = bf(formula),
    family = family,
    data = data_subset,
    prior = prior,
    file = file_path,
    !!!brm_settings
  )
}

# 3. Main Models ####

# M1: Species richness change (logRR_S)
m1  <- fit_brm_model(update(f_rhs, logRR_S ~ .), gaussian(), priors_gaussian, brms_name('m1_lnRR_S'))
m1i <- fit_brm_model(update(f_rhs_i, logRR_S ~ .), gaussian(), priors_gaussian, brms_name('m1i_lnRR_S'))

# M2: Compositional dissimilarity (Jaccard)
m2  <- fit_brm_model(update(f_rhs, jaccard ~ .), zero_one_inflated_beta(), priors_zoib, brms_name('m2_jaccard'))
m2i <- fit_brm_model(update(f_rhs_i, jaccard ~ .), zero_one_inflated_beta(), priors_zoib, brms_name('m2i_jaccard'))

# M3: Mean Nutrient Indicator change (ΔCM_N)
m3  <- fit_brm_model(update(f_rhs, diff.CM_N ~ .), gaussian(), priors_gaussian, brms_name('m3_CM_N'))
m3i <- fit_brm_model(update(f_rhs_i, diff.CM_N ~ .), gaussian(), priors_gaussian, brms_name('m3i_CM_N'))

# M4: Mean Moisture Indicator change (ΔCM_M)
m4  <- fit_brm_model(update(f_rhs, diff.CM_M ~ .), gaussian(), priors_gaussian, brms_name('m4_CM_M'))
m4i <- fit_brm_model(update(f_rhs_i, diff.CM_M ~ .), gaussian(), priors_gaussian, brms_name('m4i_CM_M'))

# 4. Sensitivity Models ####

# Broader protection definition models
m1_pabd <- fit_brm_model(update(f_pabd, logRR_S ~ .), gaussian(), priors_gaussian, brms_name('m1_lnRR_S_pabd'))
m2_pabd <- fit_brm_model(update(f_pabd, jaccard ~ .), zero_one_inflated_beta(), priors_zoib, brms_name('m2_jaccard_pabd'))
m3_pabd <- fit_brm_model(update(f_pabd, diff.CM_N ~ .), gaussian(), priors_gaussian, brms_name('m3_CM_N_pabd'))
m4_pabd <- fit_brm_model(update(f_pabd, diff.CM_M ~ .), gaussian(), priors_gaussian, brms_name('m4_CM_M_pabd'))

# Post-1975 subset models
dat_post75 <- dat %>% filter(historic_sampling_year >= 1975)
m1_post75 <- fit_brm_model(update(f_rhs, logRR_S ~ .), gaussian(), priors_gaussian, brms_name('m1_lnRR_S_post75'))
m2_post75 <- fit_brm_model(update(f_rhs, jaccard ~ .), zero_one_inflated_beta(), priors_zoib, brms_name('m2_jaccard_post75'))
m3_post75 <- fit_brm_model(update(f_rhs, diff.CM_N ~ .), gaussian(), priors_gaussian, brms_name('m3_CM_N_post75'))
m4_post75 <- fit_brm_model(update(f_rhs, diff.CM_M ~ .), gaussian(), priors_gaussian, brms_name('m4_CM_M_post75'))

# M5: Mean Temperature Indicator change (ΔCM_T)
m5 <- fit_brm_model(update(f_rhs, diff.CM_T ~ .), gaussian(), priors_gaussian, brms_name('m5_CM_T'))
m5i <- fit_brm_model(update(f_rhs_i, diff.CM_T ~ .), gaussian(), priors_gaussian, brms_name('m5i_CM_T'))

# Prior sensitivity: wider priors Normal(0, 2.5)
priors_gaussian_wide <- c(prior(normal(0, 2.5), class = 'b'))
priors_zoib_wide <- c(prior(normal(0, 2.5), class = 'b'))
m1_wide <- fit_brm_model(update(f_rhs, logRR_S ~ .), gaussian(), priors_gaussian_wide, brms_name('m1_lnRR_S_wideprior'))
m2_wide <- fit_brm_model(update(f_rhs, jaccard ~ .), zero_one_inflated_beta(), priors_zoib_wide, brms_name('m2_jaccard_wideprior'))
m3_wide <- fit_brm_model(update(f_rhs, diff.CM_N ~ .), gaussian(), priors_gaussian_wide, brms_name('m3_CM_N_wideprior'))
m4_wide <- fit_brm_model(update(f_rhs, diff.CM_M ~ .), gaussian(), priors_gaussian_wide, brms_name('m4_CM_M_wideprior'))