################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
################################################################################
# Description: Plot figures of brms models (main and SI figures)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(bayestestR)
  library(brms)
  library(cowplot)
})
source('src/utils.R')

if(!dir.exists('fig')){dir.create('fig')}

# load data
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

# load main models
pth2models <- 'data/output_brms_models'
m1 <- read_rds(brms_name('m1_lnRR_S.rds'))
m2 <- read_rds(brms_name('m2_jaccard.rds'))
m3 <- read_rds(brms_name('m3_CM_N.rds'))
m4 <- read_rds(brms_name('m4_CM_M.rds'))
m1i <- read_rds(brms_name('m1i_lnRR_S.rds'))
m2i <- read_rds(brms_name('m2i_jaccard.rds'))
m3i <- read_rds(brms_name('m3i_CM_N.rds'))
m4i <- read_rds(brms_name('m4i_CM_M.rds'))

# 1. Response Distributions ####

panel_labels <- c(
  'Species Richness Change (logRRS)' = 'Species~Richness~Change~(logRR[S])',
  'Compositional Dissimilarity (Jaccard)' = 'Compositional~Dissimilarity~(Jaccard)',
  'Mean Nutrient Indicator Change (ΔCMN)' = 'Mean~Nutrient~Indicator~Change~(Delta*CM[N])',
  'Mean Moisture Indicator Change (ΔCMM)' = 'Mean~Moisture~Indicator~Change~(Delta*CM[M])'
)

hist_data_raw <- bind_rows(
  dat %>% transmute(panel = 'Species Richness Change (logRRS)', value = logRR_S),
  dat %>% transmute(panel = 'Compositional Dissimilarity (Jaccard)', value = jaccard),
  dat %>% transmute(panel = 'Mean Nutrient Indicator Change (ΔCMN)', value = diff.CM_N),
  dat %>% transmute(panel = 'Mean Moisture Indicator Change (ΔCMM)', value = diff.CM_M)
) %>%
  mutate(panel = factor(panel, levels = names(panel_labels)))

raw_stats <- hist_data_raw %>%
  group_by(panel) %>%
  summarise(mean_ci(value), .groups = 'drop')

pHist_raw <- ggplot(hist_data_raw, aes(value)) +
  geom_histogram(fill = 'grey85', col = 'grey70', bins = 20) +
  geom_rect(data = raw_stats, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = '95% CI'),
            inherit.aes = FALSE, alpha = 0.35, color = NA) +
  geom_vline(data = raw_stats, aes(xintercept = mean, color = 'Mean'), linewidth = 1) +
  facet_wrap(~panel, scales = 'free', ncol = 2,
             labeller = labeller(panel = as_labeller(panel_labels, label_parsed))) +
  scale_x_continuous(breaks = scales::breaks_pretty(4)) +
  scale_color_manual(NULL, values = c('Mean' = '#2166AC')) +
  scale_fill_manual(NULL, values = c('95% CI' = '#74ADD1')) +
  labs(x = NULL, y = 'Count (No. Plots)') +
  theme_bw() +
  theme(
    strip.text = element_text(face = 'bold'),
    legend.position = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha('white', 0.8), color = NA)
  )

pHist_raw

ggsave('fig/fig.2_hist_raw.pdf', pHist_raw, width = 6, height = 4.75, device = cairo_pdf)

# 2. Main Model Coefficients ####
m1r2 <- round(bayes_R2(m1), 2)
m2r2 <- round(bayes_R2(m2), 2)
m3r2 <- round(bayes_R2(m3), 2)
m4r2 <- round(bayes_R2(m4), 2)

pg_simplemodels <- plot_grid(
  plot_brms_coefs(m1, title = expression(Species~Richness~Change~(logRR[S]))) + 
    labs(subtitle = bquote(R^2 == .(m1r2[1, 'Estimate']) %+-% .(m1r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m2, title = expression(Compositional~Dissimilarity~(Jaccard))) + 
    labs(subtitle = bquote(R^2 == .(m2r2[1, 'Estimate']) %+-% .(m2r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m3, title = expression(Mean~Nutrient~Indicator~Change~(Delta*CM[N]))) + 
    labs(subtitle = bquote(R^2 == .(m3r2[1, 'Estimate']) %+-% .(m3r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m4, title = expression(Mean~Moisture~Indicator~Change~(Delta*CM[M]))) + 
    labs(subtitle = bquote(R^2 == .(m4r2[1, 'Estimate']) %+-% .(m4r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none')
)
pg_simplemodels

ggsave('fig/fig.3_mainmodels_raw.pdf', pg_simplemodels, width = 9, height = 6, device = cairo_pdf)

# 3. Conditional Interactions ####

piA <- plot_brms_condeffect(
  model = m1i, effect = 'timespan:gs_temp_change_decade', dat = dat,
  x_var = 'gs_temp_change_decade', group_var = 'timespan',
  x_label = 'Warming 1979-2021 (°C/decade)', y_label = expression(logRR[S]),
  group_label = 'Timespan (years)', new_pal = c('#BAE4B3', '#008f3a', '#003315'),
  alpha.point = .4, title = 'Species Richness Change', subtitle = 'Warming × Timespan'
) + geom_hline(yintercept = 0, linetype = 2, color = 'grey50', linewidth = 0.4)

piB <- plot_brms_condeffect(
  model = m2i, effect = 'timespan:protected_restricted', dat = dat,
  x_var = 'timespan', group_var = 'protected_restricted',
  x_label = 'Timespan (years)', y_label = 'Jaccard index',
  group_label = 'Protected areas', title = 'Compositional Dissimilarity',
  subtitle = 'Protected areas × Timespan', alpha.point = .4
)

piC <- plot_brms_condeffect(
  model = m4i, effect = 'timespan:protected_restricted', dat = dat,
  x_var = 'timespan', group_var = 'protected_restricted',
  x_label = 'Timespan (years)', y_label = expression(Delta * CM[M]),
  group_label = 'Protected areas', title = 'Mean Moisture Indicator Change',
  subtitle = 'Protected areas × Timespan', alpha.point = .4
) + geom_hline(yintercept = 0, linetype = 2, color = 'grey50', linewidth = 0.4)

pi_grd <- plot_grid(piA, piB, piC, labels = LETTERS[1:3], nrow = 1)
pi_grd
ggsave('fig/fig.4_interaction_scatterplot.pdf', f4pi, device = cairo_pdf, width = 10.5, height = 4)

# 4. Supplementary Figures ####

# SI: Interaction Models Overview
m1ir2 <- round(bayes_R2(m1i), 2)
m2ir2 <- round(bayes_R2(m2i), 2)
m3ir2 <- round(bayes_R2(m3i), 2)
m4ir2 <- round(bayes_R2(m4i), 2)

pgsimplemodels_interaction <- plot_grid(
  plot_brms_coefs(m1i, title = expression(Species~Richness~Change~(logRR[S]))) + 
    labs(subtitle = bquote(R^2 == .(m1ir2[1, 'Estimate']) %+-% .(m1ir2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m2i, title = expression(Compositional~Dissimilarity~(Jaccard))) + 
    labs(subtitle = bquote(R^2 == .(m2ir2[1, 'Estimate']) %+-% .(m2ir2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m3i, title = expression(Mean~Nutrient~Indicator~Change~(Delta*CM[N]))) + 
    labs(subtitle = bquote(R^2 == .(m3ir2[1, 'Estimate']) %+-% .(m3ir2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m4i, title = expression(Mean~Moisture~Indicator~Change~(Delta*CM[M]))) + 
    labs(subtitle = bquote(R^2 == .(m4ir2[1, 'Estimate']) %+-% .(m4ir2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  ncol = 2
)
pgsimplemodels_interaction
ggsave('fig/SI_interaction.pdf', pgsimplemodels_interaction, width = 11, height = 9, device = cairo_pdf)

# SI: Models with data sampled post 1975 ####
m1_post75 <- read_rds(brms_name('m1_lnRR_S_post75.rds'))
m2_post75 <- read_rds(brms_name('m2_jaccard_post75.rds'))
m3_post75 <- read_rds(brms_name('m3_CM_N_post75.rds'))
m4_post75 <- read_rds(brms_name('m4_CM_M_post75.rds'))

m1_p75_r2 <- round(bayes_R2(m1_post75), 2)
m2_p75_r2 <- round(bayes_R2(m2_post75), 2)
m3_p75_r2 <- round(bayes_R2(m3_post75), 2)
m4_p75_r2 <- round(bayes_R2(m4_post75), 2)

pg_post75 <- plot_grid(
  plot_brms_coefs(m1_post75, title = expression(Species~Richness~Change~(logRR[S]))) + 
    labs(subtitle = bquote(R^2 == .(m1_p75_r2[1, 'Estimate']) %+-% .(m1_p75_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m2_post75, title = expression(Compositional~Dissimilarity~(Jaccard))) + 
    labs(subtitle = bquote(R^2 == .(m2_p75_r2[1, 'Estimate']) %+-% .(m2_p75_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m3_post75, title = expression(Mean~Nutrient~Indicator~Change~(Delta*CM[N]))) + 
    labs(subtitle = bquote(R^2 == .(m3_p75_r2[1, 'Estimate']) %+-% .(m3_p75_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m4_post75, title = expression(Mean~Moisture~Indicator~Change~(Delta*CM[M]))) + 
    labs(subtitle = bquote(R^2 == .(m4_p75_r2[1, 'Estimate']) %+-% .(m4_p75_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  ncol = 2
)
pg_post75
ggsave('fig/SI_post1975_models.pdf', pg_post75, width = 11, height = 9, device = cairo_pdf)

# SI: Models with Broadly defined protected areas ####
m1_pabd <- read_rds(brms_name('m1_lnRR_S_pabd.rds'))
m2_pabd <- read_rds(brms_name('m2_jaccard_pabd.rds'))
m3_pabd <- read_rds(brms_name('m3_CM_N_pabd.rds'))
m4_pabd <- read_rds(brms_name('m4_CM_M_pabd.rds'))

m1_pabd_r2 <- round(bayes_R2(m1_pabd), 2)
m2_pabd_r2 <- round(bayes_R2(m2_pabd), 2)
m3_pabd_r2 <- round(bayes_R2(m3_pabd), 2)
m4_pabd_r2 <- round(bayes_R2(m4_pabd), 2)

pg_pabd <- plot_grid(
  plot_brms_coefs(m1_pabd, title = expression(Species~Richness~Change~(logRR[S]))) + 
    labs(subtitle = bquote(R^2 == .(m1_pabd_r2[1, 'Estimate']) %+-% .(m1_pabd_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m2_pabd, title = expression(Compositional~Dissimilarity~(Jaccard))) + 
    labs(subtitle = bquote(R^2 == .(m2_pabd_r2[1, 'Estimate']) %+-% .(m2_pabd_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m3_pabd, title = expression(Mean~Nutrient~Indicator~Change~(Delta*CM[N]))) + 
    labs(subtitle = bquote(R^2 == .(m3_pabd_r2[1, 'Estimate']) %+-% .(m3_pabd_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  plot_brms_coefs(m4_pabd, title = expression(Mean~Moisture~Indicator~Change~(Delta*CM[M]))) + 
    labs(subtitle = bquote(R^2 == .(m4_pabd_r2[1, 'Estimate']) %+-% .(m4_pabd_r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none'),
  ncol = 2
)
pg_pabd
ggsave('fig/SI_pabd_models.pdf', pg_pabd, width = 11, height = 9, device = cairo_pdf)

# SI: Ellenberg Temperature Change (M5)
m5 <- read_rds(brms_name('m5_CM_T.rds'))
m5i <- read_rds(brms_name('m5i_CM_T.rds'))

pHist_temperature <- ggplot(dat, aes(diff.CM_T)) +
  geom_histogram(fill = 'grey85', col = 'grey70', bins = 20) +
  geom_rect(data = mean_ci(dat$diff.CM_T), aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = '95% CI'),
            inherit.aes = FALSE, alpha = 0.35, color = NA) +
  geom_vline(data = mean_ci(dat$diff.CM_T), aes(xintercept = mean, color = 'Mean'), linewidth = 1) +
  scale_x_continuous(breaks = scales::breaks_pretty(4)) +
  scale_color_manual(NULL, values = c('Mean' = '#2166AC')) +
  scale_fill_manual(NULL, values = c('95% CI' = '#74ADD1')) +
  labs(title = expression(bold(Delta*CM[`T`])), x = expression(Delta*CM[`T`]), y = 'Count') +
  theme_bw()

m5r2  <- round(bayes_R2(m5), 2)
m5ir2 <- round(bayes_R2(m5i), 2)
m5p1 <- plot_brms_coefs(m5, title = 'Main Model') +
  labs(subtitle = bquote(R^2 == .(m5r2[1, 'Estimate']) %+-% .(m5r2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none')
m5p2 <- plot_brms_coefs(m5i, title = 'Timespan Interaction Model') +
  labs(subtitle = bquote(R^2 == .(m5ir2[1, 'Estimate']) %+-% .(m5ir2[1, 'Est.Error']) ~ '(SE)')) + theme(legend.position = 'none')

m5_grid_title <- ggdraw() + draw_label(expression(bold(Mean~Temperature~Indicator~Change~(Delta*CM[`T`]))), fontface = 'bold', x = 0.5, hjust = 0.5, size = 16)
m5_plots_row  <- plot_grid(pHist_temperature, m5p1, m5p2, nrow = 1, labels = LETTERS[1:3], rel_widths = c(1, 1, 1.3))
m5p_final     <- plot_grid(m5_grid_title, m5_plots_row, ncol = 1, rel_heights = c(0.1, 1))

ggsave('fig/SI_EllenbergTemperature.pdf', m5p_final, width = 12.5, height = 4)