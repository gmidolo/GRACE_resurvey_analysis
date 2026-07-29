################################################################################
# Script: utils.R
# Author: Gabriele Midolo
# Description: Helper functions for posterior plotting, custom themes, 
#              summary statistics, and conditional effect visualizations.
################################################################################


# Coefficient plot for brms models (forest plot style)

plot_brms_coefs <- function(model,
                            exclude_intercept = TRUE,
                            exclude_pattern = NULL,
                            rename_labels = NULL,
                            order_by_effect = FALSE,
                            ref_line = 0,
                            widths = c(0.50, 0.95),
                            point_size = 2.5,
                            title = NULL) {
  
  # Posterior summaries
  summ <- as_draws_df(model) %>%
    select(starts_with('b_')) %>%
    pivot_longer(everything(), names_to = 'parameter', values_to = 'value') %>%
    group_by(parameter) %>%
    summarise(
      median    = median(value),
      lo_wide   = quantile(value, (1 - widths[2]) / 2),
      hi_wide   = quantile(value, 1 - (1 - widths[2]) / 2),
      lo_narrow = quantile(value, (1 - widths[1]) / 2),
      hi_narrow = quantile(value, 1 - (1 - widths[1]) / 2),
      pd        = max(mean(value > 0), mean(value < 0)),
      .groups   = 'drop'
    ) %>%
    mutate(
      evidence = case_when(
        pd >= 0.95 ~ 'strong (pd ≥ 0.95)',
        # pd >= 0.90 ~ 'moderate (pd ≥ 0.90)',
        TRUE       ~ 'weak (pd < 0.90)'
      )
    )
  
  # Use str_detect so ZOIB sub-model intercepts (b_mu_Intercept, etc.) are
  # also excluded when exclude_intercept = TRUE
  if (exclude_intercept)
    summ <- summ %>% filter(!str_detect(parameter, 'Intercept'))
  
  if (!is.null(exclude_pattern))
    summ <- summ %>% filter(!str_detect(parameter, exclude_pattern))
  
  #Label cleaning 
  term_lut <- c(
    'timespan'           = 'Timespan',
    'log_area_ratio'     = 'ΔPlot size',
    'historic_plot_size_log' = 'Base plot size',
    'dist_meters'        = 'Relocation distance',
    # ecological predictors:
    'elevation' = 'Elevation',
    'delta_CA'           = 'ΔCore Area',
    'delta_NP'           = 'ΔNo. Patches',
    'protected_restrictedprotected' = 'Protected site',
    'protectedprotected' = 'Protected site\n(all designations)',
    'gs_temp_change_decade' = 'Warming',
    'countryCZSK'        = 'Country (CZE-SVK)'
  )
  
  clean_param <- function(p) {
    p <- str_remove(p, '^b_')
    p <- str_remove(p, '^(mu_|phi_|zoi_|coi_)')   # ZOIB sub-model prefix
    p <- str_remove_all(p, 'scale')
    sides <- str_split(p, ':')[[1]]
    sides <- ifelse(sides %in% names(term_lut), term_lut[sides], sides)
    # for interactions involving timespan, always put timespan first
    if (length(sides) > 1 && 'timespan' %in% sides)
      sides <- c('timespan', sides[sides != 'timespan'])
    paste(sides, collapse = ' × ')
  }
  
  summ <- summ %>%
    mutate(label = vapply(parameter, clean_param, character(1)))
  
  if (!is.null(rename_labels))
    summ <- summ %>%
    mutate(label = ifelse(label %in% names(rename_labels),
                          rename_labels[label], label))
  
  # Label Ordering
  canonical <- c(
    unname(term_lut),
    'Timespan × ΔCore Area',
    'Timespan × ΔNo. Patches',
    'Timespan × Protected site',
    'Timespan × Warming'
  )
  
  if (order_by_effect) {
    summ <- summ %>% arrange(median)
  } else {
    summ <- summ %>%
      mutate(
        sort_key = match(label, canonical),
        # labels absent from canonical sink to the bottom
        sort_key = ifelse(is.na(sort_key), length(canonical) + row_number(), sort_key)
      ) %>%
      arrange(desc(sort_key)) %>%   # desc -> ggplot shows lowest sort_key at top
      select(-sort_key)
  }
  summ <- summ %>% mutate(label = factor(label, levels = unique(label)))
  
  # Plot
  cols <- c(
    'strong (pd ≥ 0.95)'   = '#2166AC',
    'moderate (pd ≥ 0.90)' = '#74ADD1',
    'weak (pd < 0.90)'          = '#BABABA'
  )
  
  ggplot(summ, aes(y = label, color = evidence)) +
    geom_vline(xintercept = ref_line, linetype = 'dashed',
               color = 'grey50', linewidth = 0.4) +
    # Wide interval (thin line)
    geom_linerange(aes(xmin = lo_wide, xmax = hi_wide), linewidth = 0.7) +
    # Narrow interval (thick line)
    geom_linerange(aes(xmin = lo_narrow, xmax = hi_narrow), linewidth = 2) +
    # Median point with white fill so it sits on top of intervals
    geom_point(aes(x = median), size = point_size,
               shape = 21, fill = 'white', stroke = 1.2) +
    scale_color_manual(
      values = cols,
      name   = 'Evidence',
      guide  = guide_legend(override.aes = list(linewidth = 1.5, size = 3))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 10),
      legend.position    = 'bottom',
      plot.title         = element_text(face = 'bold')
    ) +
    labs(
      x     = paste0('Posterior median [',
                     widths[1] * 100, '% and ',
                     widths[2] * 100, '% HDI]'),
      y     = NULL,
      title = title
    )
}

# conditional Effects Visualization Helper
plot_brms_condeffect <- function(model,
                                 effect,
                                 dat,
                                 x_var,
                                 group_var = NULL,
                                 x_label = NULL,
                                 y_label = NULL,
                                 group_label = NULL,
                                 group_quantiles = c(0.10, 0.50, 0.90),
                                 title = NULL,
                                 subtitle = NULL,
                                 new_pal = NULL,
                                 alpha.point = .5, size.point = 1.5,
                                 add_bg_points = TRUE) {
  
  x_nm <- x_var
  grp_nm <- if (!is.null(substitute(group_var))) group_var else NULL
  
  int_cond          <- NULL
  grp_is_continuous <- !is.null(grp_nm) && is.numeric(dat[[grp_nm]])
  
  if (grp_is_continuous) {
    q_vals   <- quantile(dat[[grp_nm]], group_quantiles, na.rm = TRUE)
    x_seq    <- seq(range(dat[[x_nm]], na.rm = TRUE)[1],
                    range(dat[[x_nm]], na.rm = TRUE)[2],
                    length.out = 200)
    int_cond           <- list()
    int_cond[[grp_nm]] <- as.numeric(q_vals)
    int_cond[[x_nm]]   <- x_seq
  }
  
  ce    <- conditional_effects(model, effects = effect,
                               int_conditions = int_cond,
                               plot = FALSE)
  ce_df <- as.data.frame(ce[[effect]])
  
  # convert continuous moderator to labelled factor
  if (grp_is_continuous) {
    q_vals   <- sort(unique(round(ce_df[[grp_nm]], 8)))
    q_labels <- paste0(c("Low", "Medium", "High")[seq_along(q_vals)],
                       " (", round(q_vals, 2), ")")
    ce_df[[grp_nm]] <- factor(round(ce_df[[grp_nm]], 8),
                              levels = q_vals, labels = q_labels)
  }
  
  # colour palette
  n_grp <- if (!is.null(grp_nm)) nlevels(factor(ce_df[[grp_nm]])) else 0
  if (is.null(new_pal)) {
    # pal <- c("#D73027", "#4DAC26", "#4575B4")[seq_len(max(n_grp, 1))]
    pal <- c("orange", "royalblue", "green")[seq_len(max(n_grp, 1))]
  } else {
    pal <- new_pal
  }
  
  # base plot
  p <- ggplot(ce_df, aes(x = .data[[x_nm]], y = estimate__))
  
  # background points
  if (add_bg_points) {
    y_nm  <- as.character(as.formula(model$formula$formula)[[2]])
    bg_df <- dat
    
    if (grp_is_continuous) {
      # assign each point to nearest quantile group label for colours (assign to highest closeness to quantile group)
      q_vals_raw      <- quantile(dat[[grp_nm]], group_quantiles, na.rm = TRUE)
      bg_df[[grp_nm]] <- factor(
        q_labels[max.col(-abs(outer(dat[[grp_nm]], q_vals_raw, "-")))],
        levels = q_labels)
      
      p <- p +
        geom_point(data = bg_df, aes(x = .data[[x_nm]], y = .data[[y_nm]],
                                     colour = .data[[grp_nm]]),
                   alpha = alpha.point, size = size.point, inherit.aes = FALSE)
    } else if (!is.null(grp_nm)) {
      p <- p +
        geom_point(data = bg_df, aes(x = .data[[x_nm]], y = .data[[y_nm]],
                                     colour = .data[[grp_nm]]),
                   alpha = alpha.point, size = size.point, inherit.aes = FALSE)
    } else {
      p <- p +
        geom_point(data = bg_df, aes(x = .data[[x_nm]], y = .data[[y_nm]]),
                   colour = "grey50", alpha = alpha.point, size = size.point, inherit.aes = FALSE)
    }
  }
  
  if (!is.null(grp_nm)) {
    p <- p +
      aes(color = .data[[grp_nm]], fill = .data[[grp_nm]]) +
      geom_ribbon(data = ce_df, aes(ymin = lower__, ymax = upper__),
                  alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = pal,
                         name = if (!is.null(group_label)) group_label else grp_nm) +
      scale_fill_manual(values = pal, guide = "none")
  } else {
    p <- p +
      geom_ribbon(aes(ymin = lower__, ymax = upper__),
                  alpha = 0.15, fill = "steelblue", colour = NA) +
      geom_line(linewidth = 1.1, color = "steelblue")
  }
  
  p  +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom") +
    labs(x     = if (!is.null(x_label)) x_label else x_nm,
         y     = y_label,
         title = title,
         subtitle = subtitle)
}


# get response variable name
response_var <- function(model) {
  as.character(formula(model)$formula[[2]])
}

# get conditional effect line 
ce_line <- function(model, var) {
  conditional_effects(model, effects = var, plot = FALSE)[[var]]
}

# get 95% confidence interval
mean_ci <- function(x) {
  x <- na.omit(x)
  m <- mean(x)
  se <- sd(x) / sqrt(length(x))
  ci <- 1.96 * se
  tibble(mean = m, lo = m - ci, hi = m + ci, sd = sd(x))
}

# get file name to save
brms_name <- function(model_name, where = pth2models){
  file.path(pth2models, model_name)
}
