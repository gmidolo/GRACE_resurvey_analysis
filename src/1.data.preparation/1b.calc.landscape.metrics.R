##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.01.2026
##################################################################################

# Description: calculate lanscape metrics

##################################################################################

### load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
  library(landscapemetrics)
})
source('src/0.utils.R')

# version name to append to the file name (corresponds to EVA&GRACE_....csv files)
version_name = '20260202'

# get raster names
plot_ids_rasters <- list.files('landscape/buffer_rasters/', pattern = '\\.tif$') %>%
  str_remove_all('.tif')

# calc ls metrics for each landscape
land_metric_res <- list()
# message progress every n buffer processed
cat_to = 50
st = Sys.time()
for (i in plot_ids_rasters) {
  # load raster
  ri <- paste0('landscape/buffer_rasters/',i,'.tif') %>%
    rast()
  
  # get class table to merge with landscape metrics results
  class_table <- list()
  for (clsl in 1:length(names(ri))) {
    cls = names(ri)[clsl]
    class_table[[cls]] <- as.data.frame(levels(ri[[cls]])[[1]])
    names(class_table[[cls]]) <- c('class','class_name')
    class_table[[cls]]$layer <- clsl
    class_table[[cls]]$layer_name <- cls
  }
  class_table <- bind_rows(class_table)
  
  # calculate class level metics
  class_level_metrics <- calculate_lsm(
    landscape = ri,
    what = c('lsm_c_ca',   # total class area (in hectares)
             'lsm_c_mesh', # effective mesh size: if I randomly pick two points in the landscape, how large is the area within which those two points are likely to be connected (i.e. in the same class)?
             'lsm_c_ai',   # aggregation index: equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
             'lsm_c_np',   # number of patches (equal to patch density, but the latter is standardized per area)
             'lsm_c_ed'    # equals ED = 0 if only one patch is present (and the landscape boundary is not included) and increases, without limit, as the landscapes becomes more patchy
             ) 
  ) 
  
  # merge table
  class_level_metrics <- class_level_metrics %>%
    left_join(class_table, by = join_by(layer, class)) %>%
    select(contains('layer'), contains('class'), everything()) %>%
    select(-id)
  
  # store results for each buffer's plot_id
  land_metric_res[[i]] <- class_level_metrics
  
  # show progress
  serie_i <- which(plot_ids_rasters %in% i)
  if(serie_i %% cat_to == 0){
    cat(serie_i,' - ')
    et = Sys.time()-st
    message(paste0(
      serie_i,
      ' out of ',
      length(plot_ids_rasters),
      ' - time elapsed: ',
      paste(round(et, 2)),
      ' ',
      attr(et, 'units')
    ))
  }
}
print(Sys.time() - st) # around 3 minutes

# bind results
land_metric_res <- land_metric_res %>%
  bind_rows(.id = 'plot_id')

# # grassland results
# land_metric_grsld <- land_metric_res %>%
#   filter(class_name == 'settlement') %>%
#   select(-level, -class_name, -class, - layer) %>%
#   unite('metric', layer_name, metric) %>%
#   spread(metric, value) %>%
#   mutate(
#     # absolute difference
#     diff_ca = class_2023_ca - class_1840_ca,
#     diff_mesh = class_2023_mesh - class_1840_mesh,
#     diff_ai = class_2023_ai - class_1840_ai,
#     diff_np = class_2023_np - class_1840_np,
#     diff_pd = class_2023_pd - class_1840_pd,
#     # lnRR
#     lnRR_ca = log(class_2023_ca / class_1840_ca),
#     lnRR_mesh = log(class_2023_mesh / class_1840_mesh),
#     lnRR_ai = log(class_2023_ai / class_1840_ai),
#     lnRR_np = log(class_2023_np / class_1840_np),
#     lnRR_pd = log(class_2023_pd / class_1840_pd)
#     
#   )
# 
# land_metric_grsld %>% 
#   select(plot_id, contains('diff_')) %>%
#   gather('k','v', contains('diff_')) %>%
#   ggplot(aes(v))+
#   geom_histogram(color='black')+
#   facet_wrap(~k, scales = 'free_x')+
#   theme_bw()
# 
# land_metric_grsld %>% 
#   select(plot_id, contains('lnRR_')) %>%
#   gather('k','v', contains('lnRR_')) %>%
#   ggplot(aes(v))+
#   geom_histogram(color='black')+
#   facet_wrap(~k, scales = 'free_x')+
#   theme_bw()
# 
# # buffer centroids
# 
# land_metric_grsld %>% 
#   select(plot_id, contains('lnRR_')) %>%
#   gather('k','v', contains('lnRR_')) %>%
#   left_join(
#     'data/EVA&GRACE_header_20251114.csv' %>%
#       read_csv(show_col_types = F)  %>% 
#       select(plot_id, EVA_country) %>% 
#       mutate(plot_id = as.character(plot_id))
#   ) %>%
#   filter(EVA_country !='Slovak Republic') %>%
#   ggplot(aes(EVA_country, v, fill=EVA_country))+
#   geom_boxplot(color='black', notch = T)+
#   facet_wrap(~k, scales = 'free_y')+
#   theme_bw()
# 
# land_metric_grsld %>% 
#   select(plot_id, contains('class_1840_')) %>%
#   gather('k','v', contains('class_1840_')) %>%
#   left_join(
#     'data/EVA&GRACE_header_20251114.csv' %>%
#       read_csv(show_col_types = F)  %>% 
#       select(plot_id, EVA_country) %>% 
#       mutate(plot_id = as.character(plot_id))
#   ) %>%
#   filter(EVA_country !='Slovak Republic') %>%
#   ggplot(aes(EVA_country, v, fill=EVA_country))+
#   geom_boxplot(color='black', notch = T)+
#   facet_wrap(~k, scales = 'free_y')+
#   theme_bw()
# 
# land_metric_grsld %>% 
#   select(plot_id, contains('class_2023_')) %>%
#   gather('k','v', contains('class_2023_')) %>%
#   left_join(
#     'data/EVA&GRACE_header_20251114.csv' %>%
#       read_csv(show_col_types = F)  %>% 
#       select(plot_id, EVA_country) %>% 
#       mutate(plot_id = as.character(plot_id))
#   ) %>%
#   filter(EVA_country !='Slovak Republic') %>%
#   ggplot(aes(EVA_country, v, fill=EVA_country))+
#   geom_boxplot(color='black', notch = T)+
#   facet_wrap(~k, scales = 'free_y')+
#   theme_bw()
# 
# 
# bc <- bc %>%
#   left_join(land_metric_grsld %>% mutate(plot_id = as.numeric(plot_id)))

# pivot wider
land_metric_wider <- land_metric_res %>%
  select(-level, -class, -layer) %>%
  mutate(class_name = str_replace(class_name, ' ', '.'),
         layer_name = str_remove(layer_name, 'class_')) %>%
  unite('metric', class_name, metric, layer_name) %>%
  spread(metric, value)
  

# export
land_metric_wider %>%
  write_csv(
    paste0('env.predictors/class.level_landscapemetrics_',version_name,'.csv')
    )


# check relationships ####
land_metric_wider <- read_csv(paste0('env.predictors/class.level_landscapemetrics_',version_name,'.csv'), show_col_types = F)
# Calculate landscape metrics change
land_metric_wider$ed_change <- land_metric_wider$grassland_ed_2023 - land_metric_wider$grassland_ed_1870
land_metric_wider$ai_change <- land_metric_wider$grassland_ai_2023 - land_metric_wider$grassland_ai_1870
land_metric_wider$np_change <- land_metric_wider$grassland_np_2023 - land_metric_wider$grassland_np_1870
land_metric_wider$ca_change <- land_metric_wider$grassland_ca_2023 - land_metric_wider$grassland_ca_1870

# change im core area of other cover classes (settlments and arable land)
hist(land_metric_wider$settlement_ca_2023 - land_metric_wider$settlement_ca_1870)
hist(land_metric_wider$arable.field_ca_2023 - land_metric_wider$arable.field_ca_1870)
hist(land_metric_wider$vineyard_ca_2023 - land_metric_wider$vineyard_ca_1870)


par(mfrow=c(2,2))
hist(land_metric_wider$ca_change, main = 'Core area change')
hist(land_metric_wider$ed_change, main = 'Edge density change')
hist(land_metric_wider$ai_change, main = 'Aggregation index change')
hist(land_metric_wider$np_change, main = 'No. patches change')

# Examples of landscapes
par(mfrow=c(3,2))
# large increase in AI
sel_plts <- arrange(land_metric_wider, desc(ai_change)) %>% slice(1:3) %>% select(plot_id, contains('grassland_ai'))
sel_plts
for (i in sel_plts$plot_id) {
  rend <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_2023
  rsrt <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_1870
  plot_consistent_landscape(rend, "1870")
  plot_consistent_landscape(rsrt, "2023")
}
# large decrease in AI
sel_plts <- arrange(land_metric_wider, ai_change) %>% slice(1:3) %>% select(plot_id, contains('grassland_ai'))
sel_plts
for (i in sel_plts$plot_id) {
  rend <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_2023
  rsrt <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_1870
  plot_consistent_landscape(rend, "1870")
  plot_consistent_landscape(rsrt, "2023")
}
# large increase in AI
sel_plts <- arrange(land_metric_wider, desc(ed_change)) %>% slice(1:3) %>% select(plot_id, contains('grassland_ed'))
sel_plts
for (i in sel_plts$plot_id) {
  rend <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_2023
  rsrt <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_1870
  plot_consistent_landscape(rend, "1870")
  plot_consistent_landscape(rsrt, "2023")
}
# large decrease in AI
sel_plts <- arrange(land_metric_wider, ed_change) %>% slice(1:3) %>% select(plot_id, contains('grassland_ed'))
sel_plts
for (i in sel_plts$plot_id) {
  rend <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_2023
  rsrt <- rast(paste0('landscape/buffer_rasters/', i, '.tif'))$class_1870
  plot_consistent_landscape(rend, "1870")
  plot_consistent_landscape(rsrt, "2023")
}


#Explore indices ####
# see correlations
land_metric_wider %>% 
  select(contains('grassland')) %>% 
  select(contains('2023')) %>% 
  cor() %>% 
  corrplot::corrplot.mixed(upper='ellipse')

# --- 1. EDGE DENSITY (ED) ---
par(mfrow=c(3,2))

plt_min_ed <- land_metric_wider %>% filter(grassland_ed_2023 == min(grassland_ed_2023, na.rm=T)) %>% pull(plot_id)
plt_max_ed <- land_metric_wider %>% filter(grassland_ed_2023 == max(grassland_ed_2023, na.rm=T)) %>% pull(plot_id)

r1_ed <- rast(paste0('landscape/buffer_rasters/', plt_min_ed[1], '.tif'))$class_2023
r2_ed <- rast(paste0('landscape/buffer_rasters/', plt_max_ed[1], '.tif'))$class_2023

plot_consistent_landscape(r1_ed, "Min Edge Density")
plot_consistent_landscape(r2_ed, "Max Edge Density")


# --- 2. AGGREGATION INDEX (AI) ---

plt_min_ai <- land_metric_wider %>% filter(grassland_ai_2023 == min(grassland_ai_2023, na.rm=T)) %>% pull(plot_id)
plt_max_ai <- land_metric_wider %>% filter(grassland_ai_2023 == max(grassland_ai_2023, na.rm=T)) %>% pull(plot_id)

r1_ai <- rast(paste0('landscape/buffer_rasters/', plt_min_ai[1], '.tif'))$class_2023
r2_ai <- rast(paste0('landscape/buffer_rasters/', plt_max_ai[1], '.tif'))$class_2023

plot_consistent_landscape(r1_ai, "Min Aggregation Index")
plot_consistent_landscape(r2_ai, "Max Aggregation Index")


# --- 3. PATCH DENSITY (NP) ---

plt_min_np <- land_metric_wider %>% filter(grassland_np_2023 == min(grassland_np_2023, na.rm=T)) %>% sample_n(1) %>% pull(plot_id)
plt_max_np <- land_metric_wider %>% filter(grassland_np_2023 == max(grassland_np_2023, na.rm=T)) %>% sample_n(1) %>% pull(plot_id)

r1_np <- rast(paste0('landscape/buffer_rasters/', plt_min_np, '.tif'))$class_2023
r2_np <- rast(paste0('landscape/buffer_rasters/', plt_max_np, '.tif'))$class_2023

plot_consistent(r1_np, "Min Patch Density")
plot_consistent(r2_np, "Max Patch Density")

# Reset layout
par(mfrow=c(1,1))