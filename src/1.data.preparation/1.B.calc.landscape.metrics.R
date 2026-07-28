##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 26.01.2026
##################################################################################

# Description: Calculate various landscape metrics indices

##################################################################################

### load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
  library(landscapemetrics)
})
# source('src/0.utils.R')

# version name to append to the file name (corresponds to EVA&GRACE_....csv files)
version_name <- '20260202'

# path to save predictors data
pth2predictors <- 'data/predictors/'
if (!dir.exists(pth2predictors)) {dir.create(pth2predictors)}

# 1. Calculate landscape metrics using rasterized landscape buffers ####

# get raster names
plot_ids_rasters <- list.files('data/landscape/buffer_rasters/', pattern = '\\.tif$') %>%
  str_remove_all('.tif')

# calc ls metrics for each landscape
land_metric_res <- list()
# message progress every n buffer processed
cat_to = 50
st = Sys.time()
for (i in plot_ids_rasters) {
  # load raster
  ri <- paste0('data/landscape/buffer_rasters/',i,'.tif') %>%
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
    what = c('lsm_c_ca', # total class area (in hectares)
             'lsm_c_mesh', # effective mesh size: if I randomly pick two points in the landscape, how large is the area within which those two points are likely to be connected (i.e. in the same class)?
             'lsm_c_ai', # aggregation index: equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
             'lsm_c_np', # number of patches (equal to patch density, but the latter is standardized per area)
             'lsm_c_ed' # equals ED = 0 if only one patch is present (and the landscape boundary is not included) and increases, without limit, as the landscapes becomes more patchy
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
    paste0(pth2predictors, 'class.level_landscapemetrics_',version_name,'.csv.gz')
    )


# 2. Visualize results ####
land_metric_wider <- read_csv(paste0(pth2predictors, 'class.level_landscapemetrics_',version_name,'.csv.gz'), show_col_types = F)

# Calculate landscape metrics change
land_metric_wider$ed_change <- land_metric_wider$grassland_ed_2023 - land_metric_wider$grassland_ed_1870
land_metric_wider$ai_change <- land_metric_wider$grassland_ai_2023 - land_metric_wider$grassland_ai_1870
land_metric_wider$np_change <- land_metric_wider$grassland_np_2023 - land_metric_wider$grassland_np_1870
land_metric_wider$ca_change <- land_metric_wider$grassland_ca_2023 - land_metric_wider$grassland_ca_1870

library(corrplot)
cor(land_metric_wider %>% select(contains('_change'))) %>%
  corrplot::corrplot.mixed(upper = 'ellipse')

# change im core area of other cover classes (settlments and arable land)
hist(land_metric_wider$settlement_ca_2023 - land_metric_wider$settlement_ca_1870)
hist(land_metric_wider$arable.field_ca_2023 - land_metric_wider$arable.field_ca_1870)
hist(land_metric_wider$vineyard_ca_2023 - land_metric_wider$vineyard_ca_1870)

par(mfrow=c(2,2))
hist(land_metric_wider$ca_change, main = 'Core area change')
hist(land_metric_wider$ed_change, main = 'Edge density change')
hist(land_metric_wider$ai_change, main = 'Aggregation index change')
hist(land_metric_wider$np_change, main = 'No. patches change')