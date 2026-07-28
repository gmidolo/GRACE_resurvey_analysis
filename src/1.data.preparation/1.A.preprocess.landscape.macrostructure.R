##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 25.01.2026
##################################################################################

# Description: Process 1.5km buffer landscape data into united polygons and rasters

##################################################################################

# 1. load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
})

getwd() # must be './Dropbox/GRACE_Project/analysis/gabri'

# # 2. load landscape buffer data ####
# land_change_buffers_raw <- 'data/landscape/sites_lc_changes_fin.shp' %>%
#   read_sf()
# 
# # plot example
# land_change_buffers_raw %>%
#   filter(plot_id == 56580) %>%
#   select(class_2023, geometry) %>%
#   plot(key.width = lcm(1))
# 
# # 3. get union of buffers inside landscape class categories and periods ####
# st=Sys.time()
# land_change_buffers <- land_change_buffers_raw %>% 
#   gather('key_class','val_class', contains('class_')) %>%
#   group_by(plot_id, key_class, val_class) %>% 
#   summarize(geometry = st_union(geometry))
# Sys.time()-st # about 4 minutes
# 
# # export
# land_change_buffers %>%
#  write_rds('data/landscape/buffer_union_shape.rds', compress = 'gz') 
# 
# # export total union (buffer shape for each id)
# st=Sys.time()
# land_change_buffers %>%
#   group_by(plot_id) %>%
#   summarise(geometry = st_union(geometry)) %>%
#   write_rds('data/landscape/buffer_total.union_shape.rds', compress = 'gz')
# Sys.time()-st # about 1 minute

# 4. buffer reassignment to plot ids ####
buf_shap <- 'data/landscape/buffer_union_shape.rds' %>%
  read_rds()

buf_shap.tot <- 'data/landscape/buffer_total.union_shape.rds' %>%
  read_rds()

# slightly buffer them (by 10 meters)...
buf_shap.tot <- buf_shap.tot %>% st_buffer(10)

# new ids
buf_ids <- buf_shap.tot %>% 
  st_drop_geometry() %>% 
  mutate(buffer_id = 1:nrow(.))
buf_shap <- buf_shap %>% 
  left_join(buf_ids, by = join_by(buffer_id, plot_id)) %>% 
  select(-plot_id)
buf_shap.tot <- buf_shap.tot %>% 
  left_join(buf_ids, by = join_by(plot_id)) %>% 
  select(-plot_id)

# load resurvey data
plots_sel <- 'data/raw_veg_data/EVA&GRACE_header_20260202.csv' %>%
  read_csv(show_col_types = F)  %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84') %>%
  st_transform(st_crs(buf_shap.tot))

# # # explore missing plots
# plots_miss <- plots_sel %>%
#    anti_join(buf_shap.tot %>% st_drop_geometry(), 'plot_id')
# mapview(buf_shap.tot) +
#    mapview(plots_miss, col.regions = 'red')

# get buffer centroids
centr_shap.tot <- buf_shap.tot %>%
  st_centroid()

# for each point in EVA/Grace to the buffers, get index of nearest centroid buffer
idx_near <- st_nearest_feature(plots_sel, centr_shap.tot) 

# new plots assigned
buf_new_assigned <- buf_shap.tot[idx_near, ] %>%
  mutate(resurvey = plots_sel$resurvey,
         plot_id = plots_sel$plot_id, # store proper EVA plot id
         .before = 1 ) %>%
  select(buffer_id, everything())

# compute distances only for the nearest pairs, and add it to the sf
buf_new_assigned$dist_min <- st_distance(plots_sel, centr_shap.tot[idx_near, ], by_element = TRUE)
buf_new_assigned$dist_min <- as.numeric(buf_new_assigned$dist_min)

# # have a look at plots located to a higher distance of 1 km to the centroid of each buffer
buf_new_assigned_high.dist <- buf_new_assigned %>% filter(dist_min > 1000)
plots_sel %>% semi_join(buf_new_assigned_high.dist %>% st_drop_geometry(),
                        'plot_id') %>%
  mapview(col.regions = 'red') +
  mapview(buf_new_assigned_high.dist, col.regions = 'yellow')

# there are plot pairs of survey-resurvey that falls into different buffers..
buf_new_assigned %>% st_drop_geometry() %>% select(plot_id, buffer_id) %>% unique() %>% group_by(plot_id) %>% filter(n()>1) 

#therefore we will only use coordinates from the past (resurvey == 0)..
buf_shap_new <- buf_shap %>% 
  semi_join(
    buf_new_assigned %>% st_drop_geometry() %>% filter(resurvey==0), by = 'buffer_id'
  ) %>% 
  left_join(
    buf_new_assigned %>% st_drop_geometry() %>% filter(resurvey==0), by = join_by(buffer_id, resurvey, dist_min), relationship = 'many-to-many'
  ) %>%
  select(buffer_id, resurvey, plot_id, dist_min, everything()) 

# overwrite buffer shape file
buf_shap_new %>%
  write_rds('data/landscape/buffer_union_shape.rds', compress = 'gz') 


#5. rasterize landscapes ####

buf_shap_new <- read_rds('data/landscape/buffer_union_shape.rds')

# set up nested list of buffers
buffers_l <- buf_shap_new %>%
  split(.$plot_id) %>%
  map(~ split(.x, .x$key_class)) 

# loop to rasterize each plot_id/key_class combination
st=Sys.time()
raster_res <- list()
for (p in names(buffers_l)) {
  
  for (k in names(buffers_l[[p]])) {
    
    vi <- buffers_l[[p]][[k]] %>% 
      vect()
    
    ri <- rast(ext = ext(vi), crs = crs(vi), resolution = 10)
    
    class_factor <- as.factor(vi$val_class)
    vi$val_id <- as.integer(class_factor)      
    lut <- data.frame(value = 1:nlevels(class_factor),
                      label = levels(class_factor))

    zi <- rasterize(vi, ri, field = 'val_id', background = NA, touches = TRUE)
    levels(zi) <- lut
    names(zi) <- 'val_class'
    raster_res[[p]][[k]] <- zi
  }
  
}
Sys.time()-st # about 1 minutes

# export rasters res
pth.for.rasters <- 'data/landscape/buffer_rasters/'
if (!dir.exists(pth.for.rasters)) {dir.create(pth.for.rasters)}

st=Sys.time()
for (i in names(raster_res)) {
  cat(which(i == names(raster_res)), ' - ')
  rast(raster_res[[i]]) %>%
    writeRaster(
      paste0(pth.for.rasters, i, '.tif'),
      overwrite = T
    )
}
Sys.time()-st # ~ 7 minutes
