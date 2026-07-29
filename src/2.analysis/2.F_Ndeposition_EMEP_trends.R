################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
################################################################################
# Description: For supplementary materials only
# Analyze the EMEP Nitrogen Deposition patterns across AT vs CZ-SK (1990-2024)

# Data source: https://www.emep.int/mscw/mscw_moddata.html
################################################################################

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(rnaturalearth)
  library(tidyverse)
  library(broom)
  library(viridis)
  library(cowplot)
})

# 1. Setup data ####
if(!dir.exists('data/EMEP_data')){dir.create('data/EMEP_data')}

# download country boundaries via rnaturalearth
countries <- ne_countries(scale = 50, returnclass = 'sf')

# use iso_a3
at <- countries %>% filter(iso_a3 == 'AUT')
cz <- countries %>% filter(iso_a3 == 'CZE')
sk <- countries %>% filter(iso_a3 == 'SVK')
de <- countries %>% filter(iso_a3 == 'DEU')

# union CZ and SK into a single region
czsk <- st_union(cz, sk)

# project to metric FIRST to safely handle distance operations
crs_metric <- 3035
at_proj <- st_transform(st_geometry(at), crs_metric)
cz_proj <- st_transform(st_geometry(cz), crs_metric)
sk_proj <- st_transform(st_geometry(sk), crs_metric)

# bridge micro-gaps: Buffer Austria by 100 meters to ensure it overlaps Czechia
at_proj_buff <- st_buffer(at_proj, dist = 100)

# extract the shared border overlap
at_cz_border_proj <- st_intersection(at_proj_buff, cz_proj)
at_sk_border_proj <- st_intersection(at_proj_buff, sk_proj)

# create the true 50 km buffer from the border
buffer_50km_proj <- st_buffer(at_cz_border_proj, dist = 50000)

# transform back to WGS84 (EPSG:4326) to match the EMEP raster data
crs_wgs84 <- 4326
buffer_50km <- st_transform(buffer_50km_proj, crs_wgs84)
print(buffer_50km)
plot(buffer_50km)

# 2. Download EMEP helper function ####
get_emep_url <- function(year) {
  base_url <- 'https://thredds.met.no/thredds/fileServer/data/EMEP/2025_Reporting/'
  if (year <= 2022) {
    file_name <- paste0('EMEP01_rv5.6_year.', year, 'met_', year, 'emis_rep2025.nc')
  } else if (year == 2023) {
    file_name <- 'EMEP01_rv5.6_year.2023met_2023emis.nc'
  } else if (year == 2024) {
    file_name <- 'EMEP01_rv5.6_year.2024met_2023emis.nc'
  }
  return(paste0(base_url, file_name))
}


# 3. Cell-country assignment ####
# we download the 1990 file first to establish the spatial grid and classify cells.
year1 <- 1990
url1 <- get_emep_url(year1)
dest1 <- file.path('data/EMEP_data', basename(url1))

if (!file.exists(dest1)) {
  message('Downloading template year (1990) for spatial setup...')
  download.file(url1, dest1, mode = 'wb', quiet = T)
}

# load the raster template
template_nc <- rast(dest1, subds = 'DDEP_OXN_m2Grid') 

# crop and strictly Mask to the buffer
template_crop <- crop(template_nc, buffer_50km)
template_mask <- mask(template_crop, vect(buffer_50km))

# extract cell centroids
cell_points <- as.points(template_mask, na.rm = T)
cell_points_sf <- st_as_sf(cell_points) %>% 
  mutate(cell_id = 1:n()) # Assign a unique ID to each cell

# classify by Region based on Centroid intersection
points_in_at <- st_intersects(cell_points_sf, at, sparse = F)[, 1]
points_in_czsk <- st_intersects(cell_points_sf, czsk, sparse = F)[, 1]
points_in_de <- st_intersects(cell_points_sf, de, sparse = F)[, 1]

# assign region and filter out cells in Germany or unclassified cells
classified_cells <- cell_points_sf %>% 
  mutate(
    region = case_when(
      points_in_de ~ 'Exclude',
      points_in_at ~ 'AT',
      points_in_czsk ~ 'CZSK',
      TRUE ~ 'Exclude'
    )
  ) %>% 
  filter(region != 'Exclude')

# grid plot assignment
p_sanity <- ggplot() +
  geom_sf(data = at, fill = NA, color = 'gray50') +
  geom_sf(data = czsk, fill = NA, color = 'gray50') +
  geom_sf(data = buffer_50km, fill = NA, color = 'black', linetype = 'dashed', linewidth = 1) +
  geom_sf(data = classified_cells, aes(color = region), size = 1.5, alpha = 0.8) +
  #scale_color_manual(values = c('AT' = '#d95f02', 'CZSK' = '#1b9e77')) +
  coord_sf(
    xlim = st_bbox(buffer_50km)[c('xmin', 'xmax')], 
    ylim = st_bbox(buffer_50km)[c('ymin', 'ymax')],
    expand = F
  ) +
  theme_minimal() +
  labs(title = 'Centroid-Based Grid Classification')

print(p_sanity)


# 4. Download loop ####
years <- 1990:2024
results_list <- list()

# convert sf points to SpatVector for fast terra extraction
cells_vect <- vect(classified_cells)

message('Starting extraction...')
for (yr in years) {
  url <- get_emep_url(yr)
  dest <- file.path('data/EMEP_data', basename(url))
  
  if (!file.exists(dest)) {
    message(paste('Downloading', yr, '...'))
    download.file(url, dest, mode = 'wb', quiet = T)
  }
  
  # load the 4 required subdatasets (Updated for correct 'Grid' names)
  vars <- c('DDEP_OXN_m2Grid', 'DDEP_RDN_m2Grid', 'WDEP_OXN', 'WDEP_RDN')
  r_yr <- rast(dest, subds = vars)
  
  # calculate Total N Deposition (sum layers)
  # EMEP unit: mg N / m2 / yr; convert target: kg N / ha / yr
  total_n_raster <- sum(r_yr, na.rm = T) / 100
  names(total_n_raster) <- 'total_N_dep'
  
  # extract values specifically at our pre-classified centroid locations
  extracted_vals <- terra::extract(total_n_raster, cells_vect)
  
  # bind to results list
  results_list[[as.character(yr)]] <- tibble(
    year = yr,
    cell_id = classified_cells$cell_id,
    region = classified_cells$region,
    total_N_dep = extracted_vals$total_N_dep
  )
}

# combine into single long-format data frame
df_long <- bind_rows(results_list)

# 5. Stats and results ####
# compute cumulative N deposition and slope per cell
df_stats <- df_long %>% 
  group_by(cell_id, region) %>% 
  summarize(
    cum_dep = sum(total_N_dep, na.rm = T),
    # Extract slope from linear model (N_dep ~ year)
    slope = coef(lm(total_N_dep ~ year))[2],
    .groups = 'drop'
  )

# wilcoxon rank-sum tests
test_cum <- wilcox.test(cum_dep ~ region, data = df_stats)
test_slope <- wilcox.test(slope ~ region, data = df_stats)

# regional summary table
summary_tbl <- df_stats %>% 
  group_by(region) %>% 
  summarize(
    cell_count = n(),
    mean_cum_dep = mean(cum_dep),
    median_slope = median(slope)
  )
print(summary_tbl)

# Wilcoxon Test: cumulative Deposition
print(test_cum)

# Wilcoxon Test: trend (reg slope)
print(test_slope)


# 6. Plots ####

# time Series Plot with error bands (Mean +/- SE)
df_ts <- df_long %>% 
  group_by(year, region) %>% 
  summarize(
    mean_dep = mean(total_N_dep, na.rm = T),
    sd_dep = sd(total_N_dep, na.rm = T),
    se_dep = sd(total_N_dep, na.rm = T) / sqrt(n()),
    .groups = 'drop'
  )

p_ts <- ggplot(df_ts, aes(x = year, y = mean_dep, color = region, fill = region)) +
  geom_ribbon(aes(ymin = mean_dep - se_dep, ymax = mean_dep + se_dep), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  #scale_color_manual(values = c('AT' = '#d95f02', 'CZSK' = '#1b9e77')) +
  #scale_fill_manual(values = c('AT' = '#d95f02', 'CZSK' = '#1b9e77')) +
  theme_minimal() +
  labs(
    title = 'Annual Total Nitrogen Deposition (1990-2024)',
    subtitle = 'Mean +/- Std. Err. across ~11km cells',
    x = 'Year',
    y = 'Total N Deposition (kg N / ha / yr)',
    color = 'Region', fill = 'Region'
  )
p_ts

# boxplot of cumulative deposition
p_box <- ggplot(df_stats, aes(x = region, y = cum_dep, fill = region)) +
  geom_boxplot(alpha = 0.8, notch = T) +
  #scale_fill_manual(values = c('AT' = '#d95f02', 'CZSK' = '#1b9e77')) +
  theme_minimal() +
  labs(
    title = 'Cumulative Nitrogen Deposition (1990-2024)',
    x = 'Region',
    y = 'Cumulative N Deposition (kg N / ha)',
    fill = 'Region'
  )
p_box

# Spatial Map: cells colored by cumulative Deposition

res_deg <- res(template_nc)[1] # EMEP grid resolution (typically 0.1 degrees)
map_grid <- classified_cells %>% 
  left_join(df_stats %>% select(cell_id, cum_dep), by = 'cell_id') %>% 
  st_buffer(dist = res_deg / 2, endCapStyle = 'SQUARE') # make points into grid squares
# extract X and Y coordinates from the geometry for geom_tile
map_data <- classified_cells %>%
  left_join(df_stats %>% select(cell_id, cum_dep), by = 'cell_id') %>%
  mutate(
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  )

# get EMEP grid resolution for tile dimensions
res_deg <- res(template_nc)[1]

at_trimmed <- st_crop(at, buffer_50km)
czsk_trimmed <- st_crop(czsk, buffer_50km)

# add veg. plots data to the map
dat_sf <- 'data/raw_veg_data/EVA&GRACE_header_20260202.csv' %>%
  read_csv(show_col_types = F)  %>%
  filter(resurvey == 1) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84') 

p_map <- ggplot() +
  
  geom_sf(data = at_trimmed, fill = 'gray95', color = 'gray70') +
  geom_sf(data = czsk_trimmed, fill = 'gray95', color = 'gray70') +
  
  geom_tile(data = map_data, aes(x = x, y = y, fill = cum_dep), 
            width = res_deg, height = res_deg, color = NA) +
  
  geom_sf(data = buffer_50km, fill = NA, color = 'black', 
          linetype = 'dashed', linewidth = 0.8) +
  
  geom_sf(data = st_transform(at_cz_border_proj, 4326),
          color = 'gray70', linewidth = 1) +
  geom_sf(data = st_transform(at_sk_border_proj, 4326),
          color = 'gray70', linewidth = 1) +
  
  scale_fill_viridis_c(option = 'viridis', name = 'Cum. N Dep.\n(kg N/ha)') +
  geom_sf(data=dat_sf, col='brown3', alpha=.8, size=1.25) + 
  theme_minimal() +
  labs(
    title = 'Spatial Pattern of Cumulative N Deposition*',
    subtitle = '* = sum of reduced and oxidized N (wet and dry)',
    caption = 'Red points: Vegetation plots; Gray line: AT-CZSK Border; Dashed line: 50km Buffer'
  ) +
  coord_sf(
    xlim = st_bbox(buffer_50km)[c('xmin', 'xmax')], 
    ylim = st_bbox(buffer_50km)[c('ymin', 'ymax')],
    expand = T
  )+
  theme(axis.title = element_blank())

print(p_map)

# export plots
fplot <- plot_grid(
  p_map, plot_grid(p_box, p_ts, ncol=2, labels = LETTERS[2:3]), ncol=1, labels = LETTERS[1]
)
fplot
ggsave('fig/SI_Ndep_trends.png', width = 9, height = 6, dpi = 500, bg='white')