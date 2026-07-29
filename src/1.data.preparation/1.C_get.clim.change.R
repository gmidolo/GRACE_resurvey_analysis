##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
##################################################################################

# Description: Download CHELSA v2.1 monthly temperature (tas, 1980–2018) for a
#              Central European bounding box, compute annual growing-season means
#              (Apr–Sep), and derive a per-cell warming trend raster (°C/year)

##################################################################################


#### 1. Configuration ####


# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
})

# paths
out_dir <- 'data/climate'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
cache_clim_dir <- 'data/climate/chelsa_annual_crop'
if (!dir.exists(cache_clim_dir)) {dir.create(cache_clim_dir)}
chelsa_tas_url <- function(yr, mo) {
  sprintf(
    '/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/tas/%d/CHELSA_tas_%s_%d_V.2.1.tif',
    yr, mo, yr
  )
}

# study area
bbox      <- ext(13.33, 17.57, 47.90, 49.50)
years     <- 1979:2021 # CHELSA v2.1 monthly covers 1979–2018
gs_months <- sprintf('%02d', 4:9) # Apr–Sep growing season


#### 2. Download + cache annual growing-season mean per year ####

ann_paths <- setNames(
  file.path(cache_clim_dir, sprintf('tas_gs_%d.tif', years)),
  years
)

for (yr in years) {
  fp <- ann_paths[[as.character(yr)]]
  if (file.exists(fp)) { cat('.'); next }
  cat('\nyr - ', yr)
  monthly <- map(gs_months, function(mo) {
    crop(rast(chelsa_tas_url(yr, mo)), bbox) - 273.15   # convert to celsius degree
  })
  writeRaster(mean(rast(monthly)), fp, overwrite = TRUE)
}
cat('\n')

# stack results
ann_stack <- rast(ann_paths)
names(ann_stack) <- as.character(years)


#### 3. Grid cell linear trend (°C / year) ####

trend_r <- app(ann_stack, function(v) {
  if (anyNA(v)) return(NA_real_)
  n   <- length(v)
  x_c <- seq_len(n) - (n + 1L) / 2
  coef(lm.fit(cbind(1, x_c), v))[2L]
}, cores = 4)

names(trend_r) <- 'trend_degC_per_year'

#### 4. Export ####
writeRaster(ann_stack, file.path(out_dir, 'chelsa_tas_gs_annual_1980_2018.tif'),
            overwrite = T)
writeRaster(trend_r,   file.path(out_dir, 'chelsa_tas_trend_1980_2018.tif'),
            overwrite = T)

#### 5. Visualize ####
out_dir <- 'data/climate'
ann_stack <- file.path(out_dir, 'chelsa_tas_gs_annual_1980_2018.tif') %>% rast()
trend_r <- file.path(out_dir, 'chelsa_tas_trend_1980_2018.tif') %>% rast()
plot(trend_r, main = 'Growing-season T trend (°C/year), 1979-2021')
plot(mean(ann_stack), main = 'Mean growing-season temperature (°C), 1979-2021')
hist(values(trend_r), breaks = 40,
     main = 'Distribution of per-pixel warming trends',
     xlab = '°C / year')


#### 6. Plotting - country-level patters (optional) ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(viridis)
  library(cowplot)
})

# 6.1 Regions shape #
countries <- ne_countries(scale = 50, returnclass = 'sf')

at <- countries %>% filter(iso_a3 == 'AUT')
cz <- countries %>% filter(iso_a3 == 'CZE')
sk <- countries %>% filter(iso_a3 == 'SVK')
de <- countries %>% filter(iso_a3 == 'DEU')

czsk <- st_union(cz, sk)

crs_metric <- 3035
at_proj <- st_transform(st_geometry(at), crs_metric)
cz_proj <- st_transform(st_geometry(cz), crs_metric)
sk_proj <- st_transform(st_geometry(sk), crs_metric)

at_proj_buff <- st_buffer(at_proj, dist = 100)
at_cz_border_proj <- st_intersection(at_proj_buff, cz_proj)
at_sk_border_proj <- st_intersection(at_proj_buff, sk_proj)
buffer_50km_proj <- st_buffer(at_cz_border_proj, dist = 50000)

crs_wgs84 <- 4326
buffer_50km <- st_transform(buffer_50km_proj, crs_wgs84)

buffer_vect <- vect(buffer_50km)

trend_buff <- trend_r %>% crop(buffer_vect) %>% mask(buffer_vect)
ann_buff   <- ann_stack %>% crop(buffer_vect) %>% mask(buffer_vect)

# 6.2 Regions rasterization #
de_r   <- rasterize(vect(de),   trend_buff, field = 0, background = NA)
at_r   <- rasterize(vect(at),   trend_buff, field = 1, background = NA)
czsk_r <- rasterize(vect(czsk), trend_buff, field = 2, background = NA)

region_r <- cover(de_r, cover(at_r, czsk_r))
region_r[region_r == 0] <- NA  # drop DE / unclassified cells
names(region_r) <- 'region_code'

region_labels <- c('1' = 'AT', '2' = 'CZSK')

trend_buff <- mask(trend_buff, region_r)

# 6.3 Long format annual data #
df_clim_wide <- c(ann_buff, region_r) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  as_tibble() %>%
  mutate(cell_id = row_number(),
         region  = region_labels[as.character(region_code)]) %>%
  filter(!is.na(region)) %>%
  select(-region_code)

year_cols <- setdiff(names(df_clim_wide), c('x', 'y', 'cell_id', 'region'))

df_clim_long <- df_clim_wide %>%
  pivot_longer(
    cols      = all_of(year_cols),
    names_to  = 'year',
    values_to = 'temp_gs'
  ) %>%
  mutate(year = as.integer(year)) %>%
  select(cell_id, x, y, region, year, temp_gs)


# 6.4 Time series: mean +/- SE growing-season temperature  #
df_clim_ts <- df_clim_long %>%
  group_by(year, region) %>%
  summarize(
    mean_temp = mean(temp_gs, na.rm = TRUE),
    sd_temp   = sd(temp_gs, na.rm = TRUE) ,
    se_temp   = sd(temp_gs, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p_ts <- ggplot(df_clim_ts, aes(x = year, y = mean_temp, color = region, fill = region)) +
  geom_ribbon(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  #scale_color_manual(values = c('AT' = '#d95f02', 'CZSK' = '#1b9e77')) +
 # scale_fill_manual(values = c('AT' = '#d95f02', 'CZSK' = '#1b9e77')) +
  theme_minimal() +
  labs(
    title = 'Growing-Season Mean Temperature (1979-2021)',
    subtitle = 'Mean +/- Std. Dev. across 1km cells',
    x = 'Year',
    y = 'Growing-season mean T (°C)',
    color = 'Region', fill = 'Region'
  )
print(p_ts)

# 6.5 Spatial map:  Warming trend across the buffer #
df_trend_map <- trend_buff %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  as_tibble()

at_trimmed   <- st_crop(at, buffer_50km)
czsk_trimmed <- st_crop(czsk, buffer_50km)

# add veg. plots data to the map
dat_sf <- 'data/raw_veg_data/EVA&GRACE_header_20260202.csv' %>%
  read_csv(show_col_types = F)  %>%
  filter(resurvey == 1) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84') 

p_map <- ggplot() +
  geom_sf(data = at_trimmed, fill = 'gray95', color = 'gray70') +
  geom_sf(data = czsk_trimmed, fill = 'gray95', color = 'gray70') +
  
  geom_raster(data = df_trend_map, aes(x = x, y = y, fill = trend_degC_per_year)) +
  
  geom_sf(data = buffer_50km, fill = NA, color = 'black',
          linetype = 'dashed', linewidth = 0.8) +
  
  geom_sf(data = st_transform(at_cz_border_proj, 4326),
          color = 'gray70', linewidth = 1) +
  geom_sf(data = st_transform(at_sk_border_proj, 4326),
          color = 'gray70', linewidth = 1) +
  
  geom_sf(data=dat_sf, col='brown3', alpha=.8) + 

  
  scale_fill_viridis_c(option = 'viridis', name = 'Trend\n(°C/year)') +
  theme_minimal() +
  labs(
    title = 'Spatial patterns - Growing-Season Warming',
    subtitle = '50km buffer of the AT-CZ border, CHELSA 1km, 1979-2021',
    caption = 'Red points: Vegetation plots; Gray line: AT-CZSK Border; Dashed line: 50km Buffer'
  ) +
  coord_sf(
    xlim = st_bbox(buffer_50km)[c('xmin', 'xmax')],
    ylim = st_bbox(buffer_50km)[c('ymin', 'ymax')],
    expand = TRUE
  ) +
  theme(axis.title = element_blank())

print(p_map)

# 6.6 Combine and export #
fplot <- plot_grid(
  p_map, p_ts, ncol=1, labels = LETTERS[1:2]
)
fplot

ggsave('fig/SI_warming_trend.png', width = 6, height = 8, dpi = 600, bg='white')