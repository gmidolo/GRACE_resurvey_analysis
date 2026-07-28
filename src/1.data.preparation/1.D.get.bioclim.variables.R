##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 27.01.2026
##################################################################################

# Description: Put together plot-level bioclim variables and export a csv file

##################################################################################

### load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
  library(corrplot)
})

# path to save predictors data
pth2predictors <- 'data/predictors/'
if (!dir.exists(pth2predictors)) {dir.create(pth2predictors)}

# version name to append to the file name (corresponds to EVA&GRACE_....csv files)
version_name = '20260202'

### load resurvey metadata
hea <-  paste0('data/raw_veg_data/EVA&GRACE_header_',version_name,'.csv') %>%
  read_csv(show_col_types = F)

# spatial data
hea_sp <- hea %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84', remove = F)

#load chelsa raster data
chelsapth <- paste0('C:/Users/', Sys.getenv('USERNAME'), '/OneDrive - CZU v Praze/brno/brno_postdoc/climdata/chelsa')
chelsar <- list.files(chelsapth, full.names = T) %>%
  rast()

# simplify names
names(chelsar) <- names(chelsar) %>% 
  str_remove('_1981-2010_V.2.1') %>% 
  str_remove('CHELSA_')

# extract bioclim
chelsa_dat <- extract(chelsar, hea_sp)
chelsa_dat <- chelsa_dat %>%
  select(-ID) %>%
  mutate_all(round, 1) %>%
  mutate(resurvey = hea_sp$resurvey, plot_id = hea_sp$plot_id, .before = 1)

#load elevation
ele <- paste0('C:/Users/', Sys.getenv('USERNAME'), '/OneDrive - CZU v Praze/brno/brno_postdoc/climdata/topodata/elevation/Copernicus_GLO.90_DEM_Europe.tif') %>%
  rast()
ele_dat <- extract(ele, hea_sp)
ele_dat <- ele_dat %>%
  select(-ID) %>%
  setNames('elevation') %>%
  mutate_all(round, 1) %>%
  mutate(resurvey = hea_sp$resurvey, plot_id = hea_sp$plot_id, .before = 1)

#load temperature change Growing-season T trend (°C/year), 43 years, 1979-2021 (obtained from `1c.get.clim.change.R`)
Tchange <- 'data/climate/chelsa_tas_trend_1980_2018.tif' %>%
  rast()

Tchange_dat <- terra::extract(Tchange, hea_sp)
Tchange_dat$trend_degC_per_year %>% hist
Tchange_dat <- Tchange_dat %>%
  select(-ID) %>%
  setNames('gs_temp_change_decade') %>%
  mutate_all(round, 5) %>%
  mutate(gs_temp_change_decade = gs_temp_change_decade * 10, resurvey = hea_sp$resurvey, plot_id = hea_sp$plot_id, .before = 1)

# biocilm variables
bioclim <- ele_dat %>%
  left_join(Tchange_dat, c('resurvey','plot_id')) %>%
  left_join(chelsa_dat, c('resurvey','plot_id'))

# visualize correlations
corrplot.mixed(cor(bioclim[,-c(1,2)]), upper='ellipse')

# export
bioclim %>%
  write_csv(paste0(pth2predictors, 'bioclim_', version_name, '.csv.gz'))