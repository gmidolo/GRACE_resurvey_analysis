##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
##################################################################################

# Description: Prepare protected area status at sites

##################################################################################

### load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
  library(rnaturalearth)
})


#### 1. Extract data ####

# path to save predictors data
pth2predictors <- 'data/predictors/'
if (!dir.exists(pth2predictors)) {dir.create(pth2predictors)}

# path to save protected data metadata
pth2pas <- 'data/protected_areas/'
if (!dir.exists(pth2pas)) {dir.create(pth2pas)}

# version name to append to the file name (corresponds to EVA&GRACE_....csv files)
version_name = '20260202'

### load resurvey metadata
hea <- paste0('data/raw_veg_data/EVA&GRACE_header_', version_name, '.csv') %>%
  read_csv(show_col_types = F)

# spatial data
hea_sp <- hea %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84', remove = F) %>%
  st_transform(
    'data/landscape/buffer_union_shape.rds' %>%
      read_rds() %>%
      st_crs() # use CRS in shape file (ETRS89-extended / LAEA Europe)
  )

hea_buf <- hea_sp %>%
  mutate(geometry = st_buffer(geometry, 10000))

# load pas
pa.file <- paste0(
  'C:/Users/', Sys.getenv('USERNAME'), '/OneDrive - CZU v Praze/czu/wdpa.emma/data/wetransfer_wdpa_eu_041224_2026-01-20_0840/wdpa_eu_041224/wdpa_eu_041224.shp'
)

pas_crop <- pa.file %>%
  read_sf() %>%
  st_crop(st_transform(hea_buf, crs = st_crs(.)))

plot(select(pas_crop, geometry, IUCN_CA))

### export pas crop layer
pas_crop %>%
  write_rds(paste0(pth2pas, 'protected_areas.rds'), compress = 'gz')

# Subset of stricter PAs
pa.inclusion.check <- read_csv(paste0(pth2pas, 'pa.file_GRACE_edit.csv'), show_col_types = F) %>%
  select(1:3, exclude) %>%
  mutate(exclude = ifelse(is.na(exclude), 0, 1)) %>%
  rename(IUCN_CA = `IUCN CATEGORY`)

# PAs to keep
strict_pa_filter <- pa.inclusion.check %>%
  filter(exclude == 1) %>%
  select(-IUCN_CA) %>%
  unique()

# polygon PAs excluded
pas_crop_ex <- pas_crop %>%
  semi_join(strict_pa_filter, by = c('DESIG', 'DESIG_E'))

# polygon PAs to keep
pas_crop_ss <- pas_crop %>%
  anti_join(strict_pa_filter, by = c('DESIG', 'DESIG_E'))

## extract PA at site (All vs Restricted)
extract_pa_sites <- function(sp_data, pa_subset) {
  res_intersect <- st_intersection(sp_data, pa_subset) %>% 
    select(resurvey, plot_id, geometry, IUCN_CA)
  
  res_na <- sp_data %>% 
    anti_join(st_drop_geometry(res_intersect), c('resurvey', 'plot_id')) %>%
    select(resurvey, plot_id, geometry) %>%
    mutate(IUCN_CA = as.character(NA))
  
  res_intersect %>%
    bind_rows(res_na) %>%
    st_drop_geometry()
}

hea_sp_trans <- hea_sp %>% st_transform(crs = st_crs(pas_crop))

hea_pa <- extract_pa_sites(hea_sp_trans, pas_crop)
table(hea_pa$IUCN_CA)
glimpse(hea_pa)

hea_pa_restricted <- extract_pa_sites(hea_sp_trans, pas_crop_ss)
table(hea_pa_restricted$IUCN_CA)
glimpse(hea_pa_restricted)

### final protected areas
pa_data <- hea_pa %>%
  setNames(c('resurvey', 'plot_id', 'IUCN')) %>%
  mutate(pa_site = ifelse(!is.na(IUCN), 1, 0)) %>%
  left_join(
    hea_pa_restricted %>%
      setNames(c('resurvey', 'plot_id', 'IUCN')) %>%
      mutate(pa_site_restrict = ifelse(!is.na(IUCN), 1, 0)) %>%
      select(-IUCN),
    by = c('resurvey', 'plot_id')
  )
glimpse(pa_data)

# export
pa_data %>%
  write_csv(
    paste0(pth2predictors, 'pa_', version_name, '.csv.gz')
  )


#### 2. Plotting ####

# 1. Base country polygons and border lines
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

at_trimmed   <- st_crop(at, buffer_50km)
czsk_trimmed <- st_crop(czsk, buffer_50km)

pas_crop_ss_tr <- st_transform(pas_crop_ss, crs_wgs84)
pas_crop_tr    <- st_transform(pas_crop, crs_wgs84)

# 2. Group DESIG_E into broader, meaningful categories
categorize_pa <- function(df) {
  df %>%
    st_transform(4326) %>%
    mutate(pa_group = case_when(
      str_detect(DESIG_E, 'National Park|National Nature Reserve|National Nature Monument|Nature Reserve|Nature Monument') ~ 'Strict / Core Protected',
      str_detect(DESIG_E, 'Landscape Protection Area|Nature Park|Biosphere|Protected Landscape') ~ 'Landscape Protection (incl. Zone II)',
      str_detect(DESIG_E, 'Habitats Directive|Birds Directive|Community Importance|Special Protection Area|Ramsar') ~ 'Natura 2000 / International',
      TRUE ~ 'Other / Local Designations'
    ))
}

pas_combined <- bind_rows(
  categorize_pa(pas_crop_ss) %>% mutate(Panel = 'A) Protected Areas (Strict Protection; main analyses)'),
  categorize_pa(pas_crop)    %>% mutate(Panel = 'B) Protected Areas (All Designations; sensitivity)')
)

pa_colors <- c(
  'Strict / Core Protected' = 'forestgreen',
  'Landscape Protection (incl. Zone II)' = 'orange',
  'Natura 2000 / International' = '#41b6c4',
  'Other / Local Designations' = '#c51b8a'
)

# 3. Single faceted plot
fplot_pa <- ggplot() +
  geom_sf(data = at_trimmed, fill = 'gray95', color = 'gray70', linewidth = 0.4) +
  geom_sf(data = czsk_trimmed, fill = 'gray95', color = 'gray70', linewidth = 0.4) +
  geom_sf(data = pas_combined, aes(fill = pa_group), color = NA, alpha = 0.6) +
  geom_sf(data = buffer_50km, fill = NA, color = 'black', linetype = 'dashed', linewidth = 0.6) +
  geom_sf(data = st_transform(at_cz_border_proj, crs_wgs84), color = 'gray50', linewidth = 0.8) +
  geom_sf(data = st_transform(at_sk_border_proj, crs_wgs84), color = 'gray50', linewidth = 0.8) +
  scale_fill_manual(name = 'Protection Category', values = pa_colors) +
  facet_wrap(~ Panel, ncol = 1) +
  theme_minimal() +
  coord_sf(
    xlim = st_bbox(buffer_50km)[c('xmin', 'xmax')],
    ylim = st_bbox(buffer_50km)[c('ymin', 'ymax')],
    expand = FALSE, crs = crs_wgs84
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = 'bold'),
    legend.direction = 'vertical',
    strip.text = element_text(size = 9.5, face = 'bold', hjust = 0)
  )
fplot_pa

# 4. Save final figure output
ggsave('fig/SI_protected_areas_comparison.png', plot = fplot_pa, width = 4.5, height = 5.8, dpi = 600, bg = 'white')