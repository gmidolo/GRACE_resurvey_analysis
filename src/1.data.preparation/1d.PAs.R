##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 27.01.2026
##################################################################################

# Description: Prepare protected area at sites

##################################################################################

### load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
})

# version name to append to the file name (corresponds to EVA&GRACE_....csv files)
version_name = '20260202'

### load resurvey metadata
hea <-  paste0('./data/EVA&GRACE_header_',version_name,'.csv') %>%
  read_csv(show_col_types = F)

# spatial data
hea_sp <- hea %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84', remove = F) %>%
  st_transform(
    'landscape/buffer_union_shape.rds' %>%
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

# pas_crop %>% 
#   st_drop_geometry %>%
#   select(DESIG, DESIG_E, IUCN_CA) %>%
#   unique %>%
#   select(IUCN_CA, everything()) %>%
#   arrange(DESIG) %>%
#   write_excel_csv('pa.file_GRACE.csv')

mapview(pas_crop)

# Subset of stricter PAs
pa.inclusion.check <- read_csv('./pa.file_GRACE_edit.csv', show_col_types = F)
pa.inclusion.check <- pa.inclusion.check %>% select(1:3, `Milan exclude`) %>%
  mutate(`Milan exclude`=ifelse(is.na(`Milan exclude`), 0, 1)) %>%
  rename(IUCN_CA=`IUCN CATEGORY`)

pas_crop_ex <- pas_crop %>%
  semi_join(pa.inclusion.check %>% filter(`Milan exclude`==1) %>% select(-IUCN_CA) %>% unique(), by = c('DESIG', 'DESIG_E'))

pas_crop_ss <- pas_crop %>%
  anti_join(pa.inclusion.check %>% filter(`Milan exclude`==1) %>% select(-IUCN_CA) %>% unique(), by = c('DESIG', 'DESIG_E'))

## extract PA at site
hea_sp <- hea_sp %>% st_transform(crs = st_crs(pas_crop))
hea_sp_pa <- st_intersection(hea_sp,pas_crop) %>% 
  select(resurvey, plot_id, geometry, IUCN_CA) 

hea_sp_pa.na <- hea_sp %>% 
  anti_join(st_drop_geometry(hea_sp_pa), c('resurvey', 'plot_id')) %>%
  select(resurvey, plot_id, geometry) %>%
  mutate(IUCN_CA = NA)
hea_sp_pa.na$IUCN_CA <- as.character(hea_sp_pa.na$IUCN_CA) 
hea_sp_pa <- hea_sp_pa %>%
  bind_rows(hea_sp_pa.na)

table(hea_sp_pa$IUCN_CA)

hea_pa <- st_drop_geometry(hea_sp_pa)

glimpse(hea_pa)

## extract PA at site RESTRICTED
hea_sp_pa_restricted <- st_intersection(hea_sp,pas_crop_ss) %>% 
  select(resurvey, plot_id, geometry, IUCN_CA) 

hea_sp_pa_restricted.na <- hea_sp %>% 
  anti_join(st_drop_geometry(hea_sp_pa_restricted), c('resurvey', 'plot_id')) %>%
  select(resurvey, plot_id, geometry) %>%
  mutate(IUCN_CA = NA)
hea_sp_pa_restricted.na$IUCN_CA <- as.character(hea_sp_pa_restricted.na$IUCN_CA) 
hea_sp_pa_restricted <- hea_sp_pa_restricted %>%
  bind_rows(hea_sp_pa_restricted.na)

table(hea_sp_pa_restricted$IUCN_CA)

hea_pa_restricted <- st_drop_geometry(hea_sp_pa_restricted)

glimpse(hea_pa_restricted)

# 
# hea_sp_pa %>% semi_join(data.frame(IUCN_CA = c('Ia', 'II'))) %>% filter(resurvey==1)

## extract area of pas in 1 km buffer

# # rasterize pas
# vi <- vect(pas_crop)
# ri <- rast(ext = ext(vi), crs = crs(vi), resolution = 10)
# class_factor <- as.factor(vi$IUCN_CA)
# vi$val_id <- as.integer(class_factor)      
# lut <- data.frame(value = 1:nlevels(class_factor),
#                   label = levels(class_factor))
# zi <- rasterize(vi, ri, field = 'val_id', background = NA, touches = TRUE)
# levels(zi) <- lut
# names(zi) <- 'val_class'
# 
# # extract areas
# st=Sys.time()
# ext <- extract(
#   zi,
#   hea_sp %>% mutate(geometry = st_buffer(geometry, 1000)),
#   table = TRUE
# )
# Sys.time()-st
# 
# ext <- ext %>%
#   left_join(
#     data.frame(
#       ID = 1:length(unique(ext$ID)),
#       resurvey =  hea_sp$resurvey,
#       plot_id = hea_sp$plot_id
#     ), 'ID'
#   )
# 
# ext_count <- ext %>% 
#   group_by(resurvey, plot_id, val_class) %>%
#   summarise(n = n()) %>%
#   mutate(n = ifelse(is.na(val_class), 0, n)) %>%
#   mutate(area_m2 = n*100) %>%
#   mutate(area_ha = area_m2 /1e4) %>%
#   ungroup()
# 
# ext_count_tot <- ext_count %>%
#   select(-val_class) %>%
#   group_by(resurvey, plot_id) %>%
#   summarise_all(sum) %>%
#   select(-area_m2)
# 
# names(ext_count_tot) <- c('resurvey','plot_id',
#                           'pa_1kmbuff_no.100m2.cells',
#                           'pa_1kmbuff_area_ha')
# 
# ext_count_tot$pa_1kmbuff_area_ha <- round(ext_count_tot$pa_1kmbuff_area_ha, 1)

### final protected areas
# pa_data <- hea_pa %>%
#   setNames(c('resurvey', 'plot_id', 'pa_cat')) %>%
#   mutate(pa_site_iucn = ifelse(!is.na(pa_site_iucn_cat), 1, 0)) %>%
  # left_join(
  #   ext_count_tot, by = c('resurvey','plot_id')
  # )

pa_data <- hea_pa %>%
  setNames(c('resurvey', 'plot_id', 'IUCN')) %>%
  mutate(pa_site = ifelse(!is.na(IUCN), 1, 0)) %>%
  left_join(
    hea_pa_restricted %>%
      setNames(c('resurvey', 'plot_id', 'IUCN'))%>%
      mutate(pa_site_restrict = ifelse(!is.na(IUCN), 1, 0)) %>%
      select(-IUCN)
  )

diff = pa_data %>% 
  filter(pa_site!=pa_site_restrict)

pa_data %>% select(resurvey, plot_id, pa_site) %>%
  spread(resurvey, pa_site) %>%
  filter(`0`!=`1`)

# export
pa_data %>%
  write_csv(
    paste0('env.predictors/pa_', version_name, '.csv')
  )

#### Plotting ####

library(sf)
library(rnaturalearth)

# 1. Base country polygons and border lines
countries <- ne_countries(scale = 50, returnclass = "sf")
at <- countries %>% filter(iso_a3 == "AUT")
cz <- countries %>% filter(iso_a3 == "CZE")
sk <- countries %>% filter(iso_a3 == "SVK")
de <- countries %>% filter(iso_a3 == "DEU")

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
    mutate(pa_group = case_when(
      str_detect(DESIG_E, "National Park|National Nature Reserve|National Nature Monument|Nature Reserve|Nature Monument") ~ "Strict / Core Protected",
      str_detect(DESIG_E, "Landscape Protection Area|Nature Park|Biosphere|Protected Landscape") ~ "Landscape Protection (incl. Zone II)",
      str_detect(DESIG_E, "Habitats Directive|Birds Directive|Community Importance|Special Protection Area|Ramsar") ~ "Natura 2000 / International",
      TRUE ~ "Other / Local Designations"
    ))
}

pas_crop_ss_tr <- categorize_pa(pas_crop_ss_tr)
pas_crop_tr    <- categorize_pa(pas_crop_tr)

# 3. Combine into a single sf dataframe with a panel tag
pas_combined <- bind_rows(
  pas_crop_ss_tr %>% mutate(Panel = "A) Protected Areas (Strict Protection; main analyses)"),
  pas_crop_tr    %>% mutate(Panel = "B) Protected Areas (All Designations; sensitivity)")
)

pa_colors <- c(
  "Strict / Core Protected"        = "forestgreen",
  "Landscape Protection (incl. Zone II)"   = "orange",
  "Natura 2000 / International"    = "#41b6c4",
  "Other / Local Designations"     = "#c51b8a"
)

# 4. Single faceted plot
fplot_pa <- ggplot() +
  geom_sf(data = at_trimmed, fill = "gray95", color = "gray70", linewidth = 0.4) +
  geom_sf(data = czsk_trimmed, fill = "gray95", color = "gray70", linewidth = 0.4) +
  geom_sf(data = pas_combined, aes(fill = pa_group), color = NA, alpha = 0.6) +
  geom_sf(data = buffer_50km, fill = NA, color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_sf(data = st_transform(at_cz_border_proj, crs_wgs84), color = "gray50", linewidth = 0.8) +
  geom_sf(data = st_transform(at_sk_border_proj, crs_wgs84), color = "gray50", linewidth = 0.8) +
  scale_fill_manual(name = "Protection Category", values = pa_colors) +
  facet_wrap(~ Panel, ncol = 1) +
  theme_minimal() +
  coord_sf(
    xlim = st_bbox(buffer_50km)[c("xmin", "xmax")],
    ylim = st_bbox(buffer_50km)[c("ymin", "ymax")],
    expand = FALSE, crs = crs_wgs84
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    legend.direction = "vertical",
    strip.text = element_text(size = 9.5, face = "bold", hjust = 0)
  )

#print(fplot_pa)

# 5. Save final figure output
ggsave('fig/SI_protected_areas_comparison.png', plot = fplot_pa, width = 4.5, height = 5.8, dpi = 600, bg = 'white')
