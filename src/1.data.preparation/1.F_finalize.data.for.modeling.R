##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
##################################################################################

# Description: Finalize data for modeling and analyses

##################################################################################

#### 0. Package & Header Data Setup ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(betapart)
  library(geosphere)
  library(readxl)
})

# where predictor data is stored
pth2predictors <- 'data/predictors/'

# where to export final data
pt2modeldata <- 'data/model_data/'
if (!dir.exists(pt2modeldata)) {dir.create(pt2modeldata)}

# load header data
hea <- 'data/raw_veg_data/EVA&GRACE_header_20260202.csv' %>%
  read_csv(show_col_types = FALSE)

# assign proper country subdivision
hea$country <- ifelse(hea$EVA_country == 'Austria', 'AT', 'CZSK')

# load species data
spd <- 'data/raw_veg_data/EVA&GRACE_species_20260202.csv' %>%
  read_csv(show_col_types = FALSE) %>%
  semi_join(hea, by = 'plot_id') %>%
  filter(taxongroup == 'Vascular plant')

#### 1. Species Richness Dataset ####

# calculate species richness change and join metadata
dat <- spd %>%
  select(resurvey, plot_id, species) %>%
  distinct() %>%
  group_by(resurvey, plot_id) %>%
  summarise(S = n(), .groups = 'drop') %>%
  left_join(hea, by = c('resurvey', 'plot_id')) %>%
  mutate(
    resurvey = ifelse(resurvey == 1, 'resurvey', 'historic'),
    sampling_year = as.numeric(str_sub(recording_date, -4, -1))
  ) %>%
  select(plot_id, resurvey, S, lon, lat, EVA_dataset, country, plot_size, sampling_year, surveyor_name) %>%
  rename_with(~ str_remove_all(.x, 'EVA_'))

# add surveyor IDs
surveyor_name_ids <- 'data/raw_veg_data/surveyor_name_IDs.csv' %>% 
  read_csv(show_col_types = FALSE) %>% 
  select(1, 3) %>% 
  distinct() %>%
  rename(surveyor_name = name)

dat <- dat %>% 
  left_join(surveyor_name_ids, by = 'surveyor_name') %>%
  arrange(plot_id)

# add environmental predictors
prdcrts <- pth2predictors %>%
  list.files(full.names = TRUE) %>%
  map(read_csv, show_col_types = FALSE) %>%
  map(function(x) {
    if ('resurvey' %in% colnames(x)) {
      x$resurvey <- ifelse(x$resurvey == 1, 'resurvey', 'historic')
    }
    return(x)
  })

dat <- dat %>% 
  left_join(prdcrts[[1]], by = c('plot_id', 'resurvey')) %>%
  left_join(
    prdcrts[[2]] %>% 
      select(plot_id, contains('grassland')) %>% 
      select(-contains('_mesh_'), -contains('_ed_')), 
    by = 'plot_id'
  ) %>%
  left_join(prdcrts[[3]], by = c('plot_id', 'resurvey')) %>%
  rename(survey_type = resurvey)


#### 2. Calculate other metrics too ####

# calculate Jaccard and Bray-Curtis dissimilarity
plot_ids <- unique(spd$plot_id)
res_list <- vector('list', length(plot_ids))

for (i in seq_along(plot_ids)) {
  pid <- plot_ids[i]
  
  x <- spd %>% 
    filter(plot_id == pid) %>%
    filter(!str_detect(species, ' species')) # remove genus-only records
  
  # Bray–Curtis (cover)
  bray <- x %>%
    group_by(resurvey, species) %>%
    summarise(cover = sum(cover), .groups = 'drop') %>%
    pivot_wider(names_from = species, values_from = cover, values_fill = 0) %>%
    select(-resurvey) %>%
    vegdist(method = 'bray') %>%
    as.numeric()
  
  # Jaccard (presence–absence)
  pa_matrix <- x %>%
    select(resurvey, species) %>%
    distinct() %>%
    mutate(pa = 1) %>%
    pivot_wider(names_from = species, values_from = pa, values_fill = 0)
  
  pa_numeric <- pa_matrix %>% select(-resurvey)
  
  jaccard <- vegdist(pa_numeric, method = 'jaccard', binary = TRUE) %>%
    as.numeric()
  
  # betapart decomposition: Jaccard family
  bp <- beta.pair(pa_numeric, index.family = 'jaccard')
  
  res_list[[i]] <- tibble(
    plot_id = pid,
    jaccard = round(jaccard, 4),
    bray.curtis = round(bray, 4),
    beta.jtu = round(as.numeric(bp$beta.jtu), 4),  # turnover
    beta.jne = round(as.numeric(bp$beta.jne), 4)   # nestedness
  )
}

diss_res <- bind_rows(res_list)

# split header metadata for contrasts
hist_hea <- hea %>% 
  filter(resurvey == 0) %>% 
  rename(historic_lon = lon, historic_lat = lat, historic_plot_size = plot_size, historic_recording_date = recording_date, historic_surveyor_name = surveyor_name)

rsrv_hea <- hea %>% 
  filter(resurvey == 1) %>% 
  rename(rsrv_lon = lon, rsrv_lat = lat, rsrv_plot_size = plot_size, rsrv_recording_date = recording_date, rsrv_surveyor_name = surveyor_name)

# extract historic & resurvey metrics directly from memory
hist_S <- dat %>% 
  filter(survey_type == 'historic') %>% 
  select(plot_id, S, sampling_year) %>%
  left_join(hist_hea, by = 'plot_id') %>%
  rename(historic_S = S, historic_sampling_year = sampling_year) %>% 
  select(-contains('recording_date'), -resurvey)

rsrv_S <- dat %>% 
  filter(survey_type == 'resurvey') %>% 
  select(plot_id, S, sampling_year) %>%
  left_join(rsrv_hea, by = 'plot_id') %>%
  rename(rsrv_S = S, rsrv_sampling_year = sampling_year) %>% 
  select(-contains('recording_date'), -contains('EVA_'), -country, -resurvey)

rsrv_env <- dat %>% 
  filter(survey_type == 'resurvey') %>% # use data at resurvey plot coordinates
  select(plot_id, contains('elevation'), contains('gs_temp_change_'), contains('bio'), contains('grassland_'), contains('pa_'))

# merge temporal metrics
final_merge <- rsrv_S %>% 
  left_join(hist_S, by = 'plot_id') %>% 
  left_join(rsrv_env, by = 'plot_id') %>% 
  left_join(diss_res, by = 'plot_id') %>% 
  mutate(logRR_S = round(log10(rsrv_S / historic_S), 4)) %>%
  select(
    plot_id, rsrv_S, historic_S, logRR_S, jaccard, beta.jtu, beta.jne, country, 
    contains('_lon'), contains('_lat'), contains('_sampling_year'), contains('_plot_size'), 
    contains('_surveyor_name'), everything()
  ) %>%
  select(-contains('bio')) %>%
  mutate(
    diff_plot_size = rsrv_plot_size - historic_plot_size,
    timespan = rsrv_sampling_year - historic_sampling_year
  ) %>%
  relocate(diff_plot_size, timespan, .before = rsrv_plot_size)

# calculate geographic distances
final_merge <- final_merge %>%
  mutate(
    dist_meters = distHaversine(
      matrix(c(historic_lon, historic_lat), ncol = 2),
      matrix(c(rsrv_lon, rsrv_lat), ncol = 2)
    ),
    dist_meters = round(dist_meters, 2),
    dist_meters = ifelse(dist_meters == 0, 0.1, dist_meters)
  ) %>%
  relocate(dist_meters, .after = historic_lat)

#### 3. Community Mean EIVs ####

# load & format country-specific EIVs
tichyEIVs_raw <- 'data/raw_veg_data/Indicator.values-tables-2022-11-07-Zenodo.v2.xlsx' %>%
  read_excel(sheet = 2) %>%
  select(Source, Taxon, L:N) 

tichyEIVs <- tichyEIVs_raw %>%
  drop_na() %>%
  semi_join(data.frame(Source = c('Czech Republic', 'Austria')), by = 'Source') %>% 
  group_by(Taxon, Source) %>% 
  summarise(across(everything(), \(x) median(x, na.rm = TRUE)), .groups = 'drop_last') %>%
  ungroup() %>% 
  select(-Source) %>%
  group_by(Taxon) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
  rename(species = Taxon)

# species synonym mapping
species_map <- tibble::tribble(
  ~species,                         ~new_species,
  'Cerastium fontanum aggr.',      'Cerastium fontanum subsp. vulgare',
  'Festuca stricta',               'Festuca stricta subsp. sulcata',
  'Festuca stricta',               'Festuca stricta subsp. stricta',
  'Festuca stricta',               'Festuca stricta subsp. trachyphylla',
  'Euphrasia rostkoviana',         'Euphrasia rostkoviana aggr.',
  'Asplenium adiantum-nigrum',     'Asplenium adiantum-nigrum subsp. serpentini',
  'Salix cinerea',                 'Salix cinerea subsp. cinerea',
  'Asparagus officinalis',         'Asparagus officinalis aggr.',
  'Juniperus communis',            'Juniperus communis subsp. nana',
  'Leopoldia comosa',              'Muscari comosum',
  'Leopoldia tenuiflora',          'Muscari tenuiflorum',
  'Ornithogalum orthophyllum',     'Ornithogalum kochii',
  'Senecio squalidus',             'Senecio rupestris',
  'Helianthemum nummularium',      'Helianthemum nummularium subsp. obscurum',
  'Pulsatilla vulgaris',           'Pulsatilla vulgaris subsp. Grandis',
  'Armeria maritima',              'Armeria maritima subsp. elongata'
)

corrMiss <- tichyEIVs %>%
  inner_join(species_map, by = 'species') %>%
  mutate(species = new_species) %>%
  select(-new_species)

tichyEIVs <- tichyEIVs %>% bind_rows(corrMiss)

# join EIVs to plot species list
spdEIV <- spd %>%
  filter(!str_detect(species, ' species')) %>%
  select(plot_id, species) %>% 
  distinct() %>%
  count(species, name = 'frq') %>%
  left_join(tichyEIVs, by = 'species') %>%
  distinct()

miss_spdEIV <- spdEIV %>% anti_join(drop_na(spdEIV), by = 'species')

# integrate EIVs from full raw table for missing taxa
tichyEIVs_raw_miss <- tichyEIVs_raw %>% 
  semi_join(miss_spdEIV, by = c('Taxon' = 'species')) %>%
  select(-Source) %>%
  group_by(Taxon) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
  drop_na() %>%
  rename(species = Taxon) %>%
  left_join(spdEIV %>% select(species, frq), by = 'species')

spdEIV <- spdEIV %>% 
  anti_join(tichyEIVs_raw_miss, by = 'species') %>%
  bind_rows(tichyEIVs_raw_miss)

# calculate community mean EIVs
spdCMEIV <- spd %>%
  left_join(spdEIV, by = 'species') %>%
  select(resurvey, plot_id, species, L, `T`, M, R, N) %>%
  distinct() %>%
  group_by(resurvey, plot_id) %>%
  summarise(
    across(c(L, `T`, M, R, N),
           list(
             mean = ~mean(.x, na.rm = TRUE),
             n = ~mean(!is.na(.x))
           ),
           .names = '{fn}.{col}'
    ),
    .groups = 'drop'
  )

# calculate temporal EIV differences per plot
diff_EIV_plot <- spdCMEIV %>%
  pivot_wider(
    id_cols = plot_id,
    names_from = resurvey,
    values_from = starts_with('mean')
  ) %>%
  mutate(
    diff.CM_L = mean.L_1 - mean.L_0,
    diff.CM_T = mean.T_1 - mean.T_0,
    diff.CM_M = mean.M_1 - mean.M_0,
    diff.CM_R = mean.R_1 - mean.R_0,
    diff.CM_N = mean.N_1 - mean.N_0
  ) %>%
  select(plot_id, starts_with('diff.CM'))

# merge EIV differences into contrast dataset
data_final <- diff_EIV_plot %>%
  mutate(across(starts_with('diff.CM'), \(x) round(x, 4))) %>%
  left_join(final_merge, by = 'plot_id') %>%
  select(plot_id, rsrv_S, historic_S, logRR_S, jaccard, beta.jtu, beta.jne, everything()) %>%
  select(-bray.curtis, -beta.jtu, -beta.jne)

# glimpse
glimpse(data_final)

# export final dataset
data_final %>%
  write_csv(paste0(pt2modeldata,'model_input_data.csv'))