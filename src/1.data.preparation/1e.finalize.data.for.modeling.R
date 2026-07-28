##################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 27.01.2025
##################################################################################

# Description: Finalize data for modeling and analyses

##################################################################################

####0. header data set up ####

# load packages
library(tidyverse)

# load data
getwd()

# load header data
hea <-  './data/EVA&GRACE_header_20260202.csv' %>%
  read_csv(show_col_types = F)
# hea <-  './data/old/EVA&GRACE_header_20251114.csv' %>%
#   read_csv(show_col_types = F)

# hea$plot_size <- ifelse(hea$plot_size == 125, 12.5, hea$rsrv_plot_size)

# # calculate timespan in years
# hea$timespan <- as.numeric(str_sub(hea$rsrv_recording_date, -4, -1)) - as.numeric(str_sub(hea$EVA_recording_date, -4, -1))

# assing proper country subdivision
hea$country <- ifelse(hea$EVA_country == 'Austria', 'AT', 'CZSK')

# > str(hea)
# spc_tbl_ [822 × 10] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
# $ resurvey      : num [1:822] 0 0 0 0 0 0 0 0 0 0 ...
# $ plot_id       : num [1:822] 1401631 1401630 1401362 1389586 1385622 ...
# $ lon           : num [1:822] 17 17 17 17 17 ...
# $ lat           : num [1:822] 48.2 48.2 48.3 48.2 48.2 ...
# $ EVA_dataset   : chr [1:822] "Slovakia_nvd" "Slovakia_nvd" "Slovakia_nvd" "Slovakia_nvd" ...
# $ EVA_country   : chr [1:822] "Slovak Republic" "Slovak Republic" "Slovak Republic" "Slovak Republic" ...
# $ plot_size     : num [1:822] 25 25 25 30 25 30 10 20 15 16 ...
# $ recording_date: chr [1:822] "03.07.1964" "02.07.1964" "16.07.2004" "08.07.1960" ...
# $ surveyor_name : chr [1:822] "M Kaleta" "M Kaleta" "K Hegedušová" "V Zajacová" ...
# $ country       : chr [1:822] "CZSK" "CZSK" "CZSK" "CZSK" ...
# - attr(*, "spec")=
#   .. cols(
#     ..   resurvey = col_double(),
#     ..   plot_id = col_double(),
#     ..   lon = col_double(),
#     ..   lat = col_double(),
#     ..   EVA_dataset = col_character(),
#     ..   EVA_country = col_character(),
#     ..   plot_size = col_double(),
#     ..   recording_date = col_character(),
#     ..   surveyor_name = col_character()
#     .. )
# - attr(*, "problems")=<externalptr> 
#   

# load species data
spd <- './data/EVA&GRACE_species_20260202.csv' %>%
  read_csv() %>%
  semi_join(hea, 'plot_id') %>%
  filter(taxongroup=='Vascular plant')
# > str(spd)
# spc_tbl_ [23,657 × 6] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
# $ resurvey  : num [1:23657] 0 0 0 0 0 0 0 0 0 0 ...
# $ plot_id   : num [1:23657] 50332 50332 50332 50332 50332 ...
# $ species   : chr [1:23657] "Acer campestre" "Achillea millefolium aggr." "Arrhenatherum elatius" "Astragalus cicer" ...
# $ taxongroup: chr [1:23657] "Vascular plant" "Vascular plant" "Vascular plant" "Vascular plant" ...
# $ layer     : num [1:23657] 6 6 6 6 6 6 6 6 6 6 ...
# $ cover     : num [1:23657] 0.2 1 37 1 1 1 37 3 1 15 ...


####1. species richness dataset #####

# calculate species richness change and join data
dat <- spd %>%
  select(resurvey, plot_id, species) %>%
  unique() %>%
  group_by(resurvey, plot_id) %>%
  summarise(S = n()) %>%
  ungroup() %>%
  left_join(hea)

# resurvey period (as categorical)
dat$resurvey <- ifelse(dat$resurvey == 1, 'resurvey', 'historic')

# # reassign the plot size
# dat$plot_size <- dat$EVA_plot_size
# dat$plot_size <- ifelse(dat$resurvey == 'resurvey', dat$rsrv_plot_size, dat$plot_size)
# 
# # reassign survey years
# dat$sampling_year <- as.numeric(str_sub(dat$EVA_recording_date, -4, -1))
# dat$sampling_year <- ifelse(dat$resurvey == 'resurvey', as.numeric(str_sub(dat$rsrv_recording_date, -4, -1)), dat$sampling_year)

# calc sampling year
dat$sampling_year <- as.numeric(str_sub(dat$recording_date, -4, -1))

# filter variables needed
dat <- dat %>% select(
  plot_id, resurvey, S, lon, lat, EVA_dataset, country, plot_size, sampling_year, surveyor_name, surveyor_name
) 
names(dat) <- str_remove_all(names(dat), 'EVA_')

# # rearrange surveyor names
# dat$surveyor_name <- ifelse(dat$resurvey == 'historic', dat$surveyor_name, dat$rsrv_surveyor_name)
# dat <- dat %>% select(-rsrv_surveyor_name)
# 
# # abbreviate resurveyors names
# dat <- dat %>%
#   mutate(rsrv_surveyor_name = ifelse(str_detect(rsrv_surveyor_name, ' Barbora Chmelov'), 'Petr Keil',  rsrv_surveyor_name )) %>%
#   mutate(resurveyor_id =
#            abbreviate(stringi::stri_trans_general(rsrv_surveyor_name, 'Latin-ASCII'), 4, named = F) %>%
#            str_replace_all('_','') %>% str_replace_all('-','') %>% str_replace_all('\\.','') %>%
#            tolower()) %>%
#   select(-rsrv_surveyor_name)
surveyor_name_ids <- './data/surveyor_name_IDs.csv' %>% 
  read_csv(show_col_types = F) %>% 
  select(1,3) %>% 
  unique() %>%
  rename(surveyor_name=name)

dat <- dat %>% left_join(surveyor_name_ids)

# arrange order
dat <- dat %>%
  arrange(plot_id)

# add environmental predictors
prdcrts <- 'env.predictors' %>%
  list.files(full.names = T) %>%
  map(read_csv, show_col_types = F) %>%
  map(\(x){
    if (any(colnames(x)=='resurvey')){
      y=x
      y$resurvey <- ifelse(y$resurvey == 1, 'resurvey', 'historic')
      return(y)
    } else {
      x
    }
  }
  )

dat <- dat %>% left_join(prdcrts[[1]], by = c('plot_id','resurvey'))
dat <- dat %>% left_join(prdcrts[[2]] %>% 
                           select(plot_id, contains('grassland')) %>% # merge only indeces fro grasslands
                           select(-contains('_mesh_'), -contains('_ed_'))
                         , by = 'plot_id')
dat <- dat %>% left_join(prdcrts[[3]], by = c('plot_id','resurvey'))

dat <- dat %>% rename(survey_type = resurvey)
 
# # export
dat %>%
  write_csv('data/model_data/species.richness_data.csv')


#### 2. Contrast-based dataset #####

# calculate Jaccard and bray curtis
library(vegan)
library(betapart)
plot_ids <- unique(spd$plot_id)

res_list <- list()
for (i in seq_along(plot_ids)) {
  
  pid <- plot_ids[i]
  
  x <- spd %>% 
    filter(plot_id == pid) %>%
    filter(!str_detect(species, ' species')) # remove genus-species only
  
  # Bray–Curtis (cover)
  bray <- x %>%
    group_by(resurvey, species) %>%
    summarise(cover = sum(cover), .groups = "drop") %>%
    spread(species, cover, fill = 0) %>%
    as.data.frame() %>%
    select(-resurvey) %>%
    vegdist(method = "bray") %>%
    as.numeric()
  
  # Jaccard (presence–absence) — build the matrix once, reuse for betapart
  pa_matrix <- x %>%
    select(resurvey, species) %>%
    distinct() %>%
    mutate(pa = 1) %>%
    spread(species, pa, fill = 0) %>%
    as.data.frame()
  
  # store resurvey labels then drop for vegdist/betapart
  pa_numeric <- pa_matrix %>% select(-resurvey)
  
  jaccard <- vegdist(pa_numeric, method = "jaccard", binary = TRUE) %>%
    as.numeric()
  
  # Betapart decomposition: Jaccard family
  # beta.pair returns three dist objects: beta.jtu, beta.jne, beta.jac
  # For a 2-row matrix, each is a single value
  bp <- beta.pair(pa_numeric, index.family = "jaccard")
  
  res_list[[as.character(pid)]] <- tibble(
    plot_id    = pid,
    jaccard    = round(jaccard, 4),
    bray.curtis = round(bray, 4),
    beta.jtu   = round(as.numeric(bp$beta.jtu), 4),  # turnover component
    beta.jne   = round(as.numeric(bp$beta.jne), 4)   # nestedness-resultant component
  )
}
diss_res <- bind_rows(res_list) %>% mutate_all(as.numeric)
diss_res

# split dat 
hist_hea <- hea %>% filter(resurvey==0) %>% rename(historic_lon=lon, historic_lat=lat, historic_plot_size = plot_size, historic_recording_date=recording_date, historic_surveyor_name = surveyor_name)
rsrv_hea <- hea %>% filter(resurvey==1) %>% rename(rsrv_lon=lon, rsrv_lat=lat, rsrv_plot_size = plot_size, rsrv_recording_date=recording_date, rsrv_surveyor_name = surveyor_name)

hist_S <- read_csv('data/model_data/species.richness_data.csv', show_col_types = F) %>% 
  filter(survey_type=='historic') %>% 
  mutate(resurvey = 0, .before = 1) %>% 
  select(resurvey,plot_id,survey_type, S, sampling_year) %>%
  left_join(hist_hea) %>%
  select(-resurvey, -survey_type) %>%
  rename(historic_S = S, historic_sampling_year = sampling_year) %>% 
  select(-contains('recording_date'))

hist_env <- read_csv('data/model_data/species.richness_data.csv', show_col_types = F) %>% 
  filter(survey_type=='historic') %>%
  select(plot_id, contains('elevation'), contains('gs_temp_change_'), contains('bio'), contains('grassland_'), contains('pa_'))
#names(hist_env)[-1] <- paste0('historic_', names(hist_env)[-1])

rsrv_env <- read_csv('data/model_data/species.richness_data.csv', show_col_types = F) %>% 
  filter(survey_type=='resurvey') %>%
  select(plot_id, contains('elevation'), contains('gs_temp_change_'), contains('bio'), contains('grassland_'), contains('pa_'))
#names(rsrv_env)[-1] <- paste0('rsrv_', names(rsrv_env)[-1])

rsrv_S <- read_csv('data/model_data/species.richness_data.csv', show_col_types = F) %>% 
  filter(survey_type=='resurvey') %>% 
  mutate(resurvey = 1, .before = 1) %>% 
  select(resurvey,plot_id,survey_type, S, sampling_year) %>%
  left_join(rsrv_hea) %>%
  select(-resurvey, -survey_type) %>%
  rename(rsrv_S = S, rsrv_sampling_year = sampling_year) %>% 
  select(-contains('recording_date'))
rsrv_S <- rsrv_S %>% select(-contains('EVA_'), - country) 

final_merge <- rsrv_S %>% 
  left_join(hist_S, 'plot_id') %>% 
  # USE historical OR resurvey
  left_join(rsrv_env, 'plot_id') 

final_merge <- final_merge %>%
  left_join(diss_res) %>% # add dissimilarity results
  mutate(logRR_S = round(log10(rsrv_S/historic_S), 4)) 

final_merge <- final_merge %>%
  select(
    plot_id, rsrv_S, historic_S, logRR_S, jaccard, beta.jtu, beta.jne, country, contains('_lon'), contains('_lat'), contains('_sampling_year'), contains('_plot_size'), contains('_surveyor_name'), everything()
  ) %>%
  #remove CHELSA bio
  select(-contains('bio'))

# add plot size difference and timespan
final_merge <- final_merge %>% 
  mutate(diff_plot_size = rsrv_plot_size - historic_plot_size, .before = rsrv_plot_size ) %>% 
  mutate(timespan = rsrv_sampling_year - historic_sampling_year, .before = rsrv_plot_size)

# calculate geographic distances
library(geosphere)
final_merge <- final_merge %>%
  mutate(dist_meters = distHaversine(
    matrix(c(historic_lon, historic_lat), ncol = 2), # Start points
    matrix(c(rsrv_lon, rsrv_lat), ncol = 2)          # End points
  ), .after = historic_lat)

final_merge$dist_meters <- round(final_merge$dist_meters, 2)
# force a minimum distance error (10 cm)
final_merge$dist_meters <- ifelse(final_merge$dist_meters == 0, 0.1, final_merge$dist_meters)

# # export
# final_merge %>%
#   write_csv('data/model_data/contrast.based_data.csv')

#3. CM EIVs####
tichyEIVs <- "data/Indicator.values-tables-2022-11-07-Zenodo.v2.xlsx" %>%
  readxl::read_excel(sheet = 2) %>%
  select(Source, Taxon, L:N) %>%
  drop_na() %>%
  semi_join(
    data.frame(Source=c('Czech Republic', 'Austria'))
  ) %>% 
  group_by(Taxon, Source) %>% # median by country
  summarise_all(
    median, na.rm=T
  ) %>%
  ungroup() %>%  # median across country
  select(-Source) %>%
  group_by(Taxon) %>%
  summarise_all(
    mean, na.rm=T
  ) %>%
  rename(species=Taxon)

species_map <- tibble::tribble(
  ~species,                         ~new_species,
  "Cerastium fontanum aggr.",      "Cerastium fontanum subsp. vulgare",
  "Festuca stricta",               "Festuca stricta subsp. sulcata",
  "Festuca stricta",               "Festuca stricta subsp. stricta",
  "Festuca stricta",               "Festuca stricta subsp. trachyphylla",
  "Euphrasia rostkoviana",         "Euphrasia rostkoviana aggr.",
  "Asplenium adiantum-nigrum",     "Asplenium adiantum-nigrum subsp. serpentini",
  "Salix cinerea",                 "Salix cinerea subsp. cinerea",
  "Asparagus officinalis",         "Asparagus officinalis aggr.",
  "Juniperus communis",            "Juniperus communis subsp. nana",
  'Leopoldia comosa',              "Muscari comosum",
  'Leopoldia tenuiflora',          'Muscari tenuiflorum',
  'Ornithogalum orthophyllum',     'Ornithogalum kochii',
  'Senecio squalidus',             'Senecio rupestris',
  'Helianthemum nummularium',      'Helianthemum nummularium subsp. obscurum',
  'Pulsatilla vulgaris',           'Pulsatilla vulgaris subsp. Grandis',
  'Armeria maritima',              'Armeria maritima subsp. elongata'
)

corrMiss <- tichyEIVs %>%
  inner_join(species_map, by = "species") %>%
  mutate(species = new_species) %>%
  select(-new_species)

tichyEIVs <- tichyEIVs %>% bind_rows(corrMiss)

spdEIV <- spd %>%
  filter(!str_detect(species, ' species')) %>% # remove genus-species only
  select(plot_id, species) %>% unique %>%
  group_by(species) %>%
  summarise(frq=n()) %>%
  left_join(tichyEIVs) %>%
  unique

miss_spdEIV <- spdEIV %>% anti_join(drop_na(spdEIV)) %>% arrange(desc(frq))
print(miss_spdEIV, n=Inf)

# integrate other EIVs from other indicator systems
tichyEIVs_raw <- "data/Indicator.values-tables-2022-11-07-Zenodo.v2.xlsx" %>%
  readxl::read_excel(sheet = 2) %>%
  select(Source, Taxon, L:N) 

tichyEIVs_raw_miss <- tichyEIVs_raw %>% 
  semi_join(data.frame(Taxon=miss_spdEIV$species)) %>%
  select(-Source) %>%
  group_by(Taxon) %>%
  summarise_all(mean, na.rm=T) %>%
  drop_na() %>%
  rename(species=Taxon) %>%
  left_join(spdEIV %>% select(species,frq))

spdEIV <- spdEIV %>% 
  anti_join(tichyEIVs_raw_miss, by ='species') %>%
  bind_rows(tichyEIVs_raw_miss) 

miss_spdEIV <- spdEIV %>%
  anti_join(drop_na(spdEIV)) %>% 
  arrange(desc(frq))
# 
# eive <- paste0("C:/Users/", Sys.info()['user'] ,"/OneDrive - CZU v Praze/czu/intrplEU_EIV/code/data/eiv.data.nomenclature/EIVE.xlsx") %>%
#   readxl::read_excel(sheet=2) %>%
#   select(1, `EIVEres-M`, `EIVEres-N`, `EIVEres-T`,`EIVEres-L`,`EIVEres-R`)
# names(eive) <- c('species','EIVE_M','EIVE_N','EIVE_T','EIVE_L','EIVE_R')
# 
# miss_spdEIV <- miss_spdEIV %>% left_join(eive)
# miss_spdEIV
# #EIV imputation from EIVE
# mL <- lm(L~EIVE_L,data = spdEIV %>% left_join(eive))
# summary(mL)$r.squared
# miss_spdEIV$L <- ceiling(predict(mL, miss_spdEIV) * 2) / 2
# mT <- lm(`T`~EIVE_T,data = spdEIV %>% left_join(eive))
# summary(mT)$r.squared
# miss_spdEIV$`T` <- ceiling(predict(mT, miss_spdEIV) * 2) / 2
# mN <- lm(`N`~EIVE_N,data = spdEIV %>% left_join(eive))
# summary(mN)$r.squared
# miss_spdEIV$N <- ceiling(predict(mN, miss_spdEIV) * 2) / 2
# mM <- lm(`M`~EIVE_M,data = spdEIV %>% left_join(eive))
# summary(mM)$r.squared
# miss_spdEIV$M <- ceiling(predict(mM, miss_spdEIV) * 2) / 2
# mR <- lm(`R`~EIVE_R,data = spdEIV %>% left_join(eive))
# summary(mR)$r.squared
# miss_spdEIV$R <- ceiling(predict(mR, miss_spdEIV) * 2) / 2
# 
# spdEIV_new <- drop_na(spdEIV) %>%
#   bind_rows(miss_spdEIV %>% select(-contains('EIVE_'))) %>%
#   arrange(species)
# 
# tichyEIVs_raw_miss <- tichyEIVs_raw %>% semi_join(data.frame(Taxon=miss_spdEIV$species))
# 
# tichyEIVs_raw_miss <- tichyEIVs_raw_miss %>%
#   select(-Source) %>%
#   group_by(Taxon) %>%
#   summarise_all(mean, na.rm=T) %>%
#   select(Taxon, M, N) %>%
#   drop_na()
# 
# tichyEIVs_raw_miss_new = tichyEIVs_raw_miss %>%
#   rename(species=Taxon, M_tichy=M, N_tichy=N) %>%
#   left_join(miss_spdEIV)
# 
# plot(tichyEIVs_raw_miss_new$M_tichy, tichyEIVs_raw_miss_new$M)
# plot(tichyEIVs_raw_miss_new$N_tichy, tichyEIVs_raw_miss_new$N)

spdCMEIV <- spd %>%
  #left_join(spdEIV_new, by = "species") %>%
  left_join(spdEIV, by = "species") %>%
  select(resurvey, plot_id, species, L, `T`, M, R, N) %>%
  distinct() %>%
  group_by(resurvey, plot_id) %>%
  summarise(
    across(c(L, `T`, M, R, N),
           list(
             mean = ~mean(.x, na.rm = TRUE),
             n = ~mean(!is.na(.x))   # proportion of non-NA
           ),
           .names = "{fn}.{col}"
    ),
    .groups = "drop"
  )
range(spdCMEIV$n.M)
range(spdCMEIV$n.N)#... they are all above equal 70%


diff_EIV_plot <- spdCMEIV %>%
  pivot_wider(
    id_cols = plot_id,
    names_from = resurvey,
    values_from = starts_with("mean")
  ) %>%
  mutate(
    diff.CM_L = mean.L_1 - mean.L_0,
    diff.CM_T = mean.T_1 - mean.T_0,
    diff.CM_M = mean.M_1 - mean.M_0,
    diff.CM_R = mean.R_1 - mean.R_0,
    diff.CM_N = mean.N_1 - mean.N_0
  ) %>%
  select(plot_id, starts_with("diff.CM"))

diff_EIV_plot %>% select(-1) %>% cor() %>% round(2)

final_merge <- 
  diff_EIV_plot %>%
  mutate_all(round, 4) %>%
  left_join(final_merge) %>%
  select(plot_id,rsrv_S,historic_S,logRR_S,jaccard,beta.jtu,beta.jne, everything()) %>%
  select(-bray.curtis)

# export
final_merge %>%
  write_csv('data/model_data/contrast.based_data.csv')

final_merge %>% select(logRR_S:diff.CM_N) %>% GGally::ggpairs()

final_merge %>%
  summarise(
    mean_jtu_share = mean(beta.jtu / jaccard, na.rm = TRUE),
    median_jtu_share = median(beta.jtu / jaccard, na.rm = TRUE)
  )
# On average, 86.5% of Jaccard dissimilarity was attributable to species turnover rather than nestedness-resultant dissimilarity

# #4. Species matrices####
# 
# final_merge <- read_csv('data/model_data/contrast.based_data.csv')
# 
# # Remove genus-level records (e.g., 'Acer species')
# spd_sp <- spd %>%
#   filter(!str_detect(species, ' species'))
# 
# # Long-format presence/absence restricted to plots in the final modeling dataset
# plot_ids_final <- final_merge$plot_id
# 
# hist_pres <- spd_sp %>%
#   filter(resurvey == 0, plot_id %in% plot_ids_final) %>%
#   select(plot_id, species) %>%
#   distinct() %>%
#   mutate(hist = 1L)
# 
# rsrv_pres <- spd_sp %>%
#   filter(resurvey == 1, plot_id %in% plot_ids_final) %>%
#   select(plot_id, species) %>%
#   distinct() %>%
#   mutate(rsrv = 1L)
# 
# # All plot × species combinations across both time points
# full_pa <- hist_pres %>%
#   full_join(rsrv_pres, by = c('plot_id', 'species')) %>%
#   replace_na(list(hist = 0L, rsrv = 0L))
# 
# # Minimum number of plots a species must qualify in to be retained
# min_occ <- 10
# 
# # Species historically present in >= min_occ plots (eligible for extinction)
# ext_species <- hist_pres %>%
#   count(species) %>%
#   filter(n >= min_occ) %>%
#   pull(species)
# 
# # Species that colonized >= min_occ plots (observed colonization events)
# col_species <- full_pa %>%
#   filter(hist == 0L, rsrv == 1L) %>%
#   count(species) %>%
#   filter(n >= min_occ) %>%
#   pull(species)
# 
# # Extinction matrix:
# #   rows = plots, cols = species historically present in >= min_occ plots
# #   1 = present historically & absent in resurvey (extinct)
# #   0 = all other cases
# extinction_matrix <- full_pa %>%
#   filter(species %in% ext_species) %>%
#   mutate(extinct = as.integer(hist == 1L & rsrv == 0L)) %>%
#   select(plot_id, species, extinct) %>%
#   pivot_wider(names_from = species, values_from = extinct, values_fill = 0L) %>%
#   arrange(match(plot_id, plot_ids_final)) %>%
#   column_to_rownames('plot_id')
# 
# # Colonization matrix:
# #   rows = plots, cols = species that colonized >= min_occ plots
# #   1 = absent historically & present in resurvey (colonized)
# #   0 = all other cases
# colonization_matrix <- full_pa %>%
#   filter(species %in% col_species) %>%
#   mutate(colonized = as.integer(hist == 0L & rsrv == 1L)) %>%
#   select(plot_id, species, colonized) %>%
#   pivot_wider(names_from = species, values_from = colonized, values_fill = 0L) %>%
#   arrange(match(plot_id, plot_ids_final)) %>%
#   column_to_rownames('plot_id')
# 
# cat('Extinction matrix:', nrow(extinction_matrix), 'plots x', ncol(extinction_matrix), 'species\n')
# cat('Colonization matrix:', nrow(colonization_matrix), 'plots x', ncol(colonization_matrix), 'species\n')
# 
# # Save
# saveRDS(extinction_matrix, 'data/model_data/extinction_matrix.rds')
# saveRDS(colonization_matrix, 'data/model_data/colonization_matrix.rds')
# 
# # #OLD! Colonization and extinctions####
# # 
# # colext <- spd %>% 
# #   select(- taxongroup, -layer, -cover) %>%
# #   unique() %>%
# #   mutate(pres=1) %>% spread(resurvey, pres, fill = 0) %>%
# #   mutate(
# #     extinction   = ifelse(`0` == 1 & `1` == 0, 1, 0),
# #     colonization = ifelse(`0` == 0 & `1` == 1, 1, 0),
# #     stable       = ifelse(`0` == `1`, 1, 0)
# #   )
# # 
# # colext$dynamic <- ifelse(colext$extinction == 1, 'extinction', 'other')
# # colext$dynamic <- ifelse(colext$colonization == 1, 'colonization', colext$dynamic)
# # colext$dynamic <- ifelse(colext$stable == 1, 'stable', colext$dynamic)
# # 
# # colext %>% group_by(species, dynamic) %>% summarise(n=n()) %>% ungroup()
# # colext %>% group_by(species, dynamic) %>% summarise(n=n()) %>% ungroup() %>% group_by(species) %>% filter(n==max(n)) %>% filter(n()==1) 
# # 
# # # colext <- colext %>% select(plot_id, species, dynamic) 
# # 
# # colext <- colext %>% select(plot_id, species, extinction, colonization, stable, dynamic) %>%
# #   left_join(
# #     final_merge %>% select(-logRR_S, -jaccard, -bray.curtis)
# #   )
# # 
# # # export
# # colext %>%
# #   write_csv('data/model_data/ext.col.sta_data.csv')
