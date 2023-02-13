library(tidyverse)
library(lubridate)
library(fs)
library(sf)
library(sp)
library(here)
library(gstat)
library(raster)
library(terra)



# Setup -------------------------------------------------------------------

# Establish directory
here::i_am('SJV_PWP_AsSe.Rproj')

# Load script of custom functions
source('Rcustomfunctions.R')


# Load and Process GAMA Data ----------------------------------------------


# Load all GAMA GW stations
load(here::here('geospatial_data','AllGWStas.Rda')) 

# Load all GAMA GW level measurements
load(here::here('geospatial_data','AllGWlevels.Rda')) 

# Find data observation periods for each station 
ObsPeriods = 
  AllGWlevels %>%
  group_by(SITE_CODE) %>%
  summarise(startDate = min(MSMT_DATE),
            endDate = max(MSMT_DATE)) %>%
  mutate(PERIOD = (endDate - startDate)/365.25) 

Stas =
  ObsPeriods %>%
  filter(PERIOD >=10) %>%
  left_join(AllGWStas, by = c('SITE_CODE' = 'SITE_CODE')) %>%
  filter(endDate > as.Date('01/01/2000', format = "%m/%d/%Y"))

AvgGWLevels = 
  AllGWlevels %>%
  filter(is.na(WLM_QA_DESC)) %>%
  right_join(Stas, by = c('SITE_CODE' = 'SITE_CODE')) %>%
  filter(RDNG_RP < WELL_DEPTH) %>%
  mutate(WLBLS = case_when(RPE > GSE ~ abs((RPE-RDNG_RP)-(RPE-GSE)),
                           GSE >= RPE ~ abs((RPE-RDNG_RP)+(GSE-RPE)))) %>%
  group_by(SITE_CODE) %>%
  summarise(WLBLS_avg = signif(mean(WLBLS, na.rm = TRUE)/3.28, digits = 3)) %>%
  left_join(Stas, by = c('SITE_CODE' = 'SITE_CODE')) 
  
dupWells = # find duplicated wells, pull site IDs of deeper duplicates
  AvgGWLevels %>%
  group_by(LATITUDE, LONGITUDE) %>%
  filter(n() >1) %>%
  group_by(LATITUDE, LONGITUDE) %>%
  arrange(WELL_DEPTH, .by_group = TRUE) %>%
  slice(2:n()) %>%
  dplyr::select(SITE_CODE) %>% pull()

uniqueWells = 
  AvgGWLevels %>%
  filter(!(SITE_CODE %in% dupWells))

SJVwells_sf =
  uniqueWells %>%
  st_as_sf(., coords = c('LONGITUDE', 'LATITUDE'), crs = 4326) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)


# Write well data to file
st_write(SJVwells_sf, 
         here::here('geospatial_data/shapefiles/','SJVgwWells.shp'))



# Interpolate Water Table Surface via Kriging -----------------------------


# convert sf to sp
SJVwells_sp = as_Spatial(SJVwells_sf)

# create a grid of sampling points 
grd = 
  expand.grid(y = seq(-335250, -112250, by = 500), 
              x = seq(-68500, 114000, by = 500))

grd_sf =
  st_as_sf(grd, coords = c("x", "y"), crs = 3310, agr = "constant") 

grd_sp = as_Spatial(grd_sf)

# create variogram
SJV.vgm = variogram(WLBLS_avg~1, SJVwells_sp)

# fit a variogram model
SJV.fit = 
  fit.variogram(SJV.vgm, model = vgm(c('Exp', 'Gau', 'Bes'))) 
# vgm() for a list of models

plot(SJV.vgm, SJV.fit)

# Krige it
SJV.WL.krige = 
  krige((WLBLS_avg)~1, SJVwells_sp, grd_sp, model = SJV.fit)


# Convert krige results to SpatVector 
Krige_sv = 
  st_as_sf(SJV.WL.krige, coords = c('coords.x1', 'coords.x2'), crs = 3310) %>%
  dplyr::select(-var1.var) %>%
  terra::vect(.)

# Create a 500m resolution raster grid
SJV_WaterTable_cells = terra::rast(Krige_sv, res = 500)

# Fill grid with krige results, extract interpolated water 
# table surface to sample locations
SJV.WaterTable = 
  terra::rasterize(Krige_sv, SJV_WaterTable_cells, 
                   "var1.pred", touches = TRUE) %>%
  terra::extract(., 
                 read_sf(here::here('geospatial_data/shapefiles/', 
                                    'SJVCompSamps.shp'))) %>%
  setNames(., c('Rast_ID', 'WLBLS'))

# write data
save(SJV.WaterTable, 
     file = here::here('geospatial_data','SJVwaterTable.Rda'))


# gw qual -----------------------------------------------------------------

GAMA_WQ = 
  fs::dir_ls(here::here('geochemical_data/GW_data/'),
             regexp = '*.txt', recurse = FALSE) %>%
  map_dfr(~read_tsv(., col_names = TRUE) %>%
            filter(str_detect(GM_CHEMICAL_NAME, 
                              '[Aa]rsenic|[Ll]ithium|[Ss]elenium')) %>%
            filter(GM_RESULT_MODIFIER == '='),
        .id = 'source') %>%
  dplyr::mutate(across(source, 
                ~str_extract_all(.x, 
                                 'fresno|madera|mariposa|merced|kings|kern|tulare')))

GAMA_WQ_summ = 
  GAMA_WQ %>%
  select(c(GM_WELL_ID, GM_SAMP_COLLECTION_DATE, 
           GM_CHEMICAL_NAME, GM_RESULT)) %>%
  pivot_wider(id_cols = c(GM_WELL_ID, GM_SAMP_COLLECTION_DATE),
              names_from = GM_CHEMICAL_NAME,
              values_from = GM_RESULT,
              values_fn = mean) 

save(GAMA_WQ_summ, 
     file = here::here('geochemical_data/','GAMAWQData.Rda'))
