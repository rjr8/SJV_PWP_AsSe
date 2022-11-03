library(tidyverse)
library(lubridate)
library(fs)
library(sf)
library(sp)
library(gstat)
library(raster)


# Establish directory
here::i_am('SJV_Composition.Rproj')

# Load script of custom functions
source('Rcustomfunctions.R')


# Load all GAMA GW stations
load(here('geospatial_data','AllGWStas.Rda')) 

# Load all GAMA GW level measurements
load(here('geospatial_data','AllGWlevels.Rda')) 

# Find data observation periods for each station 
ObsPeriods= 
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
  select(SITE_CODE) %>% pull()

uniqueWells = 
  AvgGWLevels %>%
  filter(!(SITE_CODE %in% dupWells))

SJVwells_sf =
  uniqueWells %>%
  st_as_sf(., coords = c('LONGITUDE', 'LATITUDE'), crs = 4326) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)

st_write(SJVwells_sf, 
         here('geospatial_data/shapefiles/','SJVgwWells.shp'))


## Kriging

# convert sf to sp
SJVwells_sp = as_Spatial(SJVwells_sf)

# create a grid of sampling points 
grd = 
  expand.grid(lon = seq(-120.75, -118.75, by = 0.05), 
              lat = seq(35.0, 37.0, by = 0.05))

grd_sf =
  st_as_sf(grd, coords = c("lon", "lat"), 
           crs = 4326, agr = "constant") %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)

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

# Map it
Krige_sf =
  SJV.WL.krige %>% 
  as.data.frame %>%
  st_as_sf(., coords = c('coords.x1', 'coords.x2'), crs = 3310)
  # project to NAD_1983_California_Teale_Albers (meters)
  
# write as shapefile

st_write(Krige_sf, here('geospatial_data/shapefiles/','WaterTable.shp'))