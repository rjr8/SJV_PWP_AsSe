library(sf)
library(tidyverse)
library(raster)
library(terra)
library(here)

# Establish directory
here::i_am('SJV_PWP_AsSe.Rproj')

# Load script of custom functions
source('Rcustomfunctions.R')


# Clean/Standardize Sample Location Data ----------------------------------

### Calculate Sampling Point Centroids on per Global_ID Basis

# This code should not be rerun, these transformations have already been
# applied to the SJVPondsDat_raw.Rda file. They are present here for 
# illustrative purposes only. If replicating code begin @ line 56

load(here::here('geochemical_data','SJVPondsDat_raw.Rda'))

# Check for missing coordinates
CoordCheck =
  SJVPondsDat.raw %>%
  filter(is.na(Latitude) | is.na(Longitude)) %>%
  dplyr::select(Global_ID, Oil_Field, Facility_Name) %>%
  distinct()

# Calculate centroids
FacCentroids_Tulare = 
  SJVPondsDat.raw %>%
  group_by(Global_ID) %>%
  dplyr::summarise(across(Latitude:Longitude, ~mean(.x, na.rm = T))) %>%
  dplyr::rename(Lat_cent = Latitude,
                Lon_cent = Longitude)

# Replace with centroids
SJVPondsDat.raw = 
  SJVPondsDat.raw %>%
  dplyr::select(-(Latitude:Longitude)) %>%
  left_join(FacCentroids_Tulare, by = c('Global_ID' = 'Global_ID')) %>%
  filter(!(is.na(Lat_cent) | is.na(Lon_cent)))

# Create sf object for mapping
SJVPondsDat_sf =
  SJVPondsDat.raw %>%
  st_as_sf(., coords = c('Lon_cent', 'Lat_cent'), crs = 4326) %>%
  dplyr::select(Int_ID, geometry) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)

# Write as a shapefile
st_write(SJVPondsDat_sf, 
         here::here('geospatial_data/shapefiles/', 'SJVCompSamps.shp'))
# this shapefile was edited with the valley axis side descriptor via
# visual inspection in ArcGIS 


# Load Sample Location Data -----------------------------------------------

# BEGIN RUNNING CODE HERE
SJVPondsDat_sf = 
  read_sf(here::here('geospatial_data/shapefiles/', 'SJVCompSamps.shp')) %>%
  rowid_to_column(., 'Rast_ID') # Add a rownumber identifier to facilitate 
                                # proper joining with extracted raster values

# Process Geospatial Data -------------------------------------------------

#### SoilGrids Soil Properties

# Batch read data, calc thickness weighted average of soil properties
Soil.Properties.data =
  fs::dir_ls(here::here('geospatial_data/SoilGrids/'), 
             regexp = '*.tif', recurse = FALSE) %>%
  map_dfc(~rast(.) %>%
            # Reproject to NAD_1983_California_Teale_Albers (meters)
            terra::project(., 'EPSG:3310') %>% 
            # Extract raster values to points
            terra::extract(., SJVPondsDat_sf)) %>%
  .[!duplicated(as.list(.))] %>% # Remove duplicated ID col
  dplyr::rename(Rast_ID = ID...1) %>%
  pivot_longer(cols = 2:49,
               names_to = c('variable', 'd_start', 'd_end', 'descr'),
               names_sep = '_',
               values_to = 'rast_val') %>%
  dplyr::mutate(across(variable, ~str_replace(.x, 'PH', 'pH_soil')),
                across(starts_with('d_'), as.numeric),
                thickness = d_end - d_start) %>% 
  dplyr::select(-c(descr, d_start, d_end))  %>%
  pivot_wider(names_from = 'variable',
              values_from = 'rast_val') %>%
  # create a thickness weighted average for each
  group_by(Rast_ID) %>%
  dplyr::summarise(across(BD:SOC, ~sum(.x*thickness)/200)) %>%
  # untransform data
  dplyr::mutate(across(BD, ~.x/100),
                across(CEC:SOC, ~.x/10))

##### STATSGO2 Soil Data
STATSGO2.data = 
  st_join(SJVPondsDat_sf,
          read_sf(here::here('geospatial_data/shapefiles/', 
                       'CA_STATSGO2.shp'))) %>%
  st_drop_geometry(.) %>%
  arrange(Int_ID) %>%
  dplyr::select(any_of(c('Rast_ID', 'Int_ID', 'slopegradw', 
                         'aws0150wta', 'hydgrpdcd'))) %>%
  dplyr::rename(slope_wa = slopegradw, # give intuitive names
                aws_0.150 = aws0150wta, 
                hydrgrp = hydgrpdcd)

##### USGS Soil Chem Data

USGS.soilChem.data.raw = 
  fs::dir_ls(here::here('geospatial_data/USGS_SoilChem/'), 
             regexp = '*.tif', recurse = FALSE) %>%
  map_dfc(~rast(.) %>%
            # Reproject to NAD83 / Conus Albers (meters)
            terra::project(., 'EPSG:3310') %>% 
            # Extract raster values to points
            terra::extract(., SJVPondsDat_sf)) %>%
  .[!duplicated(as.list(.))] %>% # Remove duplicated ID col
  dplyr::rename(Rast_ID = ID...1) %>%
  st_drop_geometry(.) %>%
  arrange(Rast_ID) %>%
  dplyr::rename_with(., ~gsub('1', 'R', .x)) %>%
  dplyr::rename_with(., ~gsub('2', 'G', .x)) %>%
  dplyr::rename_with(., ~gsub('3', 'B', .x)) %>%
  unite('As_RGB', starts_with('C_As'), remove = TRUE) %>%
  unite('Fe_RGB', starts_with('C_Fe'), remove = TRUE) %>%
  unite('Mn_RGB', starts_with('C_Mn'), remove = TRUE) %>%
  unite('Se_RGB', starts_with('C_Se'), remove = TRUE) %>%
  dplyr::mutate(across(As_RGB, ~paste0('As_',.x)),
                across(Fe_RGB, ~paste0('Fe_',.x)),
                across(Mn_RGB, ~paste0('Mn_',.x)),
                across(Se_RGB, ~paste0('Se_',.x)))

# Assign concentrations based on RGB value  
vals_key = 
  lapply(USGS.soilChem.data.raw %>% 
           dplyr::select(-Rast_ID), unique) %>%
  lapply(., "length<-", max(lengths(.))) %>%
  as.data.frame() %>%
  bind_cols(
    As_conc = c(6.2, 4.35, 7.25, 5.25, 3.4, 8.5, NA, NA, NA),
    Fe_conc = c(2.495, 2.195, 1.155, 4.08, 3.315, 2.85, 1.575, 1.9, 0.455),
    Mn_conc = c(355.5, 433.5, 198.5, 530, 278.5, 119.5, 657, NA, NA),
    Se_conc = c(0.25, 0.5, 0.1, 0.35, 4.05, NA, NA, NA, NA)) %>%
  pivot_longer(everything(),
               names_to = c("set", ".value"),
               names_pattern = "(.+)_(.+)") %>%
  dplyr::select(-set)

# Clean USGS soil chem data
USGS.soilChem.data =
  USGS.soilChem.data.raw %>% 
  pivot_longer(cols = contains('RGB'),
               names_to = 'set',
               values_to = 'RGB') %>%
  left_join(vals_key, by = c('RGB' = 'RGB')) %>%
  dplyr::select(-RGB) %>%
  pivot_wider(id_cols = Rast_ID,
              names_prefix = 'C_',
              names_from = 'set', 
              values_from = 'conc') %>%
  dplyr::rename_with(., ~str_remove_all(.x, '_RGB'))


##### Climatic Data

## Pull reference ET from CVHM2 

# Read CVHM2 grid cells
SJV.CVHM2.cells = 
  st_join(SJVPondsDat_sf,
          st_transform(
            read_sf(here::here('geospatial_data/shapefiles/', 
                               'CVHM2_modelgrid.shp')), crs = 3310)) %>%
  st_drop_geometry(.) %>%
  arrange(Int_ID) %>%
  dplyr::mutate(cell = paste(ROW, COLUMN_, sep='_'))

# Read raw ET data for water years 1989-2019
raw_et = 
  fs::dir_ls(here::here('geospatial_data/maet/CVHM2_et/WY_1989_2019'),
             regexp = '*.txt', recurse = FALSE) %>%
  map_dfr(~read_csv(.,col_names = FALSE ) %>% 
            dplyr::mutate(row = row_number(), .before = everything()) %>%
            pivot_longer(cols = starts_with('X'), 
                         names_to = 'col', 
                         values_to = 'et'),
          .id = 'source')

# Clean raw ET data
clean_et =
  raw_et %>% 
  dplyr::mutate(across(source, ~str_extract_all(.x, '(?<=T)(\\d)+')),
                across(col, ~as.numeric(str_remove_all(.x, 'X'))),
                cell = paste(row, col, sep='_'),
                year = as.numeric(str_sub(source, 1,4)),
                mo = as.numeric(str_sub(source,-2,-1)),
                water_yr = case_when(mo > 9 ~ year+1,
                                     mo <= 9 ~ year)) %>%
  filter(cell %in% unique(SJV.CVHM2.cells$cell)) %>%
  group_by(water_yr, cell) %>%
  dplyr::summarise(sum_et = sum(et)) %>%
  group_by(cell) %>%
  dplyr::summarise(MAEt = mean(sum_et, na.rm = TRUE))

# Join to samples
SJV.refET = 
  left_join(SJV.CVHM2.cells, clean_et, by = c('cell' = 'cell')) %>%
  dplyr::select(c('Int_ID', 'MAEt'))
  

## Read in precipitation data
Precip.data = 
  rast(here::here('geospatial_data/PRISM_norms/', 
                  'PRISM_ppt_30yr_norm_800m.tif')) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  terra::project(., 'EPSG:3310') %>% 
  # Extract raster values to points
  terra::extract(., SJVPondsDat_sf) %>%
  setNames(., c('Rast_ID', 'MAP'))


##### GAMA groundwater Data
load(here::here('geospatial_data','SJVwaterTable.Rda')) 


#### Combine Spatial Variables for Joining

SJV.spatial.vars =
  SJVPondsDat_sf %>% 
  st_drop_geometry(.) %>%
  left_join(Soil.Properties.data, by = c('Rast_ID' = 'Rast_ID')) %>%
  left_join(SJV.refET, by = c('Int_ID' = 'Int_ID')) %>%
  left_join(Precip.data, by = c('Rast_ID' = 'Rast_ID')) %>%
  left_join(USGS.soilChem.data, by = c('Rast_ID' = 'Rast_ID')) %>%
  left_join(STATSGO2.data %>% dplyr::select(-Rast_ID), 
            by = c('Int_ID' = 'Int_ID')) %>%
  left_join(SJV.WaterTable, by = c('Rast_ID' = 'Rast_ID')) %>%
  dplyr::select(-Rast_ID) %>%
  dplyr::mutate(across(!any_of(c('Int_ID','AxisSide','hydrgrp', 'hydric')), 
                ~signif(.x, digits = 3)),
         across(where(is.character), as.factor))


# Write Joined Geospatial Data --------------------------------------------

# Write data
save(SJV.spatial.vars, 
     file = here::here('geospatial_data','SJVspatialData.Rda'))