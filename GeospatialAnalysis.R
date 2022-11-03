library(sf)
library(tidyverse)

# Establish directory
here::i_am('SJV_Composition.Rproj')

# Load script of custom functions
source('Rcustomfunctions.R')


### Calculate Sampling Point Centroids on per Global_ID Basis

# This code should not be rerun, these transformations have already been
# applied to the SJVPondsDat_raw.Rda file. They are present here for 
# illustrative purposes only

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
  left_join(FacCentroids_Tulare, by = c('Global_ID' = 'Global_ID'))

# Create sf object for mapping
SJVPondsDat_sf =
  SJVPondsDat.raw %>%
  st_as_sf(., coords = c('Lon_cent', 'Lat_cent'), crs = 4326) %>%
  dplyr::select(Int_ID, geometry) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)

# Write as a shapefile
st_write(SJVPondsDat_sf, 
         here('geospatial_data/shapefiles/', 'SJVCompSamps.shp'))

### Process Geospatial Data

#### SoilGrids Soil Properties

# Calc thickness weighted averaged of soil properties
Soil.Properties.data = 
  st_read(here('geospatial_data/shapefiles/',
               'SJVSoilGridsProps.shp'), 
          fid_column_name = 'FID') %>% #include FID as check
  st_drop_geometry(.) %>%
  arrange(Int_ID) %>%
  # add variable identifier, and interval thickness
  bind_cols(tibble(variable = rep(c('BD', 'CEC', 'cf', 'clay',
                                    'pH_soil', 'sand', 'silt', 
                                    'SOC'), 
                                  each = 6, length.out = 81216),
                   thickness = rep(c(5, 100, 15, 30, 10, 40), 
                                   each = 1, length.out = 81216))) %>%
  dplyr::select(-FID) %>%
  pivot_wider(names_from = 'variable',
              values_from = 'RASTERVALU') %>%
  # create a thickness weighted average for each
  group_by(Int_ID) %>%
  dplyr::summarise(across(BD:SOC, ~sum(.x*thickness)/200)) %>%
  # untransform data
  dplyr::mutate(across(BD, ~.x/100),
         across(CEC:SOC, ~.x/10))

##### STATSGO2 Soil Data
STATSGO2.data = 
  st_read(here('geospatial_data/shapefiles/', 'SJVSTATSGO2.shp')) %>% 
  st_drop_geometry(.) %>%
  arrange(Int_ID) %>%
  dplyr::select(any_of(c('Int_ID', 'muaggatt_s', # drop extraneous cols
                         'muaggatt_a', 'muaggatt_h'))) %>%
  dplyr::rename(slope_wa = muaggatt_s, # give intuitive names
         aws_0.150 = muaggatt_a, 
         hydrgrp = muaggatt_h)

##### USGS Soil Chem Data
USGS.soilChem.data.raw = 
  st_read(here('geospatial_data/shapefiles/','SJV_USGS_soil.shp')) %>%
  st_drop_geometry(.) %>%
  arrange(Int_ID) %>%
  dplyr::rename_with(., ~gsub('b1', 'R', .x)) %>%
  dplyr::rename_with(., ~gsub('b2', 'G', .x)) %>%
  dplyr::rename_with(., ~gsub('b3', 'B', .x))
  

# assign concentrations based on RGB value  
Asvals = 
  USGS.soilChem.data.raw %>%
  group_by(R_C_As_CA, G_C_As_CA, B_C_As_CA) %>%
  dplyr::summarise(n()) %>%
  unite('RGB', ends_with(c('_CA')), remove = FALSE) %>%
  ungroup() %>%
  bind_cols(C_As = c(4.35,3.4,5.25,6.2,8.5,7.25)) %>%
  dplyr::mutate(across(RGB, ~paste0(.x, '_As'))) %>%
  dplyr::select(c(RGB, C_As))

Fevals = 
  USGS.soilChem.data.raw %>%
  group_by(R_C_Fe_CA, G_C_Fe_CA, B_C_Fe_CA) %>%
  dplyr::summarise(n()) %>%
  unite('RGB', ends_with(c('_CA')), remove = FALSE) %>%
  ungroup() %>%
  bind_cols(C_Fe = c(0.455,1.9,1.155,1.575,2.195,2.495,4.08,3.315,2.85)) %>%
  dplyr::mutate(across(RGB, ~paste0(.x, '_Fe'))) %>%
  dplyr::select(c(RGB, C_Fe))

Mnvals = 
  USGS.soilChem.data.raw %>%
  group_by(R_C_Mn_CA, G_C_Mn_CA, B_C_Mn_CA) %>%
  dplyr::summarise(n()) %>%
  unite('RGB', ends_with(c('_CA')), remove = FALSE) %>%
  ungroup() %>%
  bind_cols(C_Mn = c(278.5,119.5,198.5,355.5,433.5,657,530)) %>%
  dplyr::mutate(across(RGB, ~paste0(.x, '_Mn'))) %>%
  dplyr::select(c(RGB, C_Mn))

Sevals = 
  USGS.soilChem.data.raw %>%
  group_by(R_C_Se_CA, G_C_Se_CA, B_C_Se_CA) %>%
  dplyr::summarise(n()) %>%
  unite('RGB', ends_with(c('_CA')), remove = FALSE) %>%
  ungroup() %>%
  bind_cols(C_Se = c(0.1,0.675,0.5,0.35,0.25)) %>%
  dplyr::mutate(across(RGB, ~paste0(.x, '_Se'))) %>%
  dplyr::select(c(RGB, C_Se))

USGS.soilChem.data =
  USGS.soilChem.data.raw %>%
  dplyr::select(!MapunitRas) %>%
  unite('RGB_As', ends_with(c('As_CA'))) %>%
  dplyr::mutate(across(RGB_As, ~paste0(.x, '_As'))) %>%
  unite('RGB_Fe', ends_with(c('Fe_CA'))) %>%
  dplyr::mutate(across(RGB_Fe, ~paste0(.x, '_Fe'))) %>%  
  unite('RGB_Mn', ends_with(c('Mn_CA'))) %>%
  dplyr::mutate(across(RGB_Mn, ~paste0(.x, '_Mn'))) %>% 
  unite('RGB_Se', ends_with(c('Se_CA'))) %>%
  dplyr::mutate(across(RGB_Se, ~paste0(.x, '_Se'))) %>%
  left_join(Asvals, by = c('RGB_As' = 'RGB')) %>%
  left_join(Fevals, by = c('RGB_Fe' = 'RGB')) %>%
  left_join(Mnvals, by = c('RGB_Mn' = 'RGB')) %>%
  left_join(Sevals, by = c('RGB_Se' = 'RGB')) %>%
  dplyr::select(!starts_with('RGB'))
  


##### Climatic Data

## Pull reference ET from CVHM2 

# Read CVHM2 grid cells
SJV.CVHM2.cells = 
  st_read(here('geospatial_data/shapefiles/','SJV_CVHM2_cells.shp')) %>%
  st_drop_geometry(.) %>%
  arrange(Int_ID) %>%
  dplyr::mutate(cell = paste(ROW, COLUMN_, sep='_'))

# Read raw ET data for water years 1989-2019
raw_et = 
  fs::dir_ls(here('geospatial_data/maet/CVHM2_et/WY_1989_2019'),
             regexp = '*.txt', recurse = FALSE) %>%
  map_dfr(~read_csv(.,col_names = FALSE ) %>% 
            dplyr::mutate(row = row_number(), .before = everything()) %>%
            pivot_longer(cols = starts_with('X'), 
                         names_to = 'col', 
                         values_to = 'et'),
          .id = 'source')
  
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

# join to samples
SJV.refET = 
  left_join(SJV.CVHM2.cells, clean_et, by = c('cell' = 'cell')) %>%
  dplyr::select(c('Int_ID', 'MAEt'))
  

## Read in precipitation data

Precip.data = 
  st_read(here('geospatial_data/shapefiles', 'SJVMAP.shp')) %>%
  st_drop_geometry(.) %>%
  dplyr::rename(MAP = RASTERVALU)


##### GAMA groundwater Data
WaterTable.data = 
  st_read(here('geospatial_data/shapefiles','SJVwaterTable.shp')) %>%
  st_drop_geometry(.) %>%
  arrange(Int_ID)


#### Combine Spatial Variables for Joining

SJV.spatial.vars =
  Soil.Properties.data %>%
  left_join(SJV.refET, by = c('Int_ID' = 'Int_ID')) %>%
  left_join(Precip.data, by = c('Int_ID' = 'Int_ID')) %>%
  left_join(USGS.soilChem.data, by = c('Int_ID' = 'Int_ID')) %>%
  left_join(STATSGO2.data, by = c('Int_ID' = 'Int_ID')) %>%
  left_join(WaterTable.data, by = c('Int_ID' = 'Int_ID')) %>%
  #dplyr::mutate(across(hydric, as_factor)) %>%
  dplyr::mutate(across(!any_of(c('Int_ID','AxisSide','hydrgrp', 'hydric')), 
                ~signif(.x, digits = 3)),
         across(where(is.character), as.factor)) %>%
  dplyr::rename(WLBLS = RASTERVALU)

# Write data
save(SJV.spatial.vars, 
     file = here('geospatial_data','SJVspatialData.Rda'))
