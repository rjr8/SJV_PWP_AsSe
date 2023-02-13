# SJV_PWP_AsSe

Within this repository are all the files needed to replicate the analysis. The  folder structure is:

/SJV_PWP_AsSe
├── geochemical_data
│ ├── CalGEM_WSD.Rda
│ ├── PHREEQC
│ │ ├── As_PHREEQCinput.csv
│ │ ├── As_SatInds.pqi
│ │ ├── As_SatInds.pqo
│ │ ├── As_SatInds.csv
│ │ ├── phreeqc.log
│ │ ├── Se_PHREEQCinput.csv
│ │ ├── Se_SatInds.pqi
│ │ ├── Se_SatInds.pqo
│ │ └── Se_SatInds.csv
│ ├── SJV_ponds_geotracker_excel.xlsx
│ ├── SJV_ponds_geotracker_data.csv
│ ├── SJVPondsDat_raw.Rda
│ └── USGS_NGDB_rock_SJV.Rda
├── geospatial_data
│ ├── AllGWlevels.Rda
│ ├── AllGWStas.Rda
│ ├── maet
│ │ └── CVHM2_et
│ │   └── WY_1989_2019
│ │     └── ETYYYYMM.txt
│ ├── PRISM_norms
│ │ ├── PRISM_ppt_30yr_norm_800m.tfw
│ │ ├── PRISM_ppt_30yr_norm_800m.tif
│ │ ├── PRISM_ppt_30yr_norm_800m.tif.aux.xml
│ │ └── PRISM_ppt_30yr_norm_800m.tif.ovr
│ ├── shapefiles
│ │ ├── Alluvial_Bnd.dbf
│ │ ├── Alluvial_Bnd.prj
│ │ ├── Alluvial_Bnd.shp
│ │ ├── Alluvial_Bnd.shp.xml
│ │ ├── Alluvial_Bnd.shx
│ │ ├── CA_GWbasins_Albers.CPG
│ │ ├── CA_GWbasins_Albers.dbf
│ │ ├── CA_GWbasins_Albers.prj
│ │ ├── CA_GWbasins_Albers.sbn
│ │ ├── CA_GWbasins_Albers.sbx
│ │ ├── CA_GWbasins_Albers.shp
│ │ ├── CA_GWbasins_Albers.shp.xml
│ │ ├── CA_GWbasins_Albers.shx
│ │ ├── CA_STATSGO2.CPG
│ │ ├── CA_STATSGO2.dbf
│ │ ├── CA_STATSGO2.prj
│ │ ├── CA_STATSGO2.sbn
│ │ ├── CA_STATSGO2.sbx
│ │ ├── CA_STATSGO2.shp
│ │ ├── CA_STATSGO2.shp.xml
│ │ ├── CA_STATSGO2.shx
│ │ ├── CVHM2_modelgrid.CPG
│ │ ├── CVHM2_modelgrid.dbf
│ │ ├── CVHM2_modelgrid.prj
│ │ ├── CVHM2_modelgrid.sbn
│ │ ├── CVHM2_modelgrid.sbx
│ │ ├── CVHM2_modelgrid.shp
│ │ ├── CVHM2_modelgrid.shp.xml
│ │ ├── CVHM2_modelgrid.shx
│ │ ├── SJV_CompSamps.dbf
│ │ ├── SJV_CompSamps.prj
│ │ ├── SJV_CompSamps.shp
│ │ ├── SJV_CompSamps.shx
│ │ ├── SJV_gwWells.dbf
│ │ ├── SJV_gwWells.prj
│ │ ├── SJV_gwWells.shp
│ │ ├── SJV_gwWells.shx.xml
│ │ └── SJV_gwWells.shx
│ ├── SJVspatialData.Rda
│ ├── SJVwaterTable.Rda
│ ├── SoilGrids
│ │ └── PROP_TOP_BOT_mean.tif
│ └── USGS_SoilChem
│   ├── C_As.tif
│   ├── C_Fe.tif
│   ├── C_Mn.tif
│   └── C_se.tif
│ └── ProximityAnalysisFunctions.R
├── GeospatialAnalysis.R
├── raw_figures
├── Rcustomfunctions.R
├── SanJoaquinValleyPonds.Rmd
├── SJV_GWLevels.R
└── SJV_PWP_AsSe.Rproj
└── tables

Due to file sizes, the "maet", "shapefiles", and "SoilGrids" folders are 7zip (7z) files. To use the contents of these folders, unzip them to the appropriate folder as shown in the folder organization above.

To replicate the entire analysis begin by running the "SJV_GWLevels.R" script. Next, run the "GeospatialAnalysis.R" script. After these scripts have been run the rest of the analysis can be replicated by running each chunk within the "SanJoaquinValleyPonds.Rmd" R markdown document. 

If you only wish to replicate the random forest regression model construction, simply run each chunk within the "SanJoaquinValleyPonds.Rmd" R markdown document. The necessary spatial variables output by the "SJV_GWLevels.R" and "GeospatialAnalysis.R" scripts will be loaded by loading the"SJVspatialData.Rda" data file.
