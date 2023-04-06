# SJV_PWP_AsSe
[![DOI](https://zenodo.org/badge/557502551.svg)](https://zenodo.org/badge/latestdoi/557502551)

Within this repository are all the files needed to replicate the analysis. The  folder structure is:

/SJV_PWP_AsSe <br />
├── geochemical_data <br />
│ ├── CalGEM_WSD.Rda <br />
│ ├── PHREEQC <br />
│ │ ├── As_PHREEQCinput.csv <br />
│ │ ├── As_SatInds.pqi <br />
│ │ ├── As_SatInds.pqo <br />
│ │ ├── As_SatInds.csv <br />
│ │ ├── phreeqc.log <br />
│ │ ├── Se_PHREEQCinput.csv <br />
│ │ ├── Se_SatInds.pqi <br />
│ │ ├── Se_SatInds.pqo <br />
│ │ └── Se_SatInds.csv <br />
│ ├── SJV_ponds_geotracker_excel.xlsx <br />
│ ├── SJV_ponds_geotracker_data.csv <br />
│ ├── SJVPondsDat_raw.Rda <br />
│ └── USGS_NGDB_rock_SJV.Rda <br />
├── geospatial_data <br />
│ ├── AllGWlevels.Rda <br />
│ ├── AllGWStas.Rda <br />
│ ├── maet <br />
│ │ └── CVHM2_et <br />
│ │   └── WY_1989_2019 <br />
│ │     └── ETYYYYMM.txt <br />
│ ├── PRISM_norms <br />
│ │ ├── PRISM_ppt_30yr_norm_800m.tfw <br />
│ │ ├── PRISM_ppt_30yr_norm_800m.tif <br />
│ │ ├── PRISM_ppt_30yr_norm_800m.tif.aux.xml <br />
│ │ └── PRISM_ppt_30yr_norm_800m.tif.ovr <br />
│ ├── shapefiles <br />
│ │ ├── Alluvial_Bnd.dbf <br />
│ │ ├── Alluvial_Bnd.prj <br />
│ │ ├── Alluvial_Bnd.shp <br />
│ │ ├── Alluvial_Bnd.shp.xml <br />
│ │ ├── Alluvial_Bnd.shx <br />
│ │ ├── CA_GWbasins_Albers.CPG <br />
│ │ ├── CA_GWbasins_Albers.dbf <br />
│ │ ├── CA_GWbasins_Albers.prj <br />
│ │ ├── CA_GWbasins_Albers.sbn <br />
│ │ ├── CA_GWbasins_Albers.sbx <br />
│ │ ├── CA_GWbasins_Albers.shp <br />
│ │ ├── CA_GWbasins_Albers.shp.xml <br />
│ │ ├── CA_GWbasins_Albers.shx <br />
│ │ ├── CA_STATSGO2.CPG <br />
│ │ ├── CA_STATSGO2.dbf <br />
│ │ ├── CA_STATSGO2.prj <br />
│ │ ├── CA_STATSGO2.sbn <br />
│ │ ├── CA_STATSGO2.sbx <br />
│ │ ├── CA_STATSGO2.shp <br />
│ │ ├── CA_STATSGO2.shp.xml <br />
│ │ ├── CA_STATSGO2.shx <br />
│ │ ├── CVHM2_modelgrid.CPG <br />
│ │ ├── CVHM2_modelgrid.dbf <br />
│ │ ├── CVHM2_modelgrid.prj <br />
│ │ ├── CVHM2_modelgrid.sbn <br />
│ │ ├── CVHM2_modelgrid.sbx <br />
│ │ ├── CVHM2_modelgrid.shp <br />
│ │ ├── CVHM2_modelgrid.shp.xml <br />
│ │ ├── CVHM2_modelgrid.shx <br />
│ │ ├── SJV_CompSamps.dbf <br />
│ │ ├── SJV_CompSamps.prj <br />
│ │ ├── SJV_CompSamps.shp <br />
│ │ ├── SJV_CompSamps.shx <br />
│ │ ├── SJV_gwWells.dbf <br />
│ │ ├── SJV_gwWells.prj <br />
│ │ ├── SJV_gwWells.shp <br />
│ │ ├── SJV_gwWells.shx.xml <br />
│ │ └── SJV_gwWells.shx <br />
│ ├── SJVspatialData.Rda <br />
│ ├── SJVwaterTable.Rda <br />
│ ├── SoilGrids <br />
│ │ └── PROP_TOP_BOT_mean.tif <br />
│ └── USGS_SoilChem <br />
│   ├── C_As.tif <br />
│   ├── C_Fe.tif <br />
│   ├── C_Mn.tif <br />
│   └── C_se.tif <br />
├── GeospatialAnalysis.R <br />
├── raw_figures <br />
├── Rcustomfunctions.R <br />
├── SanJoaquinValleyPonds.Rmd <br />
├── SJV_GWLevels.R <br />
└── SJV_PWP_AsSe.Rproj <br />
└── tables <br />

Due to file sizes, the "maet", "shapefiles", "SoilGrids", and "USGS_SoilChem" folders are 7zip (7z) files. To use the contents of these folders, unzip them to the appropriate folder as shown in the folder organization above.

To replicate the entire analysis begin by running the "SJV_GWLevels.R" script. Next, run the "GeospatialAnalysis.R" script. After these scripts have been run the rest of the analysis can be replicated by running each chunk within the "SanJoaquinValleyPonds.Rmd" R markdown document. 

If you only wish to replicate the random forest regression model construction, simply run each chunk within the "SanJoaquinValleyPonds.Rmd" R markdown document. The necessary spatial variables output by the "SJV_GWLevels.R" and "GeospatialAnalysis.R" scripts will be loaded by loading the"SJVspatialData.Rda" data file.

To replicate saturation indices calculations run the ".pqi" files in PHREEQC.
