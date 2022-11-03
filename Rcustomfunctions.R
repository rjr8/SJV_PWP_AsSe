## Setup

# Check for installed packages
setupPackages = function(req.packages) {
  
  # This function searches through your list of installed 
  # packages to see if the necessary packages are installed.  
  # If they are not, they are installed. The only thing that
  # needs defined is the list of packages you wish to use. 
  #
  # req.packages: the list of packages used. This can be user
  # 	            defined in the function call (e.g., 
  #               list.of.packages = c("x", "y")), or as a 
  #               variable prior to calling the function
  
  new.packages = req.packages[!(req.packages 
                                %in% installed.packages()
                                [, "Package"])]
  
  if(length(new.packages)) install.packages(new.packages)
  
  #load packages
  lapply(req.packages, library, character.only = TRUE)
  
}

batchReadDfs = function(path,pattern){
  
  fs::dir_ls(path, regexp = pattern, recurse = TRUE) %>%
    map_dfr(read_csv)
  
}


# Batch load dfs into one dataframe
BatchReadBindDfs = function(path, pattern, coltype) {
  
  # This function performs a batch import of files within a
  # folder, and appends them by row in a data frame, 
  # resulting in a single data frame.  
  #
  # path:     The full path from your home directory to 
  #           your results folder (e.g., "~/gis/MyFiles").
  # pattern:  The shared naming convention and extension of 
  #           the result files (e.g., [whatever].csv).
  # coltype:  The format of the data in each column 
  #           (e.g., factor, double, character).
  
  files = list.files(path, pattern, full.names = TRUE)
  lapply(files, read_csv, col_types = coltype) %>%
    bind_rows()
}


batchResultsImport = function(namingPattern, resultsFolder,
                              fullResultsPath, fileExtension, dfs){
  
  # This function performs a batch import of files in a 
  # folder, and saves them as separate data frames. The 
  # user supplied arguments define the result naming  
  # conventions, where the files live (both the folder 
  # and the full path from your home directory), and 
  # the extension of the results tables.  
  #
  # namingPattern:    The shared naming convention and 
  #                   extension of the result files 
  #                   (e.g., [whatever]ftWells.txt).
  # resultsFolder:    The folder name where the results 
  #                   live.
  # fullResultsPath:  The full path from your home 
  #                   directory to your results folder 
  #                   (e.g., "~/gis/MyResults").
  # fileExtension:    The extension of the result files
  #                   (e.g., ".txt").
  # dfs:              The name of the vector that will 
  #                   hold the names of each result 
  #                   data frame. 
  
  # Retrieve list of results file names
  file_names = 
    list.files(path = resultsFolder, 
               pattern = namingPattern,
               full.names = FALSE)
  
  # Create a vector of data frame names from file names
  dfnames = tools::file_path_sans_ext(file_names)
  
  
  for (i in dfnames){
    
    path = file.path(fullResultsPath, 
                     paste0(i, fileExtension))
    
    assign(paste0("df",i), 
           read_csv(path), 
           envir = .GlobalEnv)
  }
  
  # make this a user defined name so you can run 
  # multiple results tables at once
  assign(dfs, dfnames, envir = .GlobalEnv)
  
}


## Stats helpers 

# Calculate standard error of the mean
SEM = function(x) {
  
  # This function calculates the standard error of the 
  # mean (stdev/(# of samples)^-.5). This calculation 
  # excludes NA values.
  #
  # x: The vector you'd like to calculate SEM for.
  
  sd(x, na.rm= TRUE)/sqrt(sum(!is.na(x)))
  
}

## Geochem analysis

# Create a data frame of all non detect values
createCleanandNDdfs = function(df.raw, name, factorCols, 
                      colsToDrop, firstChemDatColname){
  
  # This function creates two data frames from one input
  # data frame. The first data frame is a clean version 
  # of the data which does not contain any non detect 
  # values. The second data frame contains only the non 
  # detect values. Prior to creating these data frames
  # user defined columns are converted to factors,  
  # columns containing extraneous information are 
  # dropped, and columns with chemical data are  
  # converted to numerics (if necessary). Output data  
  # frames contain the user defined root and ".Clean"  
  # and ".ND" suffixes.
  #
  # df.raw:       The data frame with raw chem data  
  # name:         The root name of data frames
  # factorCols:   A vector of the names of columns that 
  #               are categorical variables. If none, 
  #               use an empty vector
  # colsToDrop:   A vector of the names of columns to 
  #               drop. If none, use an empty vector
  # firstChemDatColname: The column indcies which contain chem
  #               data (e.g., 10:35)
  
  
  # create vector of all radionuclide cols
  radios = c("Gross_alpha", "Gross_beta", 
             "Ra_226", "Ra_228", "U")
  
  # Create a vector of carbonate species
  carbs = c(HCO3 = NA_real_, CO3 = NA_real_, 
            OH = NA_real_, Hardness_CaCO3 = NA_real_)
  
  # Create basic df that has all shared operations
  df.temp = 
    df.raw %>%
    mutate(across(any_of(factorCols), as.factor)) %>%
    # Drop columns that are not needed for analysis
    dplyr::select(-any_of(colsToDrop)) %>%
    # Check if HCO3, CO3, and/or OH cols are present,
    # add to data frame if not.
    add_column(., !!!carbs[!names(carbs) %in% names(.)])
  
  chemDatStart = which(colnames(df.temp) == firstChemDatColname)
  chemDataCols = c(chemDatStart:length(df.temp))
  
  df.Clean = 
    df.temp %>%
    # Remove non detections, anything < will become an NA
    mutate(across(chemDataCols, as.numeric)) %>%
    # change 0s and negatives to NAs
    mutate(across(chemDataCols, ~na_if(.x, 0)),
           across(any_of(radios), 
                  ~if_else(.x<=0, NA_real_, as.numeric(.x)))) %>%
    # Calculate missing HCO3, CO3, & OH values
    mutate(across(HCO3, 
                  ~if_else(is.na(.x), 1.22*HCO3Alk_CaCO3, 
                           as.numeric(.x))),
           across(CO3, 
                  ~if_else(is.na(.x), 0.6*CO3Alk_CaCO3,
                           as.numeric(.x))),
           across(OH,
                  ~if_else(is.na(.x), 0.34*OHAlk_CaCO3,
                           as.numeric(.x)))) %>%
#    mutate(Ca.Sr = (Ca/40.08)/(Sr/87.62),
#           Ca.Mg = (Ca/40.08)/(Mg/24.31),
#           Ba.Sr = (Ba/137.33)/(Sr/87.62),
#           Na.Cl = (Na/22.99)/(Cl/35.45)) %>%
    # calculate missing hardness values
    mutate(across(Hardness_CaCO3,
                  ~if_else(is.na(.x), (2.5*Ca)+(4.1*Mg), as.numeric(.x))))
  
  
  df.ND =
    df.temp %>%
    mutate(across(chemDataCols, 
                  ~if_else(str_detect(.x, "^<"), 
                           as.character(.x), NA_character_)),
           across(any_of(radios), 
                  ~if_else(.x<=0, "<bkg", NA_character_)))
  
  assign(paste0(name, ".Clean"), df.Clean, envir = .GlobalEnv)
  assign(paste0(name, ".ND"), df.ND, envir = .GlobalEnv)

}

# Count the number of detections for any analyte
countDetects = function(x) {
  
  # This function counts the number of detections for a
  # specified analyte. 
  #
  # x: the column/analyte you want to count the number  
  # 	 of detections of.
  
  sum(!is.na(x))
  
  }

# Calculate Charge Balance Error
calcChargeBalance = function(df, CBerrorFlag, 
                             uniqueIdentifier){
  
  # This function calculates the charge balance error 
  # for each sample. If the error exceeds a user   
  # defined threshold (e.g., 5%), the sample is 
  # flagged. 
  #
  # df:               The data frame with chemical data 
  # CBerrorFlag:      The acceptable charge balance 
  #                   error threshold
  # uniqueIdentifier: The column that is the unique 
  #                   identifier for each observation   

  
  x = c("PeriodicTable")
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
  
  # Establish user selected error flag
  flag = CBerrorFlag
  
  # Create vectors of each ion category
  pos.monovalent = c('Na', 'K')
  pos.divalent = c('Ca','Mg','Sr','Ba')
  
  neg.monovalent = c('Cl', 'HCO3', 'Br', 'NO3')
  neg.divalent = c('CO3','SO4')
  
  # Create vector of all ions
  all.ions = 
    list(pos.monovalent, pos.divalent,
         neg.monovalent, neg.divalent) %>%
    unlist(.)
  
  temp2 = 
    df %>%
    # pull unique identifier and ion columns
    dplyr::select(any_of(c(uniqueIdentifier, all.ions))) %>%
    pivot_longer(cols = !(!! sym(uniqueIdentifier)),
                 names_to = "ion",
                 values_to = "conc") %>%
    # Classify ions by type and valence
    mutate(temp.col = case_when(ion %in% pos.monovalent ~ "cat,1",
                                ion %in% neg.monovalent ~ "an,-1",
                                ion %in% pos.divalent ~ "cat,2",
                                ion %in% neg.divalent ~ "an,-2")) %>%
    # Separate classification column
    separate(temp.col, into = c("type", "valence"), 
             sep = ",", remove = TRUE, convert = TRUE) %>%
    # Find molar masses
    mutate(mass = mass(ion)) %>%
    # Cleaning up after PeriodicTable's mess...
    mutate(across(mass, 
                  ~case_when(ion == 'HCO3' ~ 61.01684,
                             ion == 'CO3' ~ 60.0089,
                             ion == 'NO3' ~ 62.0049,
                             ion == 'SO4' ~ 96.0626,
                             TRUE ~ as.numeric(.x)))) %>%
    # Convert from ppm to meq/L
    mutate(meq = (conc*valence)/mass) %>%
    group_by(!! sym(uniqueIdentifier)) %>%
    # Calculate the charge balance error
    summarise(CBE = (sum(meq, na.rm = TRUE)/sum(abs(meq), na.rm = TRUE))*100) %>%
    mutate(CBE_flag = case_when(abs(CBE) <= flag ~ 'N',
                                TRUE ~ 'Y'))
  
  df = 
    df %>% left_join(temp2, by = uniqueIdentifier)

}

# Calculate Molar Elemental Ratios
calcElementalRatios = function(df){
  
  # This function calculates molar ratios for a variety
  # of elements.   
  #
  # df: The dataframe of chemical data. This function 
  #     assumes that column names are the chemical 
  #     formula or symbol (e.g., "SO4" rather than  
  #     "sulfate", "Pb" rather than "lead").
  
  df %>%
    mutate(Ca.Sr = (Ca/40.078)/(Sr/87.620),
           Ca.Mg = (Ca/40.078)/(Mg/24.305),
           Ba.Sr = (Ba/137.327)/(Sr/87.620),
           Na.Cl = (Na/22.98977)/(Cl/35.453),
           Ba.Cl = (Ba/137.327)/(Cl/35.453),
           SO4.Cl = (SO4/96.06)/(Cl/35.453),
           Mg.Na = (Mg/24.305)/(Na/22.98977))
}
