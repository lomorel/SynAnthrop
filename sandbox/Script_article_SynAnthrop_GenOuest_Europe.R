# SynAnthrop package, draft script written by B. Bongibault
# available for minor corrections before adding portions to the main vignettes
# This script is in R format because it is executed on GenOuest's HPC.

# Load packages
library(terra) # spatial analyses
library(sf)                 #for manipulating downloaded maps
library(tidyverse)          #for tidy analysis
library(ggplot2) # draw graphs
library(dplyr) # pipe %>%, do ctrl+maj+M to compute automatically
library(magrittr) # "ceci n'est pas un pipe" double pipe %<>%
library(future.apply) # parallelize computation
library(here) # path to your files starts where the Rproj file is

plan(multicore, workers = 1)


# Load biodiversity data ################################################################

# Select the dataset to use in this instance of Synanthrop 
# datasets available are: Amphibian, Squamata, Mammal, and Aves
# Datasets are not available on GitHub, but independently through large file transfers.
# Add the datasets in the "data" folder of SynAnthrop to run the analysis.

# FOR INTERACTIVE VERSIONS, uncomment this line:
taxa <- readline(prompt = "Enter the taxa that you want, among Amphibian, Squamata, Mammal, and Aves. The string must match exactly these propositions.")
# FOR NON-INTERACTIVE VERSIONS, comment the previous line and choose directly the dataset that you want with the following line:
# taxa = "Squamata" # replace with the dataset name: Amphibian, Squamata, Mammal, or Aves.

# then execute the next line:
if (taxa == "Amphibian"){
  spOcc = vect(file.path(here(), "data", "Amphibia",
                             "Donnees_Amphibia_clean.shp"))
  
  # Tidy the dataset for analysis
  spOcc %<>% as.data.frame(geom = "XY") %>% # add coordinates as columns
    rename(X = x, Y = y) %>% # change the column names to capital letters
    # filter out erroneous data identified elsewhere:
    filter(!datasetKey %in% c("da36b9f6-5acd-4eae-af23-f49be4ed330e", 
                              "f3278405-0943-41a0-9f04-3a7733ec344d",
                              "2dd3b5c1-2a7f-4324-83d3-d05cff85676b",
                              "ba8b3250-6115-4680-8c21-6129c35ba738",
                              "bd067a4c-05b8-4def-925c-60ce83c01891",
                              "cd9e0589-187b-43c0-9692-50d8db7d0467",
                              "da36b9f6-5acd-4eae-af23-f49be4ed330e",
                              "dd915a74-f2d2-4888-a62d-4a5ebeae4e0e")) %>%
    filter(!(datasetKey == "8a863029-f435-446a-821e-275f4f641165" & countryCod == "NL"))
  
  # Sanity check: are there the expected number of observations after filtering?
  dim(spOcc)[1]
  #TODO replace the line above by the line below 
  # if(dim(spOcc)[1] == **replace this with the number of rows**) { print("Amphibian dataset ready.")} else {print("Error during Amphibian dataset setup. Try again")}
  
  
  
} else if (taxa == "Squamata"){
  spOcc = vect(file.path(here(), "data", "squamata", 
                            "Donnees_Squamata_clean.shp"))
  
  # Tidy the dataset for analysis following the same logic as amphibians
  spOcc %<>% as.data.frame(geom = "XY") %>% 
    rename(X = x, Y = y) %>% 
    filter(!datasetKey %in% c("da36b9f6-5acd-4eae-af23-f49be4ed330e",
                              "f3278405-0943-41a0-9f04-3a7733ec344d",
                              "2dd3b5c1-2a7f-4324-83d3-d05cff85676b",
                              "ba8b3250-6115-4680-8c21-6129c35ba738",
                              "bd067a4c-05b8-4def-925c-60ce83c01891",
                              "cd9e0589-187b-43c0-9692-50d8db7d0467",
                              "da36b9f6-5acd-4eae-af23-f49be4ed330e",
                              "dd915a74-f2d2-4888-a62d-4a5ebeae4e0e",
                              "24168b97-a6b0-41cf-9836-29110572efc6",
                              "45cff797-6ed9-4a25-b659-30aa4caa46f2",
                              "72ea9678-7174-44c4-9535-de9cb733a15d",
                              "d3b6cb30-0a64-4f82-91ea-0bb14637ee17",
                              "faf313a1-9ae4-43f4-bfc9-974281feac0e")) %>% 
    filter(!(datasetKey == "8a863029-f435-446a-821e-275f4f641165" & countryCod == "NL"))
    
  # Sanity check: 
  if(dim(spOcc)[1] == 605504) { print("Squamata dataset ready.")} else {print("Error during Squamata dataset setup. Try again")}
    
  
} else if (taxa == "Mammal"){
  spOcc = vect(file.path(here(), "data", "Mammal", 
                            "Donnees_Mammal_clean.shp"))
  
  # Tidy the dataset for analysis following the same logic as amphibians
  spOcc %<>% as.data.frame(geom = "XY") %>% 
    rename(X = x, Y = y) %>% 
    filter(!datasetKey %in% c("2dd3b5c1-2a7f-4324-83d3-d05cff85676b",
                              "37f56892-bd85-479a-803b-7fbfdfd6886c",
                              "584f899d-a3fd-483a-aba7-514f154feade",
                              "5e48f65e-bd59-4a14-a21b-fcd761ac380f",
                              "66e240a7-ae35-4f59-9ad9-c62687015ed8",
                              "6fc1275d-e647-4e32-835a-3cc94ee2be33",
                              "72ea9678-7174-44c4-9535-de9cb733a15d",
                              "7c4c67d2-c554-4f79-b933-89b74873e8b8",
                              "86cb4fac-86fd-415b-a8f4-1349a3a71c09",
                              "a4a74db9-fe84-40f3-8782-e4940416076a",
                              "b4c2c8f8-3526-4873-9867-b28e260a12a0",
                              "ba8b3250-6115-4680-8c21-6129c35ba738",
                              "bd067a4c-05b8-4def-925c-60ce83c01891",
                              "bdb525f6-2087-4f7f-8ab5-b5b052898996",
                              "c3b0e0ff-def0-40dd-b72d-cbe5e79c1213",
                              "c6bbb6ef-ad16-4f3c-99e2-f693760173e0",
                              "cd9e0589-187b-43c0-9692-50d8db7d0467",
                              "dd915a74-f2d2-4888-a62d-4a5ebeae4e0e",
                              "e2544278-4e9d-4c3f-8b44-a05fa07488c6",
                              "e4af8f92-af98-4649-93ef-83dc99169614",
                              "ea0b6be3-72eb-41d9-b8cc-e9ae2f53f93e",
                              "ef2ad061-ab4d-4b5a-b452-04bf4ac078c2",
                              "f3278405-0943-41a0-9f04-3a7733ec344d",
                              "f6637736-f04e-4361-87a0-541e08a8a8d3",
                              "f9b241f8-900c-4077-ab8a-cca18427fc43")) %>% 
    filter(!(datasetKey == "8a863029-f435-446a-821e-275f4f641165" & countryCod == "NL")) %>% 
    # plus, we are not interested in domestic and marine species that are also excluded
    filter(!species %in% c("Halichoerus grypus", "Delphinus delphis", "Orcinus orca",
                           "Phoca vitulina", "Phocoena phocoena", "Stenella coeruleoalba",
                           "Tursiops truncatus", "Homo sapiens", "Felis catus", "Bos taurus"))
  
  # Sanity check: 
  dim(spOcc)[1]
  #TODO replace the line above by the line below 
  # if(dim(spOcc)[1] == **replace this with the number of rows**) { print("Mammal dataset ready.")} else {print("Error during Mammal dataset setup. Try again")}
  
  
  
} else if (taxa == "Aves"){  
  spOcc = vect(file.path(here(), "data", "Aves", 
                          "Donnees_Aves_clean.shp"))
  
  # Tidy the dataset for analysis following the same logic as amphibians
  spOcc %<>% as.data.frame(geom = "XY") %>% 
    rename(X = decimalLongitude, Y = decimalLatitude,
           individual = individualCount,
           countryCod = countryCode,
           coordinate = coordinateUncertaintyInMeters) %>% 
    filter(taxonRank == "SPECIES" & occurrenceStatus == "PRESENT") %>% # keep only presences at the species level
    filter(month != 7) %>%  #TODO: Baptiste, I don't know what this is for 
    filter(!datasetKey %in% c("2dd3b5c1-2a7f-4324-83d3-d05cff85676b",
                         "37f56892-bd85-479a-803b-7fbfdfd6886c",
                         "584f899d-a3fd-483a-aba7-514f154feade",
                         "5e48f65e-bd59-4a14-a21b-fcd761ac380f",
                         "66e240a7-ae35-4f59-9ad9-c62687015ed8",
                         "6fc1275d-e647-4e32-835a-3cc94ee2be33",
                         "72ea9678-7174-44c4-9535-de9cb733a15d",
                         "7c4c67d2-c554-4f79-b933-89b74873e8b8",
                         "86cb4fac-86fd-415b-a8f4-1349a3a71c09",
                         "a4a74db9-fe84-40f3-8782-e4940416076a",
                         "b4c2c8f8-3526-4873-9867-b28e260a12a0",
                         "ba8b3250-6115-4680-8c21-6129c35ba738",
                         "bd067a4c-05b8-4def-925c-60ce83c01891",
                         "bdb525f6-2087-4f7f-8ab5-b5b052898996",
                         "c3b0e0ff-def0-40dd-b72d-cbe5e79c1213",
                         "c6bbb6ef-ad16-4f3c-99e2-f693760173e0",
                         "cd9e0589-187b-43c0-9692-50d8db7d0467",
                         "dd915a74-f2d2-4888-a62d-4a5ebeae4e0e",
                         "e2544278-4e9d-4c3f-8b44-a05fa07488c6",
                         "e4af8f92-af98-4649-93ef-83dc99169614",
                         "ea0b6be3-72eb-41d9-b8cc-e9ae2f53f93e",
                         "ef2ad061-ab4d-4b5a-b452-04bf4ac078c2",
                         "f3278405-0943-41a0-9f04-3a7733ec344d",
                         "f6637736-f04e-4361-87a0-541e08a8a8d3",
                         "f9b241f8-900c-4077-ab8a-cca18427fc43",
                         
                         "da36b9f6-5acd-4eae-af23-f49be4ed330e",
                         "f3278405-0943-41a0-9f04-3a7733ec344d",
                         "2dd3b5c1-2a7f-4324-83d3-d05cff85676b",
                         "ba8b3250-6115-4680-8c21-6129c35ba738",
                         "bd067a4c-05b8-4def-925c-60ce83c01891",
                         "cd9e0589-187b-43c0-9692-50d8db7d0467",
                         "da36b9f6-5acd-4eae-af23-f49be4ed330e",
                         "dd915a74-f2d2-4888-a62d-4a5ebeae4e0e",
                         
                         "24168b97-a6b0-41cf-9836-29110572efc6",
                         "45cff797-6ed9-4a25-b659-30aa4caa46f2",
                         "72ea9678-7174-44c4-9535-de9cb733a15d",
                         "d3b6cb30-0a64-4f82-91ea-0bb14637ee17",
                         
                         "faf313a1-9ae4-43f4-bfc9-974281feac0e",
                         "da36b9f6-5acd-4eae-af23-f49be4ed330e",
                         "f3278405-0943-41a0-9f04-3a7733ec344d",
                         "2dd3b5c1-2a7f-4324-83d3-d05cff85676b",
                         "ba8b3250-6115-4680-8c21-6129c35ba738",
                         "bd067a4c-05b8-4def-925c-60ce83c01891",
                         "cd9e0589-187b-43c0-9692-50d8db7d0467",
                         "da36b9f6-5acd-4eae-af23-f49be4ed330e",
                         "dd915a74-f2d2-4888-a62d-4a5ebeae4e0e",
                         
                         "c6bbb6ef-ad16-4f3c-99e2-f693760173e0",
                         "acf7b380-2e9d-42d2-b1f8-a802a5e9d90a",
                         "a0ee69cf-df3b-44d7-8b14-19c42a1532a0",
                         "857bce66-f762-11e1-a439-00145eb45e9a",
                         "75214a8e-445a-4785-857e-5d1e06b9150d",
                         "4921f3f7-0991-49f4-8f56-91e7b9a55bf3",
                         "f3278405-0943-41a0-9f04-3a7733ec344d",
                         "857bce66-f762-11e1-a439-00145eb45e9a",
                         "37f56892-bd85-479a-803b-7fbfdfd6886c")) %>% 
    filter(!(datasetKey == "8a863029-f435-446a-821e-275f4f641165" & countryCod == "NL")) %>% 
    
      
  #TODO: Code remnant: Baptiste I don't know what this is for; columns can be selected with:
  # select(columnName, columnName, columnName) %>% 
  # it improves reproducibility if the column order changes
  spOcc = spocc[,c(1:2,4,6:11)]

  # Sanity check: 
  dim(spOcc)[1]
  #TODO replace the line above by the line below 
  # if(dim(spOcc)[1] == **replace this with the number of rows**) { print("Aves dataset ready.")} else {print("Error during Aves dataset setup. Try again")}
  
  
} else {print("Sorry, your query did not match any of the taxa offered by SynAnthrop. Try again.")}


# Homogenize column names with expectations of SynAnthrop function
spOcc %<>% rename(Species = species,
                 Abundance = individual,
                 Year = year)

# Species with very few individuals cannot be evaluated statistically with SynAnthrop.
# Discard species with fewer than the threshold number of individuals in the dataset.
# Calculate the total number of individuals found per species (sum of "Abundance")
threshold = 30 # we want at least *threshold* individuals per species
spOcc %<>% 
  group_by(Species) %>% #group per species
  filter(sum(Abundance) >= threshold)  %>% #keep species with >= 30 individuals
  ungroup()






## Load anthropization data #####
# load the raster file
ras_raw <- raster::raster(file.path(here(), "data", "HFP_2020_europe2b.tif"))

# create an object to visualize the raster
rast_to_plot <- terra::rast(ras_raw)
terra::crs(rast_to_plot) <- "EPSG:4326" # select the correct crs
rm(ras_raw) # remove the original object




#### Faire tourner la fonction SynAnthrop

##### Fonction manuelle


 
r = rast_to_plot
data = spOcc
value = resolution = 10
sim = 100
threshold = 100

# Create objects to store the results
speciesScores <- NULL # store species scores
effSizes <- NULL # store all effect sizes calculated
samplesList <-  NULL # store all samples drawn and observed


 


  # STEP 1bis (changement echelle) | Naturalness raster file: Aggregate cells at desired resolution ####
 
# STEP 1bis | Naturalness raster file: Aggregate cells at desired resolution ####
value=resolution
# raw raster files may be too large for the intended analysis. Cells are aggregated
# to decrease the file complexity.
cat(paste(Sys.time(), "- Aggregate naturalness raster cells at resolution", value, 
          "(this step may take a few minutes)\n"))
# aggregate raster cells at defined resolution



#### Methode avec 1 SpatVector avec différents polygones emprises :

# scale=vect("E:/Region_FR.shp")
#   names=scale$DREG_L_LIB # nom des polygones

# scale=NULL
# names="Europe"

# scale=vect("/home/genouest/inra_umr0985/bbongibault/Europe_countries.shp")
# scale=vect("E:/Synanthrope/Europe_countries.shp")
# names=scale$NAME_FREN # nom des polygones

# scale=vect("D:/QGIS/Biogeo/BiogeoRegions2016.shp")
#   names=scale$short_name # nom des polygones

if(is.null(scale)==F){    
  
  
  
  if(class(scale)=="SpatVector"){
    
    
    results1 <- lapply(1:length(scale), function(nb_scale){
      # results1 <- future_lapply(1:length(scale), function(nb_scale){
        
          # results <- future_lapply(poly_list, function(nb_scale){
      
      # nb_scale=nb_scale+1
      
      cat(paste('Synanthrop processing for "',names[nb_scale], '" scale\n\n'))
      
      if (exists("r_repro")==F){
        r_repro <- terra::rast("/home/genouest/inra_umr0985/bbongibault/r_repro.tif")
        # r_repro <- terra::rast("E:/Synanthrope/r_repro.tif")
        # r_repro <- terra::project(r, "EPSG:3035")
      }else{}
      
      scale2=vect("/home/genouest/inra_umr0985/bbongibault/Europe_countries.shp")
      # scale2=vect("E:/Synanthrope/Europe_countries.shp")
      scale_curr=scale2[nb_scale]
      scale_curr<- terra::project(scale_curr, crs(r_repro))
      
      
      r_crop = crop(r_repro, scale_curr)
      r_crop = mask(r_crop, scale_curr)
      # plot(r_crop)
      
      cat(paste(Sys.time(), "- Aggregate naturalness raster cells at resolution", value, 
                "(this step may take a few minutes)\n\n"))
      
      ras_reproj <- terra::aggregate(r_crop, fact = value, fun = mean, na.rm = T)
      
      rm(r_crop)
      rm(scale_curr)
      rm(scale2)
      
      ras <- terra::project(ras_reproj, "EPSG:4326")
      names(ras)="value"
      
      
      # Extract naturalness raster values associated with each cell
      ras.value <- data.frame(value = terra::values(ras)) 
      ras.value %<>% 
        mutate(Cell = as.numeric(rownames(ras.value))) %>% 
        # remove cells without naturalness values (e.g. offshore cells)
        na.omit() 


      # STEP 2 | Dataset: calculate sampling effort per cell ####
      
      # we calculate the number of times each cell was visited per year and in total
      # in order to obtain an index of sampling effort. Sampling effort is linked with
      # species detection, i.e., cells that were heavily prospected have very reliable
      # species lists. 
      
      # aggregate observations in cells at raster resolution based on XY coordinates
      data_res <- data # create a dataset from this resolution
      data_res %<>% mutate(Cell = terra::cellFromXY(ras, data_res %>% select(X, Y)))
      
      nb_obs_1=nrow(data_res)
      
            
      ##################
      data_res=na.omit(data_res)
      
      if (nrow(data_res)==0) {
        cat(paste0("No observations for this geographical area: ",names[nb_scale],". Please verify that there are indeed no observations for this area.\n\n"))    
        return(NULL)
      }
      
      data_resb=vect(data_res, geom = c("X", "Y"), crs = "EPSG:4326")
      
      data_resb=terra::project(data_resb, crs(ras_reproj))
      
      data_resb=terra::extract(ras_reproj, data_resb, cells=F, ID=F, bind=T)
      
      
      
      
      # plot(ras_reproj)
      # points(data_resb)
      
      data_res <- terra::project(data_resb, "EPSG:4326")
      
      rm(data_resb)
      
      data_res=as.data.frame(data_res, geom="XY")
      data_res=na.omit(data_res)
      
      if (nrow(data_res)==0) {
        cat(paste0("No observations for this geographical area: ",names[nb_scale],". Please verify that there are indeed no observations for this area.\n\n"))    
        return(NULL)
      }
      
      colnames(data_res)[colnames(data_res) == 'x'] <- 'X'
      colnames(data_res)[colnames(data_res) == 'y'] <- 'Y'
      
      
      nb_obs_2=nrow(data_res)
      
      cat(paste(nb_obs_1-nb_obs_2, "observation(s) were eliminated because of mismatch with the raster file\n\n"))


      # calculate the number of years each cell has been visited
      x_visits <- data_res %>%
        group_by(Cell) %>%
        summarize(nYearVisited = n_distinct(Year))
      
      # compile a data.frame with the number of visits for all raster cells,
      # including non-visited cells (0)
      x_visits_all <- data.frame(Cell = seq(1:terra::ncell(ras))) %>%
        dplyr::full_join(x_visits, by = "Cell") %>%
        mutate(nYearVisited = ifelse(is.na(nYearVisited), 0, nYearVisited))
      
      rm(x_visits)
      
      # obtain XY coordinates of raster cells
      cellCoord <- data.frame(terra::xyFromCell(ras, x_visits_all$Cell))
      cellCoord %<>% mutate(Cell = as.numeric(row.names(cellCoord)))
      x_visits_all %<>% left_join(cellCoord, by = "Cell")
      
      # interpolate sampling effort by kernel density
      kernel_sampling <- ks::kde(x = x_visits_all %>% select(x, y),
                                 w = x_visits_all$nYearVisited)
      
      
      
      # STEP 3alpha | Dataset: creation of artificial species ####
      
   
      #Calcul des quantiles pour créer une espèce artificielle synanthrope et une espèce artificielle anthropophobe afin de fixer les valeurs extrêmes de l'index        
      mini=0.15
      maxi=0.85
      
      q10=quantile(ras.value$value, probs=mini)
      q90=quantile(ras.value$value, probs=maxi)
      
      ## Ne garder que les mailles les plus urbanisées ou les plus préservées
      ## while permet de passer le threshold pour les deux espèces artificielles
      cell_ant=subset(ras.value, value<q10)
      cell_syn=subset(ras.value, value>q90)
      
      while (nrow(cell_ant)<100) {
        mini=mini+0.05
        q10=quantile(ras.value$value, probs=mini)
        cell_ant=subset(ras.value, value<q10)
        
      }
      print(paste0("Quantile for Anthropophe species is ", mini*100,"%"))
      
      
      while (nrow(cell_syn)<100) {
        maxi=maxi-0.05
        q90=quantile(ras.value$value, probs=maxi)
        cell_syn=subset(ras.value, value>q90)
        
      }
      print(paste0("Quantile for Synanthrope species is ", maxi*100,"%"))
      
      
      ## Récupérer coordonnées des mailles avec une valeur d'anthropisation
      coordonnees=terra::extract(ras, y=c(1:(nrow(ras)*ncol(ras))), xy=T)
      coordonnees%<>% 
        mutate(Cell = as.numeric(rownames(coordonnees))) %>% 
        # remove cells without naturalness values (e.g. offshore cells)
        na.omit()  
      
      coordonnees=coordonnees[,c(1,2,4)]
      
      ## Le nombre de maille avec une observation sera la moyenne du nombre de maille des vraies espèces      
      data_res2 = droplevels(data_res)      
      abund = round(mean(summary(data_res2$Species)))
      
      if(abund<threshold){ #to force number of artificial species occurrences to be higher than thresold
        abund=threshold+1
      }
      
      if(abund>nrow(cell_ant)){ #to reduce abund to priorize the number of cell by quantile if more than 100
        abund=nrow(cell_ant)
      }
      
      
      ## Tirer un nombre de maille dans lequel les espèces articielles seront dites "présentes"
      
      E_syn=cell_syn[sample(1:nrow(cell_syn),abund),]             
      E_syn=merge(E_syn, coordonnees, by="Cell")
      E_syn$Species="Synanthrope species"
      # E_syn$SumAbundance=1
      E_syn$month=sample(1:12,nrow(E_syn), replace=T)
      E_syn$Year=sample(2015:2025,nrow(E_syn), replace=T)
      E_syn$coordinate=1
      E_syn$Abundance=1
      E_syn$countryCod=data_res2$countryCod[1]
      E_syn$datasetKey="NONUSED"
      colnames(E_syn)[colnames(E_syn) == 'value'] <- 'HFP_2020_europe2b'
      colnames(E_syn)[colnames(E_syn) == 'x'] <- 'X'
      colnames(E_syn)[colnames(E_syn) == 'y'] <- 'Y'
      
      E_ant=cell_ant[sample(1:nrow(cell_ant),abund),]             
      E_ant=merge(E_ant, coordonnees, by="Cell")
      E_ant$Species="Anthropophobe species"
      # E_ant$SumAbundance=1
      E_ant$month=sample(1:12,nrow(E_ant), replace=T)
      E_ant$Year=sample(2015:2025,nrow(E_ant), replace=T)
      E_ant$coordinate=1
      E_ant$Abundance=1
      E_ant$countryCod=data_res2$countryCod[1]
      E_ant$datasetKey="NONUSED"
      colnames(E_ant)[colnames(E_ant) == 'value'] <- 'HFP_2020_europe2b'
      colnames(E_ant)[colnames(E_ant) == 'x'] <- 'X'
      colnames(E_ant)[colnames(E_ant) == 'y'] <- 'Y'
      
      # Joindre les observations artificielles au jeu de données d'occurrences
      data_res=rbind(data_res, E_ant)
      data_res=rbind(data_res, E_syn)
      
      rm(E_ant)
      rm(E_syn)
      rm(cell_ant)
      rm(cell_syn)
      
      # STEP 3 | Dataset: identify species to evaluate ####
      
      # Species seldom detected will produce unreliable synanthropy scores. 
      # We will only evaluate the synanthropy score of species that have been detected 
      # in at least 'threshold' cells.
      
      # sum the abundances per species per cell
      data_res %<>%
        group_by(Cell, Species) %>%
        summarize(SumAbundance = sum(Abundance))
      
      # add the XY coordinates of the cell (from the raster, not the original XY coordinates)
      data_res[,c("x","y")] <- terra::xyFromCell(ras,data_res$Cell)
      
      # count the number of cells in which each species was detected
      spDetection <- data_res %>% 
        group_by(Species) %>%
        summarize(nCellsPresent = n_distinct(Cell))
      
      # species will only be evaluated if they have been detected in more than 
      # 'threshold' cells
      spDetection %<>% mutate(evaluation = ifelse(nCellsPresent > threshold, 
                                                  "evaluated", "not evaluated"))
      spEvaluated <- spDetection %>% filter(evaluation == "evaluated")
      
      # if any species has been detected (more than 2 because 2 artificials species in any case)
if (nrow(spEvaluated)<=2) {
cat("Any species can be evaluated (see threshold cells)\n\n")    
  return(NULL)
}
      
      # if (nrow(spEvaluated)==0) {
      #   cat("Any species can be evaluated (see threshold cells)\n\n")    
      #   next
      # }
      
      spEvaluated %<>% mutate(spNum = seq(1:nrow(spEvaluated)))
      
      # subset the dataset to keep only species that will be evaluated
      x_evaluated <- data_res %>% filter(Species %in% unique(spEvaluated$Species)) %>% 
        mutate(variable = "Observed")
      
      cat("Species that will be evaluated are:\n")
      cat(spEvaluated$Species, sep = "\n")
      
      
      
      # STEP 4 | Generate simulated datasets (null distributions) ####
      
      cat(paste(Sys.time(), "- Generating null distributions for:\n"))
      
      # Now we can generate a null distribution of the species repartition by randomly
      # resampling cells from the convex hull of the observed species distribution. 
      # The result is a collection of cells that represents what the species 
      # distribution could be if it was random and not linked to naturalness factors.
      
      
      # create an empty object to store all the null distributions
      nullFull <- NULL
      # sp=unique(spEvaluated$Species)[1]
      
      # loop over the species to evaluate
      for(sp in unique(spEvaluated$Species)) {
        
        # prompt species evaluated
        cat(paste(Sys.time(), sp))
        cat(paste0(" (", spEvaluated %>% filter(Species == sp) %>% 
                     select(spNum), "/", max(spEvaluated$spNum), ")\n"))
        
        
        # we draw a convex hull of the species distribution and extract the coordinates
        # of all the cells inside this hull. This list of cells will then be randomly
        # resampled to obtain null distributions of the species distribution.
        
        # create a convex hull of the species distribution
        sp_points <- sf::st_as_sf(x_evaluated %>% filter(Species == sp), 
                                  coords = c("x", "y"), crs = 4326)
        
        
        
        # Comme les buffers en degrés sont imprécis, on transforme en projection métrique
        # EPSG 3035 : Lambert Europe
        sp_points_m <- st_transform(sp_points, 3035)
        
        # Créer des buffers de 20 km
        buffers <- st_buffer(sp_points_m, 20000)  # 10 km = 10000 m
        
        rm(sp_points_m)
        # Fusionner tous les buffers en un seul objet (MULTIPOLYGON)
        buffers_union <- st_union(buffers)
        
        # S'assurer que la géométrie est valide
        buffers_union <- st_make_valid(buffers_union)
        
        # Extraire uniquement les polygones
        buffers_poly <- st_collection_extract(buffers_union, "POLYGON")
        
        rm(buffers_union)
        # Reprojeter en WGS84 pour avoir lon/lat
        buffers_poly <- st_transform(buffers_poly, 4326)
        
        # Visualiser le résultat (facultatif)
        # plot(buffers_poly, border = "blue")
        # plot(sp_points, add = TRUE, col = "red", pch = 16)
        # plot(sp_resampled, add = TRUE, col = "green", pch = 16)
        
        
        # ggplot() +
        #   # polygone pays
        #   geom_sf(data = country_map_metro, fill = "lightblue", color = "black") +
        # 
        #   # polygones tampon
        #   geom_sf(data = buffers_poly,  fill = "darkgreen") +
        # 
        #   # points
        #   geom_sf(data = sp_points, color = "red", size = 2) +
        # 
        #   # thème
        #   theme_bw()
        
        
        
        
        # attribute raster cell numbers to cells of the convex hull
        # hullCoord <- tabularaster::cellnumbers(raster::raster(ras), buffers_poly) %>%
        #   rename(Cell = "cell_")
        Cell <- terra::cells(ras, vect(buffers_poly))[,2] #%>% rename(cell = "cell_")
        hullCoord <- as.data.frame(Cell) #######################################################################
        
        # convert cell numbers to XY coordinates
        hullCoord[,c("x","y")] <- terra::xyFromCell(ras, hullCoord$Cell)
        
        # select cells within the convex hull that are in the naturalness dataset 
        # (this excludes offshore points in cases where the convex hull includes marine
        # areas) 
        hullCoord %<>% filter(Cell %in% ras.value$Cell)
        
        # select kernel weights (sampling effort for each cell within the convex hull)
        kernel_weights <- data.frame(weight = kernel_sampling$w,
                                     Cell = seq(1:length(kernel_sampling$w))) %>%
          filter(Cell %in% hullCoord$Cell)
        
        
        
        
        # create an empty object to store the null distribution for this species
        nullSp <- NULL
        
        # loop over all the simulations requested
        for(i in 1:sim) {
          
          # prompt the simulation progress
          if (i == 1 & i == sim) {
            cat(paste0("Simulation ", as.numeric(i), "/", sim, "\n"))
          } else if(i == 1){
            cat(paste0("Simulation ", as.numeric(i), "/", sim, "... "))
          } else if (i == sim) {
            cat(paste0(" ", as.numeric(i), "/", sim, "\n"))
          } else if (i %% 100 == 0) {
            cat(paste0(" ", as.numeric(i), "/", sim, "... "))
          } 
          
          
          # resample cells within the convex hull, with their associated kernel weight
          # sp_resampled <- sample_n(kernel_weights, #x_visits_all,
          #                          size = spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent),
          #                          replace = FALSE,
          #                          weight = weight) # weight = nYearVisited

          
##Probleme kernel pour espèces artificielles####          

if((sum(kernel_weights$weight > 0) < spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent))& (sp=="Anthropophobe species" | sp=="Synanthrope species")){         
          size_requested <- spEvaluated %>%
            filter(Species == sp) %>%
            pull(nCellsPresent)
          
kernel_weights_pos <- kernel_weights %>%
            filter(weight > 0)
          
size_final <- min(size_requested, nrow(kernel_weights_pos))
          
sp_resampled <- sample_n(kernel_weights_pos,
                         size = size_final,
                         replace = FALSE,
                         weight = weight)
          
}else{ 
  sp_resampled <- sample_n(kernel_weights, #x_visits_all,
                           size = spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent),
                           replace = FALSE,
                           weight = weight) # weight = nYearVisited
  
} 
          
####Fin de la modif pour l'erreur 
          
          
          # assign simulation number and species
          sp_resampled %<>% mutate(simulation = i,
                                   Species = sp,
                                   variable = "Null")
          # add sampled data from this simulation to the dataset of this species
          nullSp <- rbind(nullSp, sp_resampled)
          
        } # end of simulations for this species
        
        # concatenate sampled data from all species
        nullFull <- rbind(nullSp, nullFull)
        
      } # end of the loop for all species 
      
      
      # add XY coordinates to null distributions
      nullFull[,c("x","y")] <- terra::xyFromCell(ras, nullFull$Cell)
      
      # bind the null and observed datasets
      datasetFinal <- rbind(nullFull %>% select(-weight),#-nYearVisited),  
                            x_evaluated %>% select(-SumAbundance) %>% mutate(simulation = "no"))
      datasetFinal %<>% mutate(Resolution = value) 
      
      datasetFinal %<>% mutate(Scale = names[nb_scale]) ########################################################################
      
      # STEP 5 | Calculate effect size #### 
      
      # now that we have the observed and simulated distributions, we are going to
      # calculate the effect size of their difference in naturalness values
      
      cat(paste(Sys.time(), "- Calculating effect sizes...\n"))
      
      # create an empty object to store all effect sizes for all species
      effSizesFull <- NULL
      
      # for every simulation, calculate the effect size
      for(run in unique(nullFull$simulation)) {
        
        # prompt the calculation progress
        if (run == 1 & run == max(nullFull$simulation)) {
          cat(paste0("Simulation ", as.numeric(run), "/", max(nullFull$simulation), "\n"))
        } else if(run == 1){
          cat(paste0("Simulation ", as.numeric(run), "/", max(nullFull$simulation), "... "))
        } else if (run == sim) {
          cat(paste0(" ", as.numeric(run), "/", max(nullFull$simulation), "\n"))
        } else if (run %% 100 == 0) {
          cat(paste0(" ", as.numeric(run), "/", max(nullFull$simulation), "... "))
        } 
        
        
        # select a single simulation and the observed data, add naturalness raster values
        runN <- datasetFinal %>% dplyr::filter(simulation == run |
                                                 simulation == "no") %>% 
          inner_join(ras.value, by = "Cell") 
        # TODO: note that some points don't have naturalness values here!
        # TODO: this ends up with different sample sizes for the simulation and the observed
        # data in some cases, visible in the third result data frame.
        
        # create an empty object to store the species effect sizes
        effSizesSp <- NULL
        
        for(sp in unique(runN$Species)) {
          
          # select one species
          runNSp <- runN %>% filter(Species == sp)
          # calculate the effect size: (simulated mean - observed mean)/estimated sd
          effSizesrunNSp <- runNSp %>% 
            rstatix::cohens_d(value ~ variable) %>% # equal = TRUE ?
            mutate(Species = sp,
                   Run = run,
                   Resolution = value,
                   Scale=names[nb_scale])###########################################################################################
          # add effect size to the list
          effSizesSp <- rbind(effSizesSp, effSizesrunNSp)
          
        }
        
        # add species effect sizes to the list
        effSizesFull <- rbind(effSizesFull, effSizesSp)
        
      } #end of the loop by run
      
      
      # STEP 6 | Scoring ####
      
      # synthesize effect sizes per species
      effSizesFull_summary <- effSizesFull %>%
        group_by(Species) %>%
        summarise(mean = mean(effsize), 
                  nRun = n(),
                  Scale=Scale[1])
      
      # rescale effect sizes to obtain the synanthropy score per species
      effSizesFull_summary %<>% 
        mutate(Index = round(scales::rescale(mean, to = c(10, 1))),
               Resolution = value)
      
      
      # STEP 7 | Compile results for this resolution ####
      
      #1: species scores
      speciesScores <- rbind(speciesScores, effSizesFull_summary)
      
      #2: simulated and observed datasets
      samplesList <- rbind(samplesList, datasetFinal)
      
      #3: all effect sizes
      effSizes <- rbind(effSizes, effSizesFull)
      
      # create a list for all 3 tables
      results1 <- list("speciesScores" = speciesScores, "effSizes" = effSizes, "samplesList" = samplesList) 
      
      cat(paste(Sys.time(), "Analysis finished for resolution", value, "\n\n"))
      # end of resolution loop
      
      cat(paste(Sys.time(), "All done.\n"))    
      return(results1)
    })
    # return(results1)
    
  }else{
    
    
    for(nb_scale in 1:length(scale)){
      
      if (exists("r_repro")==F){
        r_repro <- terra::project(r, "EPSG:3035")
      }else{}
      
      
      
      cat(paste('Synanthrop processing for "',names(scale)[[nb_scale]], '" scale\n\n'))
      
      
      scale_curr=scale[[nb_scale]]
      scale_curr<- terra::project(scale_curr, crs(r_repro))
      
      
      r_crop = crop(r_repro, scale_curr)
      r_crop = mask(r_crop, scale_curr)
      # plot(r_crop)
      
      cat(paste(Sys.time(), "- Aggregate naturalness raster cells at resolution", value, 
                "(this step may take a few minutes)\n\n"))
      
      ras_reproj <- terra::aggregate(r_crop, fact = value, fun = mean, na.rm = T)
      
      ras <- terra::project(ras_reproj, "EPSG:4326")
      names(ras)="value"
      
      # Extract naturalness raster values associated with each cell
      ras.value <- data.frame(value = terra::values(ras)) 
      ras.value %<>% 
        mutate(Cell = as.numeric(rownames(ras.value))) %>% 
        # remove cells without naturalness values (e.g. offshore cells)
        na.omit() 
      
      # STEP 2 | Dataset: calculate sampling effort per cell ####
      
      # we calculate the number of times each cell was visited per year and in total
      # in order to obtain an index of sampling effort. Sampling effort is linked with
      # species detection, i.e., cells that were heavily prospected have very reliable
      # species lists. 
      
      # aggregate observations in cells at raster resolution based on XY coordinates
      data_res <- data # create a dataset from this resolution
      data_res %<>% mutate(Cell = terra::cellFromXY(ras, data_res %>% select(X, Y)))
      
      nb_obs_1=nrow(data_res)
      ##################
      data_res=na.omit(data_res)
      
      
      if (nrow(data_res)==0) {
        cat(paste0("No observations for this geographical area: ",names(scale)[[nb_scale]],". Please verify that there are indeed no observations for this area.\n\n"))    
        return(NULL)
      }
      
      data_resb=vect(data_res, geom = c("X", "Y"), crs = "EPSG:4326")
      
      data_resb=terra::project(data_resb, crs(ras_reproj))
      
      data_resb=terra::extract(ras_reproj, data_resb, cells=F, ID=F, bind=T)
      
      
      
      
      # plot(ras_reproj)
      # points(data_resb)
      
      data_res <- terra::project(data_resb, "EPSG:4326")
      
      data_res=as.data.frame(data_res, geom="XY")
      data_res=na.omit(data_res)
      
      if (nrow(data_res)==0) {
        cat(paste0("No observations for this geographical area: ",names(scale)[[nb_scale]],". Please verify that there are indeed no observations for this area.\n\n"))    
        return(NULL)
      }
      
      colnames(data_res)[colnames(data_res) == 'x'] <- 'X'
      colnames(data_res)[colnames(data_res) == 'y'] <- 'Y'
      
      
      ##################
      
      
      # missingcoord <- dim(data_res %>% filter(is.na(Cell)))[1]
      # data_res %<>% filter(!is.na(Cell))
      
      nb_obs_2=nrow(data_res)
      
      cat(paste(nb_obs_1-nb_obs_2, "observation(s) were eliminated because of mismatch with the raster file\n\n"))
      
      # calculate the number of years each cell has been visited
      x_visits <- data_res %>%
        group_by(Cell) %>%
        summarize(nYearVisited = n_distinct(Year))
      
      # compile a data.frame with the number of visits for all raster cells,
      # including non-visited cells (0)
      x_visits_all <- data.frame(Cell = seq(1:terra::ncell(ras))) %>%
        dplyr::full_join(x_visits, by = "Cell") %>%
        mutate(nYearVisited = ifelse(is.na(nYearVisited), 0, nYearVisited))
      
      # obtain XY coordinates of raster cells
      cellCoord <- data.frame(terra::xyFromCell(ras, x_visits_all$Cell))
      cellCoord %<>% mutate(Cell = as.numeric(row.names(cellCoord)))
      x_visits_all %<>% left_join(cellCoord, by = "Cell")
      
      # interpolate sampling effort by kernel density
      kernel_sampling <- ks::kde(x = x_visits_all %>% select(x, y),
                                 w = x_visits_all$nYearVisited)
      
      
      
      # STEP 3alpha | Dataset: creation of artificial species ####
      
      
      #Calcul des quantiles pour créer une espèce artificielle synanthrope et une espèce artificielle anthropophobe afin de fixer les valeurs extrêmes de l'index        
      mini=0.15
      maxi=0.85
      
      q10=quantile(ras.value$value, probs=mini)
      q90=quantile(ras.value$value, probs=maxi)
      
      ## Ne garder que les mailles les plus urbanisées ou les plus préservées
      ## while permet de passer le threshold pour les deux espèces artificielles
      cell_ant=subset(ras.value, value<q10)
      cell_syn=subset(ras.value, value>q90)
      
      while (nrow(cell_ant)<100) {
        mini=mini+0.05
        q10=quantile(ras.value$value, probs=mini)
        cell_ant=subset(ras.value, value<q10)
        
      }
      print(paste0("Quantile for Anthropophe species is ", mini*100,"%"))
      
      
      while (nrow(cell_syn)<100) {
        maxi=maxi-0.05
        q90=quantile(ras.value$value, probs=maxi)
        cell_syn=subset(ras.value, value>q90)
        
      }
      print(paste0("Quantile for Synanthrope species is ", maxi*100,"%"))
      
      
      ## Récupérer coordonnées des mailles avec une valeur d'anthropisation
      coordonnees=terra::extract(ras, y=c(1:(nrow(ras)*ncol(ras))), xy=T)
      coordonnees%<>% 
        mutate(Cell = as.numeric(rownames(coordonnees))) %>% 
        # remove cells without naturalness values (e.g. offshore cells)
        na.omit()  
      
      coordonnees=coordonnees[,c(1,2,4)]
      
      ## Le nombre de maille avec une observation sera la moyenne du nombre de maille des vraies espèces      
      data_res2 = droplevels(data_res)      
      abund = round(mean(summary(data_res2$Species)))
      
      if(abund<threshold){ #to force number of artificial species occurrences to be higher than thresold
        abund=threshold+1
      }
      
      if(abund>nrow(cell_ant)){ #to reduce abund to priorize the number of cell by quantile if more than 100
        abund=nrow(cell_ant)
      }
      
      
      ## Tirer un nombre de maille dans lequel les espèces articielles seront dites "présentes"
      
      E_syn=cell_syn[sample(1:nrow(cell_syn),abund),]             
      E_syn=merge(E_syn, coordonnees, by="Cell")
      E_syn$Species="Synanthrope species"
      # E_syn$SumAbundance=1
      E_syn$month=sample(1:12,nrow(E_syn), replace=T)
      E_syn$Year=sample(2015:2025,nrow(E_syn), replace=T)
      E_syn$coordinate=1
      E_syn$Abundance=1
      E_syn$countryCod=data_res2$countryCod[1]
      E_syn$datasetKey="NONUSED"
      colnames(E_syn)[colnames(E_syn) == 'value'] <- 'HFP_2020_europe2b'
      colnames(E_syn)[colnames(E_syn) == 'x'] <- 'X'
      colnames(E_syn)[colnames(E_syn) == 'y'] <- 'Y'
      
      E_ant=cell_ant[sample(1:nrow(cell_ant),abund),]             
      E_ant=merge(E_ant, coordonnees, by="Cell")
      E_ant$Species="Anthropophobe species"
      # E_ant$SumAbundance=1
      E_ant$month=sample(1:12,nrow(E_ant), replace=T)
      E_ant$Year=sample(2015:2025,nrow(E_ant), replace=T)
      E_ant$coordinate=1
      E_ant$Abundance=1
      E_ant$countryCod=data_res2$countryCod[1]
      E_ant$datasetKey="NONUSED"
      colnames(E_ant)[colnames(E_ant) == 'value'] <- 'HFP_2020_europe2b'
      colnames(E_ant)[colnames(E_ant) == 'x'] <- 'X'
      colnames(E_ant)[colnames(E_ant) == 'y'] <- 'Y'
      
      # Joindre les observations artificielles au jeu de données d'occurrences
      data_res=rbind(data_res, E_ant)
      data_res=rbind(data_res, E_syn)
      
      rm(E_ant)
      rm(E_syn)
      rm(cell_ant)
      rm(cell_syn)
      
      
      # STEP 3 | Dataset: identify species to evaluate ####
      
      # Species seldom detected will produce unreliable synanthropy scores. 
      # We will only evaluate the synanthropy score of species that have been detected 
      # in at least 'threshold' cells.
      
      # sum the abundances per species per cell
      data_res %<>%
        group_by(Cell, Species) %>%
        summarize(SumAbundance = sum(Abundance))
      
      # add the XY coordinates of the cell (from the raster, not the original XY coordinates)
      data_res[,c("x","y")] <- terra::xyFromCell(ras,data_res$Cell)
      
      # count the number of cells in which each species was detected
      spDetection <- data_res %>% 
        group_by(Species) %>%
        summarize(nCellsPresent = n_distinct(Cell))
      
      # species will only be evaluated if they have been detected in more than 
      # 'threshold' cells
      spDetection %<>% mutate(evaluation = ifelse(nCellsPresent > threshold, 
                                                  "evaluated", "not evaluated"))
      spEvaluated <- spDetection %>% filter(evaluation == "evaluated")
      
      # if any species has been detected
      if (nrow(spEvaluated)==0) {
        cat("Any species can be evaluated (see threshold cells\n\n")    
        return(NULL)
      }
      
      
      spEvaluated %<>% mutate(spNum = seq(1:nrow(spEvaluated)))
      
      # subset the dataset to keep only species that will be evaluated
      x_evaluated <- data_res %>% filter(Species %in% unique(spEvaluated$Species)) %>% 
        mutate(variable = "Observed")
      
      cat("Species that will be evaluated are:\n")
      cat(spEvaluated$Species, sep = "\n")
      
      
      # STEP 4 | Generate simulated datasets (null distributions) ####
      
      cat(paste(Sys.time(), "- Generating null distributions for:\n"))
      
      # Now we can generate a null distribution of the species repartition by randomly
      # resampling cells from the convex hull of the observed species distribution. 
      # The result is a collection of cells that represents what the species 
      # distribution could be if it was random and not linked to naturalness factors.
      
      
      # create an empty object to store all the null distributions
      nullFull <- NULL
      # sp=unique(spEvaluated$Species)[1]
      
      # loop over the species to evaluate
      for(sp in unique(spEvaluated$Species)) {
        
        # prompt species evaluated
        cat(paste(Sys.time(), sp))
        cat(paste0(" (", spEvaluated %>% filter(Species == sp) %>% 
                     select(spNum), "/", max(spEvaluated$spNum), ")\n"))
        
        
        # we draw a convex hull of the species distribution and extract the coordinates
        # of all the cells inside this hull. This list of cells will then be randomly
        # resampled to obtain null distributions of the species distribution.
        
        # create a convex hull of the species distribution
        sp_points <- sf::st_as_sf(x_evaluated %>% filter(Species == sp), 
                                  coords = c("x", "y"), crs = 4326)
        
        
        
        # Comme les buffers en degrés sont imprécis, on transforme en projection métrique
        # EPSG 3035 : Lambert Europe
        sp_points_m <- st_transform(sp_points, 3035)
        
        # Créer des buffers de 10 km
        buffers <- st_buffer(sp_points_m, 30000)  # 10 km = 10000 m
        
        # Fusionner tous les buffers en un seul objet (MULTIPOLYGON)
        buffers_union <- st_union(buffers)
        
        # S'assurer que la géométrie est valide
        buffers_union <- st_make_valid(buffers_union)
        
        # Extraire uniquement les polygones
        buffers_poly <- st_collection_extract(buffers_union, "POLYGON")
        
        # Reprojeter en WGS84 pour avoir lon/lat
        buffers_poly <- st_transform(buffers_poly, 4326)
        
        # Visualiser le résultat (facultatif)
        # plot(buffers_poly, border = "blue")
        # plot(sp_points, add = TRUE, col = "red", pch = 16)
        # plot(sp_resampled, add = TRUE, col = "green", pch = 16)
        
        
        # ggplot() +
        #   # polygone pays
        #   geom_sf(data = country_map_metro, fill = "lightblue", color = "black") +
        # 
        #   # polygones tampon
        #   geom_sf(data = buffers_poly,  fill = "darkgreen") +
        # 
        #   # points
        #   geom_sf(data = sp_points, color = "red", size = 2) +
        # 
        #   # thème
        #   theme_bw()
        
        
        
        
        # attribute raster cell numbers to cells of the convex hull
        # hullCoord <- tabularaster::cellnumbers(raster::raster(ras), buffers_poly) %>%
        #   rename(Cell = "cell_")
        Cell <- terra::cells(ras, vect(buffers_poly))[,2] #%>% rename(cell = "cell_")
        hullCoord <- as.data.frame(Cell) #######################################################################
        
        
        # convert cell numbers to XY coordinates
        hullCoord[,c("x","y")] <- terra::xyFromCell(ras, hullCoord$Cell)
        # select cells within the convex hull that are in the naturalness dataset 
        # (this excludes offshore points in cases where the convex hull includes marine
        # areas) 
        hullCoord %<>% filter(Cell %in% ras.value$Cell)
        
        # select kernel weights (sampling effort for each cell within the convex hull)
        kernel_weights <- data.frame(weight = kernel_sampling$w,
                                     Cell = seq(1:length(kernel_sampling$w))) %>%
          filter(Cell %in% hullCoord$Cell)
        
        
        
        
        # create an empty object to store the null distribution for this species
        nullSp <- NULL
        
        # loop over all the simulations requested
        for(i in 1:sim) {
          
          # prompt the simulation progress
          if (i == 1 & i == sim) {
            cat(paste0("Simulation ", as.numeric(i), "/", sim, "\n"))
          } else if(i == 1){
            cat(paste0("Simulation ", as.numeric(i), "/", sim, "... "))
          } else if (i == sim) {
            cat(paste0(" ", as.numeric(i), "/", sim, "\n"))
          } else if (i %% 100 == 0) {
            cat(paste0(" ", as.numeric(i), "/", sim, "... "))
          } 
          
          
          # resample cells within the convex hull, with their associated kernel weight
          sp_resampled <- sample_n(kernel_weights, #x_visits_all,
                                   size = spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent), 
                                   replace = FALSE,
                                   weight = weight) # weight = nYearVisited
          
          # assign simulation number and species
          sp_resampled %<>% mutate(simulation = i,
                                   Species = sp,
                                   variable = "Null")
          # add sampled data from this simulation to the dataset of this species
          nullSp <- rbind(nullSp, sp_resampled)
          
        } # end of simulations for this species
        
        # concatenate sampled data from all species
        nullFull <- rbind(nullSp, nullFull)
        
      } # end of the loop for all species 
      
      
      # add XY coordinates to null distributions
      nullFull[,c("x","y")] <- terra::xyFromCell(ras, nullFull$Cell)
      
      # bind the null and observed datasets
      datasetFinal <- rbind(nullFull %>% select(-weight),#-nYearVisited),  
                            x_evaluated %>% select(-SumAbundance) %>% mutate(simulation = "no"))
      datasetFinal %<>% mutate(Resolution = value) 
      
      datasetFinal %<>% mutate(Scale = names(scale)[[nb_scale]]) ########################################################################
      
      # STEP 5 | Calculate effect size #### 
      
      # now that we have the observed and simulated distributions, we are going to
      # calculate the effect size of their difference in naturalness values
      
      cat(paste(Sys.time(), "- Calculating effect sizes...\n"))
      
      # create an empty object to store all effect sizes for all species
      effSizesFull <- NULL
      
      # for every simulation, calculate the effect size
      for(run in unique(nullFull$simulation)) {
        
        # prompt the calculation progress
        if (run == 1 & run == max(nullFull$simulation)) {
          cat(paste0("Simulation ", as.numeric(run), "/", max(nullFull$simulation), "\n"))
        } else if(run == 1){
          cat(paste0("Simulation ", as.numeric(run), "/", max(nullFull$simulation), "... "))
        } else if (run == sim) {
          cat(paste0(" ", as.numeric(run), "/", max(nullFull$simulation), "\n"))
        } else if (run %% 100 == 0) {
          cat(paste0(" ", as.numeric(run), "/", max(nullFull$simulation), "... "))
        } 
        
        
        # select a single simulation and the observed data, add naturalness raster values
        runN <- datasetFinal %>% dplyr::filter(simulation == run |
                                                 simulation == "no") %>% 
          inner_join(ras.value, by = "Cell") 
        # TODO: note that some points don't have naturalness values here!
        # TODO: this ends up with different sample sizes for the simulation and the observed
        # data in some cases, visible in the third result data frame.
        
        # create an empty object to store the species effect sizes
        effSizesSp <- NULL
        
        for(sp in unique(runN$Species)) {
          
          # select one species
          runNSp <- runN %>% filter(Species == sp)
          # calculate the effect size: (simulated mean - observed mean)/estimated sd
          effSizesrunNSp <- runNSp %>% 
            rstatix::cohens_d(value ~ variable) %>% # equal = TRUE ?
            mutate(Species = sp,
                   Run = run,
                   Resolution = value,
                   Scale=names(scale)[[nb_scale]])###########################################################################################
          # add effect size to the list
          effSizesSp <- rbind(effSizesSp, effSizesrunNSp)
          
        }
        
        # add species effect sizes to the list
        effSizesFull <- rbind(effSizesFull, effSizesSp)
        
      } #end of the loop by run
      
      
      
      
      # STEP 6 | Scoring ####
      
      # synthesize effect sizes per species
      effSizesFull_summary <- effSizesFull %>%
        group_by(Species) %>%
        summarise(mean = mean(effsize), 
                  nRun = n(),
                  Scale=Scale[1])
      
      # rescale effect sizes to obtain the synanthropy score per species
      effSizesFull_summary %<>% 
        mutate(Index = round(scales::rescale(mean, to = c(10, 1))),
               Resolution = value)
      
      
      # STEP 7 | Compile results for this resolution ####
      
      #1: species scores
      speciesScores <- rbind(speciesScores, effSizesFull_summary)
      
      #2: simulated and observed datasets
      samplesList <- rbind(samplesList, datasetFinal)
      
      #3: all effect sizes
      effSizes <- rbind(effSizes, effSizesFull)
      
      # create a list for all 3 tables
      results <- list("speciesScores" = speciesScores, "effSizes" = effSizes, "samplesList" = samplesList) 
      
      cat(paste(Sys.time(), "Analysis finished for resolution", value, "\n\n"))
      # end of resolution loop
      
      cat(paste(Sys.time(), "All done.\n"))
      
    }
    
    
  }
  # return(results1)
}else{
  results1 <- lapply(1:1, function(nb_scale){
    # results1 <- future_lapply(1:1, function(nb_scale){
  if (exists("r_repro")==F){
    r_repro <- terra::rast("/home/genouest/inra_umr0985/bbongibault/r_repro.tif")
    # r_repro <- terra::rast("E:/Synanthrope/r_repro.tif")
    # r_repro <- terra::project(r, "EPSG:3035")
  }else{}
  
  scale_curr=vect("/home/genouest/inra_umr0985/bbongibault/Europe.shp")
  # scale_curr=vect("E:/Synanthrope/Europe.shp")
  scale_curr<- terra::project(scale_curr, crs(r_repro))
  
  
  r_crop = crop(r_repro, scale_curr)
  r_crop = mask(r_crop, scale_curr)
  # plot(r_crop)
  
  cat(paste(Sys.time(), "- Aggregate naturalness raster cells at resolution", value, 
            "(this step may take a few minutes)\n\n"))
  
  ras_reproj <- terra::aggregate(r_crop, fact = value, fun = mean, na.rm = T)
  
  rm(r_crop)
  rm(scale_curr)
  # rm(scale2)
  
  ras <- terra::project(ras_reproj, "EPSG:4326")
  names(ras)="value"
  
  ras <- terra::aggregate(r, fact = value, fun = mean, na.rm = T)
  
  # Extract naturalness raster values associated with each cell
  ras.value <- data.frame(value = terra::values(ras)) 
  ras.value %<>% 
    mutate(Cell = as.numeric(rownames(ras.value))) %>% 
    # remove cells without naturalness values (e.g. offshore cells)
    na.omit() 
  ras <- terra::project(ras_reproj, "EPSG:4326")
  names(ras)="value"
  
  
  # Extract naturalness raster values associated with each cell
  ras.value <- data.frame(value = terra::values(ras)) 
  ras.value %<>% 
    mutate(Cell = as.numeric(rownames(ras.value))) %>% 
    # remove cells without naturalness values (e.g. offshore cells)
    na.omit() 
  
  
  # STEP 2 | Dataset: calculate sampling effort per cell ####
  
  # we calculate the number of times each cell was visited per year and in total
  # in order to obtain an index of sampling effort. Sampling effort is linked with
  # species detection, i.e., cells that were heavily prospected have very reliable
  # species lists. 
  
  # aggregate observations in cells at raster resolution based on XY coordinates
  data_res <- data # create a dataset from this resolution
  data_res %<>% mutate(Cell = terra::cellFromXY(ras, data_res %>% select(X, Y)))
  
  nb_obs_1=nrow(data_res)
  
  
  ##################
  data_res=na.omit(data_res)
  
  if (nrow(data_res)==0) {
    cat(paste0("No observations for this geographical area. Please verify that there are indeed no observations for this area.\n\n"))    
    return(NULL)
  }
  
  data_resb=vect(data_res, geom = c("X", "Y"), crs = "EPSG:4326")
  
  data_resb=terra::project(data_resb, crs(ras_reproj))
  
  data_resb=terra::extract(ras_reproj, data_resb, cells=F, ID=F, bind=T)
  
  
  
  
  # plot(ras_reproj)
  # points(data_resb)
  
  data_res <- terra::project(data_resb, "EPSG:4326")
  
  rm(data_resb)
  
  data_res=as.data.frame(data_res, geom="XY")
  data_res=na.omit(data_res)
  
  if (nrow(data_res)==0) {
    cat(paste0("No observations for this geographical area. Please verify that there are indeed no observations for this area.\n\n"))    
    return(NULL)
  }
  
  colnames(data_res)[colnames(data_res) == 'x'] <- 'X'
  colnames(data_res)[colnames(data_res) == 'y'] <- 'Y'
  
  
  nb_obs_2=nrow(data_res)
  
  cat(paste(nb_obs_1-nb_obs_2, "observation(s) were eliminated because of mismatch with the raster file\n\n"))
  
  
  # calculate the number of years each cell has been visited
  x_visits <- data_res %>%
    group_by(Cell) %>%
    summarize(nYearVisited = n_distinct(Year))
  
  # compile a data.frame with the number of visits for all raster cells,
  # including non-visited cells (0)
  x_visits_all <- data.frame(Cell = seq(1:terra::ncell(ras))) %>%
    dplyr::full_join(x_visits, by = "Cell") %>%
    mutate(nYearVisited = ifelse(is.na(nYearVisited), 0, nYearVisited))
  
  rm(x_visits)
  
  # obtain XY coordinates of raster cells
  cellCoord <- data.frame(terra::xyFromCell(ras, x_visits_all$Cell))
  cellCoord %<>% mutate(Cell = as.numeric(row.names(cellCoord)))
  x_visits_all %<>% left_join(cellCoord, by = "Cell")
  
  # interpolate sampling effort by kernel density
  kernel_sampling <- ks::kde(x = x_visits_all %>% select(x, y),
                             w = x_visits_all$nYearVisited)
  
  
  
  # STEP 3alpha | Dataset: creation of artificial species ####
  
  
  #Calcul des quantiles pour créer une espèce artificielle synanthrope et une espèce artificielle anthropophobe afin de fixer les valeurs extrêmes de l'index        
  mini=0.15
  maxi=0.85
  
  q10=quantile(ras.value$value, probs=mini)
  q90=quantile(ras.value$value, probs=maxi)
  
  ## Ne garder que les mailles les plus urbanisées ou les plus préservées
  ## while permet de passer le threshold pour les deux espèces artificielles
  cell_ant=subset(ras.value, value<q10)
  cell_syn=subset(ras.value, value>q90)
  
  while (nrow(cell_ant)<100) {
    mini=mini+0.05
    q10=quantile(ras.value$value, probs=mini)
    cell_ant=subset(ras.value, value<q10)
    
  }
  print(paste0("Quantile for Anthropophe species is ", mini*100,"%"))
  
  
  while (nrow(cell_syn)<100) {
    maxi=maxi-0.05
    q90=quantile(ras.value$value, probs=maxi)
    cell_syn=subset(ras.value, value>q90)
    
  }
  print(paste0("Quantile for Synanthrope species is ", maxi*100,"%"))
  
  
  ## Récupérer coordonnées des mailles avec une valeur d'anthropisation
  coordonnees=terra::extract(ras, y=c(1:(nrow(ras)*ncol(ras))), xy=T)
  coordonnees%<>% 
    mutate(Cell = as.numeric(rownames(coordonnees))) %>% 
    # remove cells without naturalness values (e.g. offshore cells)
    na.omit()  
  
  coordonnees=coordonnees[,c(1,2,4)]
  
  ## Le nombre de maille avec une observation sera la moyenne du nombre de maille des vraies espèces      
  data_res2 = droplevels(data_res)      
  abund = round(mean(summary(data_res2$Species)))
  
  if(abund<threshold){ #to force number of artificial species occurrences to be higher than thresold
    abund=threshold+1
  }
  
  if(abund>nrow(cell_ant)){ #to reduce abund to priorize the number of cell by quantile if more than 100
    abund=nrow(cell_ant)
  }
  
  
  ## Tirer un nombre de maille dans lequel les espèces articielles seront dites "présentes"
  
  E_syn=cell_syn[sample(1:nrow(cell_syn),abund),]             
  E_syn=merge(E_syn, coordonnees, by="Cell")
  E_syn$Species="Synanthrope species"
  # E_syn$SumAbundance=1
  E_syn$month=sample(1:12,nrow(E_syn), replace=T)
  E_syn$Year=sample(2015:2025,nrow(E_syn), replace=T)
  E_syn$coordinate=1
  E_syn$Abundance=1
  E_syn$countryCod=data_res2$countryCod[1]
  E_syn$datasetKey="NONUSED"
  colnames(E_syn)[colnames(E_syn) == 'value'] <- 'HFP_2020_europe2b'
  colnames(E_syn)[colnames(E_syn) == 'x'] <- 'X'
  colnames(E_syn)[colnames(E_syn) == 'y'] <- 'Y'
  
  E_ant=cell_ant[sample(1:nrow(cell_ant),abund),]             
  E_ant=merge(E_ant, coordonnees, by="Cell")
  E_ant$Species="Anthropophobe species"
  # E_ant$SumAbundance=1
  E_ant$month=sample(1:12,nrow(E_ant), replace=T)
  E_ant$Year=sample(2015:2025,nrow(E_ant), replace=T)
  E_ant$coordinate=1
  E_ant$Abundance=1
  E_ant$countryCod=data_res2$countryCod[1]
  E_ant$datasetKey="NONUSED"
  colnames(E_ant)[colnames(E_ant) == 'value'] <- 'HFP_2020_europe2b'
  colnames(E_ant)[colnames(E_ant) == 'x'] <- 'X'
  colnames(E_ant)[colnames(E_ant) == 'y'] <- 'Y'
  
  # Joindre les observations artificielles au jeu de données d'occurrences
  data_res=rbind(data_res, E_ant)
  data_res=rbind(data_res, E_syn)
  
  rm(E_ant)
  rm(E_syn)
  rm(cell_ant)
  rm(cell_syn)
  
  # STEP 3 | Dataset: identify species to evaluate ####
  
  # Species seldom detected will produce unreliable synanthropy scores. 
  # We will only evaluate the synanthropy score of species that have been detected 
  # in at least 'threshold' cells.
  
  # sum the abundances per species per cell
  data_res %<>%
    group_by(Cell, Species) %>%
    summarize(SumAbundance = sum(Abundance))
  
  # add the XY coordinates of the cell (from the raster, not the original XY coordinates)
  data_res[,c("x","y")] <- terra::xyFromCell(ras,data_res$Cell)
  
  # count the number of cells in which each species was detected
  spDetection <- data_res %>% 
    group_by(Species) %>%
    summarize(nCellsPresent = n_distinct(Cell))
  
  # species will only be evaluated if they have been detected in more than 
  # 'threshold' cells
  spDetection %<>% mutate(evaluation = ifelse(nCellsPresent > threshold, 
                                              "evaluated", "not evaluated"))
  spEvaluated <- spDetection %>% filter(evaluation == "evaluated")
  
  # if any species has been detected (more than 2 because 2 artificials species in any case)
  if (nrow(spEvaluated)<=2) {
    cat("Any species can be evaluated (see threshold cells)\n\n")    
    return(NULL)
  }
  
  # if (nrow(spEvaluated)==0) {
  #   cat("Any species can be evaluated (see threshold cells)\n\n")    
  #   next
  # }
  
  spEvaluated %<>% mutate(spNum = seq(1:nrow(spEvaluated)))
  
  # subset the dataset to keep only species that will be evaluated
  x_evaluated <- data_res %>% filter(Species %in% unique(spEvaluated$Species)) %>% 
    mutate(variable = "Observed")
  
  cat("Species that will be evaluated are:\n")
  cat(spEvaluated$Species, sep = "\n")
  
  
  
  # STEP 4 | Generate simulated datasets (null distributions) ####
  
  cat(paste(Sys.time(), "- Generating null distributions for:\n"))
  
  # Now we can generate a null distribution of the species repartition by randomly
  # resampling cells from the convex hull of the observed species distribution. 
  # The result is a collection of cells that represents what the species 
  # distribution could be if it was random and not linked to naturalness factors.
  
  
  # create an empty object to store all the null distributions
  nullFull <- NULL
  # sp=unique(spEvaluated$Species)[1]
  
  # loop over the species to evaluate
  for(sp in unique(spEvaluated$Species)) {
    
    # prompt species evaluated
    cat(paste(Sys.time(), sp))
    cat(paste0(" (", spEvaluated %>% filter(Species == sp) %>% 
                 select(spNum), "/", max(spEvaluated$spNum), ")\n"))
    
    
    # we draw a convex hull of the species distribution and extract the coordinates
    # of all the cells inside this hull. This list of cells will then be randomly
    # resampled to obtain null distributions of the species distribution.
    
    # create a convex hull of the species distribution
    sp_points <- sf::st_as_sf(x_evaluated %>% filter(Species == sp), 
                              coords = c("x", "y"), crs = 4326)
    
    
    
    # Comme les buffers en degrés sont imprécis, on transforme en projection métrique
    # EPSG 3035 : Lambert Europe
    sp_points_m <- st_transform(sp_points, 3035)
    
    # Créer des buffers de 20 km
    buffers <- st_buffer(sp_points_m, 20000)  # 10 km = 10000 m
    
    rm(sp_points_m)
    # Fusionner tous les buffers en un seul objet (MULTIPOLYGON)
    buffers_union <- st_union(buffers)
    
    # S'assurer que la géométrie est valide
    buffers_union <- st_make_valid(buffers_union)
    
    # Extraire uniquement les polygones
    buffers_poly <- st_collection_extract(buffers_union, "POLYGON")
    
    rm(buffers_union)
    # Reprojeter en WGS84 pour avoir lon/lat
    buffers_poly <- st_transform(buffers_poly, 4326)
    
    # Visualiser le résultat (facultatif)
    # plot(buffers_poly, border = "blue")
    # plot(sp_points, add = TRUE, col = "red", pch = 16)
    # plot(sp_resampled, add = TRUE, col = "green", pch = 16)
    
    
    # ggplot() +
    #   # polygone pays
    #   geom_sf(data = country_map_metro, fill = "lightblue", color = "black") +
    # 
    #   # polygones tampon
    #   geom_sf(data = buffers_poly,  fill = "darkgreen") +
    # 
    #   # points
    #   geom_sf(data = sp_points, color = "red", size = 2) +
    # 
    #   # thème
    #   theme_bw()
    
    
    
    
    # attribute raster cell numbers to cells of the convex hull
    # hullCoord <- tabularaster::cellnumbers(raster::raster(ras), buffers_poly) %>%
    #   rename(Cell = "cell_")
    Cell <- terra::cells(ras, vect(buffers_poly))[,2] #%>% rename(cell = "cell_")
    hullCoord <- as.data.frame(Cell) #######################################################################
    
    # convert cell numbers to XY coordinates
    hullCoord[,c("x","y")] <- terra::xyFromCell(ras, hullCoord$Cell)
    
    # select cells within the convex hull that are in the naturalness dataset 
    # (this excludes offshore points in cases where the convex hull includes marine
    # areas) 
    hullCoord %<>% filter(Cell %in% ras.value$Cell)
    
    # select kernel weights (sampling effort for each cell within the convex hull)
    kernel_weights <- data.frame(weight = kernel_sampling$w,
                                 Cell = seq(1:length(kernel_sampling$w))) %>%
      filter(Cell %in% hullCoord$Cell)
    
    
    
    
    # create an empty object to store the null distribution for this species
    nullSp <- NULL
    
    # loop over all the simulations requested
    for(i in 1:sim) {
      
      # prompt the simulation progress
      if (i == 1 & i == sim) {
        cat(paste0("Simulation ", as.numeric(i), "/", sim, "\n"))
      } else if(i == 1){
        cat(paste0("Simulation ", as.numeric(i), "/", sim, "... "))
      } else if (i == sim) {
        cat(paste0(" ", as.numeric(i), "/", sim, "\n"))
      } else if (i %% 100 == 0) {
        cat(paste0(" ", as.numeric(i), "/", sim, "... "))
      } 
      
      
      # resample cells within the convex hull, with their associated kernel weight
      # sp_resampled <- sample_n(kernel_weights, #x_visits_all,
      #                          size = spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent),
      #                          replace = FALSE,
      #                          weight = weight) # weight = nYearVisited
      
      
      ##Probleme kernel pour espèces artificielles####          
      
      if((sum(kernel_weights$weight > 0) < spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent))& (sp=="Anthropophobe species" | sp=="Synanthrope species")){         
        size_requested <- spEvaluated %>%
          filter(Species == sp) %>%
          pull(nCellsPresent)
        
        kernel_weights_pos <- kernel_weights %>%
          filter(weight > 0)
        
        size_final <- min(size_requested, nrow(kernel_weights_pos))
        
        sp_resampled <- sample_n(kernel_weights_pos,
                                 size = size_final,
                                 replace = FALSE,
                                 weight = weight)
        
      }else{ 
        sp_resampled <- sample_n(kernel_weights, #x_visits_all,
                                 size = spEvaluated %>% filter(Species == sp) %>% pull(nCellsPresent),
                                 replace = FALSE,
                                 weight = weight) # weight = nYearVisited
        
      } 
      
      ####Fin de la modif pour l'erreur 
      
      
      # assign simulation number and species
      sp_resampled %<>% mutate(simulation = i,
                               Species = sp,
                               variable = "Null")
      # add sampled data from this simulation to the dataset of this species
      nullSp <- rbind(nullSp, sp_resampled)
      
    } # end of simulations for this species
    
    # concatenate sampled data from all species
    nullFull <- rbind(nullSp, nullFull)
    
  } # end of the loop for all species 
  
  
  # add XY coordinates to null distributions
  nullFull[,c("x","y")] <- terra::xyFromCell(ras, nullFull$Cell)
  
  # bind the null and observed datasets
  datasetFinal <- rbind(nullFull %>% select(-weight),#-nYearVisited),  
                        x_evaluated %>% select(-SumAbundance) %>% mutate(simulation = "no"))
  datasetFinal %<>% mutate(Resolution = value) 
  
  datasetFinal %<>% mutate(Scale = names) ########################################################################
  
  # STEP 5 | Calculate effect size #### 
  
  # now that we have the observed and simulated distributions, we are going to
  # calculate the effect size of their difference in naturalness values
  
  cat(paste(Sys.time(), "- Calculating effect sizes...\n"))
  
  # create an empty object to store all effect sizes for all species
  effSizesFull <- NULL
  
  # for every simulation, calculate the effect size
  for(run in unique(nullFull$simulation)) {
    
    # prompt the calculation progress
    if (run == 1 & run == max(nullFull$simulation)) {
      cat(paste0("Simulation ", as.numeric(run), "/", max(nullFull$simulation), "\n"))
    } else if(run == 1){
      cat(paste0("Simulation ", as.numeric(run), "/", max(nullFull$simulation), "... "))
    } else if (run == sim) {
      cat(paste0(" ", as.numeric(run), "/", max(nullFull$simulation), "\n"))
    } else if (run %% 100 == 0) {
      cat(paste0(" ", as.numeric(run), "/", max(nullFull$simulation), "... "))
    } 
    
    
    # select a single simulation and the observed data, add naturalness raster values
    runN <- datasetFinal %>% dplyr::filter(simulation == run |
                                             simulation == "no") %>% 
      inner_join(ras.value, by = "Cell") 
    # TODO: note that some points don't have naturalness values here!
    # TODO: this ends up with different sample sizes for the simulation and the observed
    # data in some cases, visible in the third result data frame.
    
    # create an empty object to store the species effect sizes
    effSizesSp <- NULL
    
    for(sp in unique(runN$Species)) {
      
      # select one species
      runNSp <- runN %>% filter(Species == sp)
      # calculate the effect size: (simulated mean - observed mean)/estimated sd
      effSizesrunNSp <- runNSp %>% 
        rstatix::cohens_d(value ~ variable) %>% # equal = TRUE ?
        mutate(Species = sp,
               Run = run,
               Resolution = value,
               Scale=names)###########################################################################################
      # add effect size to the list
      effSizesSp <- rbind(effSizesSp, effSizesrunNSp)
      
    }
    
    # add species effect sizes to the list
    effSizesFull <- rbind(effSizesFull, effSizesSp)
    
  } #end of the loop by run
  
  
  # STEP 6 | Scoring ####
  
  # synthesize effect sizes per species
  effSizesFull_summary <- effSizesFull %>%
    group_by(Species) %>%
    summarise(mean = mean(effsize), 
              nRun = n(),
              Scale=Scale[1])
  
  # rescale effect sizes to obtain the synanthropy score per species
  effSizesFull_summary %<>% 
    mutate(Index = round(scales::rescale(mean, to = c(10, 1))),
           Resolution = value)
  
  
  # STEP 7 | Compile results for this resolution ####
  
  #1: species scores
  speciesScores <- rbind(speciesScores, effSizesFull_summary)
  
  #2: simulated and observed datasets
  samplesList <- rbind(samplesList, datasetFinal)
  
  #3: all effect sizes
  effSizes <- rbind(effSizes, effSizesFull)
  
  # create a list for all 3 tables
  results1 <- list("speciesScores" = speciesScores, "effSizes" = effSizes, "samplesList" = samplesList) 
  
  cat(paste(Sys.time(), "Analysis finished for resolution", value, "\n\n"))
  # end of resolution loop
  
  cat(paste(Sys.time(), "All done.\n"))    
  return(results1)
  })

# return(results1)

  
}



# 
# 
# sub_effsize_res <- effSizes %>% 
#   filter(Resolution == 10) %>% 
#   left_join(speciesScores, by = c("Species", "Resolution", "Scale"))
# 
# # make the score a factor
# sub_effsize_res$Index <- as.factor(sub_effsize_res$Index)
# 
# # plot the results by scale
# 
# for (country in unique (sub_effsize_res$Scale)){
#   
#   sub_effsize_res_country=subset(sub_effsize_res, Scale==country)
#   
#   t=ggplot(sub_effsize_res_country, 
#            aes(x = reorder(Species, -effsize), 
#                y = -effsize, fill = Index)) +
#     geom_hline(yintercept = 0.0, color = "darkgrey", 
#                linewidth = 0.8, linetype = "dashed") +
#     geom_boxplot() + 
#     coord_flip() +
#     scale_fill_brewer(name = "SSI", palette = "RdYlGn", direction=-1) +
#     ylab("Effect size") +
#     xlab("Species") + labs(title = paste0("Result for Squamata species in ", country, " scale"))+
#     theme(axis.title.y = element_blank()) +
#     theme(legend.position = c(0.9, 0.2)) +
#     theme_bw()
#   
#   print(t)
#   
# } 
# 

library(data.table)

results_list <- Filter(Negate(is.null), results1)


sp_ssi <- rbindlist(lapply(results_list, `[[`, "speciesScores"),
                    use.names = TRUE, fill = TRUE)

effsize_res <- rbindlist(lapply(results_list, `[[`, "effSizes"),
                    use.names = TRUE, fill = TRUE)

points <- rbindlist(lapply(results_list, `[[`, "samplesList"),
                    use.names = TRUE, fill = TRUE)


# (sp_ssi <- results[[2]]$speciesScores)
write.table(sp_ssi, "/scratch/bbongibault/Aves_continent_100sim.csv", row.names=FALSE, sep=";",dec=".", na=" ")


# head(effsize_res <- results[[2]]$effSizes)
write.table(effsize_res, "/scratch/bbongibault/Aves_effSizes_continent_100sim.csv", row.names=FALSE, sep=";",dec=".", na=" ")



# head(points <- results[[2]]$samplesList)
write.table(points, "/scratch/bbongibault/Aves_samplesList_continent_100sim.csv", row.names=FALSE, sep=";",dec=".", na=" ")





