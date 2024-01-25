#'Generate null distributions and calculate synanthropy index
#'
#'@export
#'
#'@param r The raster to be analysed (default name: raster)
#'@param data The species dataset to be processed, containing columns named
#'species, year, X, and Y (default name: dataset)
#'@param resolution Resolution(s) to be tested, specified as the aggregation factor, 
#'i.e., the number of raster cells in each direction (horizontally and vertically) (default value: 200)
#'@param sim Number of simulations (null distributions) to be generated (default value: 500)
#'@param threshold Threshold of species occurrence under which the species will 
#'not be accounted for (default value: 30)
#'
#'@return Three data frames: 
#'`speciesScores` a short summary table (one line per species) with the mean of all effect 
#'size and the corresponding index (range between 1 to 10) for each resolution. 
#'The number of runs used to calculate the mean effect size is also specified. 
#'`effSizes` a data.frame compiling all the raw results, i.e, all the effect sizes 
#'per run, with corresponding information provided by the function cohens_d 
#'(rstatix package); n1 and n2 correspond to the number of occurrences compared 
#'(n1 for the null distribution and n2 for observed data).
#'`samplesList` a data.frame archiving all the occurrence randomly drawn to assess scores.

#'@examples
#' example <- ssi(dataset)



ssi <- function(r = raster, 
                data = dataset, 
                resolution = 200, 
                sim = 500, 
                threshold = 30) {
  
  # Create objects to store the results
  speciesScores <- NULL # store species scores
  effSizes <- NULL # store all effect sizes calculated
  samplesList <-  NULL # store all samples drawn and observed
  
  
  # Safety check in case of typos
  cat("Species found in the dataset:\n")
  cat(sort(unique(data$Species)), sep = "\n")
  if (interactive() == TRUE) {
    var = readline("\nDo you wish to continue with this species list? (Y/N)")
    if (var == "N" | var == "n") {
      stop("Execution stopped.")
    }
  }
  
  # Check if a temporal scale is available for weighing sampling effort
  if(!"Year" %in% colnames(data)) { 
    cat("No \"Year\" variable available, all observations will be considered
        from the same year \n")
  } 
  
  
  # loop over each value listed in argument "resolution"
  for(value in resolution) { 
    
    
    # STEP 1 | Naturalness raster file: Aggregate cells at desired resolution ####
    
    # raw raster files may be too large for the intended analysis. Cells are aggregated
    # to decrease the file complexity.
    cat(paste(Sys.time(), "- Aggregate naturalness raster cells at resolution", value, 
              "(this step may take a few minutes)\n"))
    # aggregate raster cells at defined resolution
    ras <- terra::aggregate(r, fact = value)
    
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
    missingcoord <- dim(data_res %>% filter(is.na(Cell)))[1]
    data_res %<>% filter(!is.na(Cell))
    cat(paste(missingcoord, "observation(s) were eliminated because of mismatch with the raster file\n"))

    # If Year is available, calculate number of years of visits
    if("Year" %in% colnames(data_res)) { 
      # calculate the number of years each cell has been visited
      x_visits <- data_res %>%
        group_by(Cell) %>%
        summarize(nYearVisited = n_distinct(Year))
    
    } else {
      # attribute one year per cell
      x_visits <- data_res %>%
        distinct(Cell) %>%
        mutate(nYearVisited = 1)
    }
      
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
                                coords = c("x", "y"), crs = sf::st_crs(ras))
      sp_hull <- sf::st_convex_hull(sf::st_union(sp_points)) 

      # attribute raster cell numbers to cells of the convex hull
      hullCoord <- tabularaster::cellnumbers(ras, sp_hull) %>%
        rename(Cell = "cell_")
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
        
        
        # resample cells within the convex hull with their associated kernel weight
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
                 Resolution = value)
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
                nRun = n())
    
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
    
    cat(paste(Sys.time(), "Analysis finished for resolution", value, "\n"))
  } # end of resolution loop
  
  cat(paste(Sys.time(), "Analysis completed."))
  
  return(results)
  
} # end of the function
