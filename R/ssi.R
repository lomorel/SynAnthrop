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
  
  
  # Sanity check in case of typos
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
  for(res in resolution) { 
    
    
    
    
    
    # STEP 1 | Naturalness raster file: Aggregate cells at desired resolution ####
    
    # raw raster files may be too large for the intended analysis. Cells are aggregated
    # to decrease the file complexity.
    cat(paste(Sys.time(), "- Aggregate naturalness raster cells at resolution", res, 
              ".\n"))
    # aggregate raster cells at defined resolution
    ras <- terra::aggregate(r, fact = res, fun = mean, na.rm = T)
    
    # Extract cell number and naturalness values for each cell
    ras.value <- data.frame(terra::values(ras)) %>% # extract values from raster
      rename(value = 1) %>% # rename column to "value"
      dplyr::mutate(cell = as.numeric(rownames(.))) %>% 
      na.omit() # remove cells without a naturalness value
    
    
    
    
    
    # STEP 2 | Dataset: calculate sampling effort per cell ####
    
    # we calculate the number of times each cell was visited per year and in total
    # in order to obtain an index of sampling effort. Sampling effort is linked with
    # species detection, i.e., cells that were heavily prospected have very reliable
    # species lists. 
    
    # aggregate observations in cells at the defined raster resolution, based on XY coordinates
    data_res <- data # create a dataset for this resolution
    data_res %<>% dplyr::mutate(cell = terra::cellFromXY(ras, data_res %>% dplyr::select(X, Y))) # add cell# based on XY coordinates of the observation
    missingcoord <- dim(data_res %>% filter(!cell %in% ras.value$cell))[1] # check how many observations are not on the spatial extent of the map
    data_res %<>% dplyr::filter(cell %in% ras.value$cell) # select only observations for which there is a naturalness score
    cat(paste(Sys.time(), missingcoord, "observation(s) were eliminated because of missing naturalness values\n"))


    # If Year is available, calculate number of years of visits
    if("Year" %in% colnames(data_res)) { 
      # calculate the number of years each cell has been visited
      x_visits <- data_res %>%
        dplyr::group_by(cell) %>%
        dplyr::summarize(nYearVisited = n_distinct(Year))
    
    } else {
      # attribute one year per cell
      x_visits <- data_res %>%
        distinct(Cell) %>%
        mutate(nYearVisited = 1)
    }
      

      # compile a data.frame with the number of visits for all raster cells, including non-visited cells (0)
      x_visits_all <- ras.value %>% # start from cells
        dplyr::full_join(x_visits, by = "cell") %>%  # join cells and visits
        dplyr::mutate(nYearVisited = ifelse(is.na(nYearVisited), 0, nYearVisited)) %>%  # replace NA (unvisited cell) by 0 
        dplyr::mutate(data.frame(terra::xyFromCell(ras, .$cell)))
      
      cat(paste(Sys.time(), "At this resolution, the landscape is divided into", length(unique(x_visits_all$x)), 
                "*", length(unique(x_visits_all$y)), 
                "cells (total =", dim(x_visits_all)[1], "cells).\n You chose to evaluate species present in more than", 
                threshold, "cells.\n"))
      
      
      # interpolate sampling effort by kernel density 
      suppressWarnings(kernel_sampling <- ks::kde(x = x_visits_all %>% select(x, y),
                                                  w = x_visits_all$nYearVisited))
      # the "Weights don't sum to sample size - they have been scaled accordingly" warning is silenced
      
      
      
    
    
    # STEP 3 | Dataset: identify species to evaluate ####
    
    # Species seldom detected will produce unreliable synanthropy scores. 
    # We will only evaluate the synanthropy score of species that have been detected 
    # in at least 'threshold' cells.
    
    # sum the abundances per species per cell
    data_res %<>%
      dplyr::group_by(cell, Species) %>%
      dplyr::summarize(SumAbundance = sum(Abundance)) %>% # abundance per species per cell
      dplyr::left_join(data.frame(cell = unique(data_res$cell), terra::xyFromCell(ras, unique(data_res$cell))), by = "cell") # add XY coordinates
    
    # count the number of cells in which each species was detected
    spDetection <- data_res %>% 
      dplyr::group_by(Species) %>%
      dplyr::summarize(nCellsPresent = dplyr::n_distinct(cell))
    
    # species will only be evaluated if they have been detected in more than 
    # 'threshold' cells
    spDetection %<>% dplyr::mutate(evaluation = ifelse(nCellsPresent > threshold, 
                                                "evaluated", "not evaluated"))
    spEvaluated <- spDetection %>% dplyr::filter(evaluation == "evaluated")
    # stop the execution if no species meets this criterion
    if (dim(spEvaluated)[1] == 0) {
      stop(paste(Sys.time(), "No species is present in the minimal requested number of cells. Consider lowering threshold and/or resolution. Execution stopped."))
    }
    
    # number the species to be evaluated
    spEvaluated %<>% dplyr::mutate(spNum = seq(1:nrow(spEvaluated)))
    
    # subset the dataset to keep only species that will be evaluated
    x_evaluated <- data_res %>% dplyr::filter(Species %in% unique(spEvaluated$Species)) %>% 
      dplyr::mutate(variable = "Observed")
    
    cat("Species that will be evaluated are:\n")
    cat(spEvaluated$Species, sep = "\n")
    
    
    
    
    
    # STEP 4 | Generate simulated datasets (null distributions) ####
    
    cat(paste("\n", Sys.time(), "- Generating null distributions for:\n"))
    
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
      cat(paste0(" (", spEvaluated %>% dplyr::filter(Species == sp) %>% 
                   select(spNum), "/", max(spEvaluated$spNum), ")\n"))
      
      
      # we draw a convex hull of the species distribution and extract the coordinates
      # of all the cells inside this hull. This list of cells will then be randomly
      # resampled to obtain null distributions of the species distribution.
      
      # create a convex hull of the species distribution
      sp_points <- sf::st_as_sf(x_evaluated %>% filter(Species == sp), 
                                coords = c("x", "y"), crs = sf::st_crs(ras))
      sp_hull <- sf::st_convex_hull(sf::st_union(sp_points)) 
      sp_hull_terra <- vect(sp_hull) # convert the convex hull to an object compatible with terra
      
      # attribute raster cell numbers to cells of the convex hull
      hullCoord <- terra::extract(ras, sp_hull_terra, cells = T) %>% 
        select(-ID) %>% dplyr::rename(value = 1) %>% 
        dplyr::left_join(data.frame(cell = .$cell, terra::xyFromCell(ras, unique(.$cell))), by = "cell") # convert cell numbers to XY coordinates

      # select cells within the convex hull that are in the naturalness dataset 
      # (this excludes offshore points in cases where the convex hull includes marine areas) 
      hullCoord %<>% dplyr::filter(cell %in% ras.value$cell)

      
      
      # select kernel weights (sampling effort for each cell within the convex hull)
      kernel_weights <- data.frame(weight = kernel_sampling$w,
                         cell = x_visits_all$cell) %>%
        dplyr::filter(cell %in% hullCoord$cell)
      
      
    
      
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
        sp_resampled <- dplyr::sample_n(kernel_weights, 
                                 size = spEvaluated %>% dplyr::filter(Species == sp) %>%
                                   dplyr::pull(nCellsPresent), 
                                 # here we use the number of cells retained by the naturalness map 
                                 # to have the same number of obs and simulated data
                                 replace = FALSE,
                                 weight = weight) # weight = nYearVisited
        
        # assign simulation number and species
        sp_resampled %<>% dplyr::mutate(simulation = i,
                                Species = sp,
                                variable = "Null")
        # add sampled data from this simulation to the dataset of this species
        nullSp <- rbind(nullSp, sp_resampled)
        
      } # end of simulations for this species
      
      # concatenate sampled data from all species
      nullFull <- rbind(nullSp, nullFull)
      
    } # end of the loop for all species 
    
    
    # add XY coordinates to null distributions
    nullFull %<>% left_join(hullCoord %>% select(-value), by = "cell")
    
    # bind the null and observed datasets
    datasetFinal <- rbind(nullFull %>% dplyr::select(-weight), 
                          x_evaluated %>% 
                            dplyr::select(-SumAbundance) %>% 
                            dplyr::mutate(simulation = "no"))
    datasetFinal %<>% dplyr::mutate(Resolution = res) 
    
    
    
    
    
    # STEP 5 | Calculate effect size #### 
    
    # now that we have the observed and simulated distributions, we are going to
    # calculate the effect size of their difference in naturalness values
    
    cat(paste("\n", Sys.time(), "- Calculating effect sizes...\n"))
    
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
        dplyr::inner_join(ras.value, by = "cell") 

      # create an empty object to store the species effect sizes
      effSizesSp <- NULL
      
      # loop over every species 
      for(sp in unique(runN$Species)) {
        
        # select one species
        runNSp <- runN %>% dplyr::filter(Species == sp)
        # calculate the effect size: (simulated mean - observed mean)/estimated sd
        effSizesrunNSp <- runNSp %>% 
          rstatix::cohens_d(value ~ variable) %>% 
          dplyr::mutate(Species = sp,
                 Run = run,
                 Resolution = res)
        # add effect size to the list
        effSizesSp <- rbind(effSizesSp, effSizesrunNSp)
        
      }
      
      # add species effect sizes to the list
      effSizesFull <- rbind(effSizesFull, effSizesSp)
      
    } #end of the loop by run
    
    
    
    
    
    # STEP 6 | Scoring ####

    # synthesize effect sizes per species
    effSizesFull_summary <- effSizesFull %>%
      dplyr::group_by(Species) %>%
      dplyr::summarise(mean = mean(effsize), 
                nRun = n())
    
    # rescale effect sizes between 1 and 10 to obtain the synanthropy score per species
    effSizesFull_summary %<>% 
      dplyr::mutate(Index = round(scales::rescale(mean, to = c(10, 1))),
             Resolution = res)
    
    
    
    
    
    # STEP 7 | Compile results for this resolution ####
    
    #1: species scores
    speciesScores <- rbind(speciesScores, effSizesFull_summary)
    
    #2: simulated and observed datasets
    samplesList <- rbind(samplesList, datasetFinal)
    
    #3: all effect sizes
    effSizes <- rbind(effSizes, effSizesFull)
    
    cat(paste(Sys.time(), "Analysis finished for resolution", res, "\n"))
  } # end of resolution loop
  

  # create a list for all 3 tables
  results <- list("speciesScores" = speciesScores, "effSizes" = effSizes, "samplesList" = samplesList) 
  cat(paste(Sys.time(), "Analysis completed."))
  
  return(results)
  
} # end of the function
