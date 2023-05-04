ssi <- function(r, x, resolution, sim, threshold) {
  
  # Create files to compile the results ####
  
  all.score_by_species <- NULL #1
  all.sp_df_null_obs <-  NULL  #2
  all.results.EF.SSI <- NULL   #3
  
  for(value in resolution) {
    
    # STEP 1 | Aggregate the raster to the defined resolution
    
    cat(paste0(Sys.time(), "  -   Aggregate raster for resolution ", value, " "))
    
    ras <- raster :: aggregate(r, fact = value)
    
    cat(paste0(Sys.time(), "  -  Run analyses for the raster of ", value, " resolution "))
    
    #to extract cell numbers of XY coordinates
    sp_by_occ_raw <- x
    sp_by_occ_raw$Cell <- raster :: cellFromXY(ras, sp_by_occ_raw[,c("X","Y")])
    
    # STEP 2 | Estimate data set sampling effort ####
    
    # transform Abundance into Occurrence
    sp_by_occ_raw$Occurrence <- sp_by_occ_raw$Abundance  
    sp_by_occ_raw$Occurrence[sp_by_occ_raw$Occurrence > 1] <- 1
    
    # estimate the number of visits by years
    sp_by_occ <- sp_by_occ_raw  %>%
      group_by(Cell, Year) %>% #add Month if want assess sampling effort with more precision
      summarize(Occurrence  = sum(Occurrence)) %>%
      as.data.frame()
    
    sp_by_occ$Occurrence[sp_by_occ$Occurrence > 1] <- 1
    
    # estimate the number of visits during the periods
    sp_by_occ <- sp_by_occ  %>%
      group_by(Cell) %>%
      summarize(nVisit  = sum(Occurrence)) %>%
      as.data.frame()
    
    # compile a data.frame with all cell 
    all.cell <- data.frame(seq(1: raster :: ncell(ras)))
    colnames(all.cell) <- "Cell"
    sp_by_occ <- merge(sp_by_occ, all.cell, by = "Cell", all = TRUE )
    sp_by_occ[is.na(sp_by_occ)] <- 0
    sp_by_occ[,c("x","y")] <- raster :: xyFromCell(ras, sp_by_occ$Cell) 
    
    # estimate sampling effort by kernel density 
    kernel_sampling <- ks :: kde(x=sp_by_occ[,c("x","y")], w= sp_by_occ$nVisit)
    
    # STEP 3 | Compile data by occurrence ####
    # prepare the table of species occurrence
    sp_by_occ <- sp_by_occ_raw %>%
      group_by(Cell, Species) %>%
      summarize(Abundance  = sum(Abundance)) %>%
      as.data.frame()
    
    sp_by_occ[,c("x","y")] <- raster :: xyFromCell(ras,sp_by_occ$Cell)  
    
    # summary table for species information
    tab_sp <- count(sp_by_occ, Species) %>%
      group_by(Species) %>%
      as.data.frame()
    
    tab_sp$Species_evaluation  <- factor(ifelse(tab_sp$n > threshold, "evaluated", "no evaluated"))
    tab_sp <- subset(tab_sp, Species_evaluation == "evaluated")
    tab_sp$Species_num <- seq(1:nrow(tab_sp))
    tab_sp$Species <- as.factor(tab_sp$Species)
    
    species_to_keep <- tab_sp[tab_sp$n > threshold,"Species"]
    
    sp_by_occ <- sp_by_occ[(sp_by_occ$Species %in% species_to_keep),]
    
    # summary table for cell raster (obs)
    cell_obs <- data.frame(Cell = unique(sp_by_occ$Cell))
    
    ras.value <- data.frame(value = raster :: getValues(ras))
    ras.value$Cell <- rownames(ras.value)
    ras.value <- na.omit(ras.value)
    
    # STEP 4 | Assess null distributions ####
    
    all.null.distri <- NULL #for all species
    
    for(Species in levels(as.factor(tab_sp$Species))) {
      
      cat(paste0(Sys.time(), "  -  Species ", tab_sp[tab_sp$Species == Species,"Species_num"],
                 " (", round(100 * tab_sp[tab_sp$Species == Species,"Species_num"]/max(tab_sp$Species_num),2), "%)\n"))
      
      # create a convex hull of the species distribution
      sp_points <- sf :: st_as_sf(sp_by_occ[sp_by_occ$Species == Species,], coords = c("x", "y"), crs = 2154)
      sp_distri <- sf :: st_convex_hull(sf :: st_union(sp_points)) 
      sp_distri <- sf :: st_sf(geom = sp_distri)
      
      cell_ras_by_distri <- data.frame(tabularaster :: cellnumbers(ras, sp_distri))
      cell_ras_by_distri <- data.frame(Cell = cell_ras_by_distri[,2])
      
      cell_ras_by_distri[,c("x","y")] <- raster :: xyFromCell(ras, cell_ras_by_distri$Cell)
      cell_ras_by_distri <- merge(cell_ras_by_distri, ras.value, by = "Cell")
      cell_ras_by_distri <- na.omit(cell_ras_by_distri)
      cell_ras_by_distri <- data.frame(Cell = cell_ras_by_distri[,c("Cell")])
      
      prob <- data.frame(weight = kernel_sampling$w)
      prob$Cell <- rownames(prob)
      
      cell_ras_by_distri <- merge(cell_ras_by_distri, prob, by = "Cell")
      
      null.distri <- NULL
      
      for(i in 1:sim) {
        
        cat(paste0(Sys.time(), "  -   Run ", as.numeric(i),
                   " (", round(100 * as.numeric(i)/length(unique(null.distri$run)),2), "%)\n"))
        
        sp_resampled <- sample_n(cell_ras_by_distri, tab_sp[tab_sp$Species == Species,"n"], 
                                 replace = FALSE, weight = weight)
        sp_resampled$run <- i
        sp_resampled$Species <- Species # to assign the null distribution to the target species
        null.distri <- rbind(null.distri, sp_resampled)
        
      }
      
      all.null.distri <- rbind(null.distri, all.null.distri) #all species with method 1 & 2.1 & 2.2
      
    } #end of the loop for all species 
    
    all.null.distri$variable <- rep("Null",nrow(all.null.distri))
    sp_by_occ$variable <- rep("Obs", nrow(sp_by_occ))
    
    sp_df_null_obs <- all.null.distri
    sp_df_null_obs[,c("x","y")] <- raster :: xyFromCell(ras, sp_df_null_obs$Cell)
    
    sp_df_null_obs <- sp_df_null_obs[c("Cell", "Species", "x", "y","variable")]  
    sp_by_occ_tobind <- sp_by_occ
    
    sp_df_null_obs <- rbind(sp_df_null_obs, sp_by_occ_tobind[,!(names(sp_by_occ_tobind) %in% "Abundance")])
    
    # STEP 5 | Calculate effect size #### 
    
    all.EF <- NULL
    
    for(run in levels(as.factor(null.distri$run))) {
      
      sub.by.run <- all.null.distri[all.null.distri$run == run, ]
      col.to.delete.df1 <- c("run", "weight")
      col.to.delete.df2 <- c("Abundance","x","y")
      
      sub.by.run <- rbind(sub.by.run[,!(names(sub.by.run) %in% col.to.delete.df1)], sp_by_occ[,!(names(sp_by_occ) %in% col.to.delete.df2)])
      sub.by.run <- merge(sub.by.run, ras.value, by = "Cell")
      
      EF <- NULL
      
      for(Species in levels(as.factor(sub.by.run$Species))) {
        
        sub <- sub.by.run[sub.by.run$Species == Species, ]
        sub$variable <- as.factor(sub$variable)
        EF_sub <- sub %>% rstatix :: cohens_d(value ~ variable) # equal = TRUE ?
        EF_sub <- add_column(EF_sub, Species, .after = 0)
        EF <- rbind(EF, EF_sub)
        
      }
      
      EF$Run <- rep(run, nrow(EF))
      all.EF <- rbind(all.EF, EF)
      
    } #end of the loop by run
    
    synth_all.EF <- data.frame(all.EF)
    synth_all.EF$effsize <- as.numeric(synth_all.EF$effsize)
    
    score_by_species <- synth_all.EF %>%
      group_by(Species) %>%
      summarise(mean = mean(effsize), nRun = n())
    
    score_by_species <- data.frame(score_by_species)
    
    # STEP 6 | Scoring ####
    
    index <- round(scales :: rescale(score_by_species$mean, to = c(10, 1)))
    index <- data.frame(cbind(Species=score_by_species$Species, Index=as.numeric(index)))
    
    # STEP 7 | Compile results ####
    
    #1
    score_by_species <- merge(score_by_species, index, by = c("Species"))
    score_by_species$Resolution <- rep(value, nrow(score_by_species)) 
    all.score_by_species <- rbind(all.score_by_species, score_by_species)
    
    #2
    sp_df_null_obs$Resolution <- rep(value, nrow(sp_df_null_obs)) 
    all.sp_df_null_obs <- rbind(all.sp_df_null_obs, sp_df_null_obs)
    
    #3
    synth_all.EF$Resolution <- rep(value, nrow(synth_all.EF)) 
    all.results.EF.SSI <- rbind(all.results.EF.SSI, synth_all.EF)
    
    all.results.SSI <- list(all.score_by_species, all.results.EF.SSI, all.sp_df_null_obs) 
    
  } # end of resolution loop
  
  return(all.results.SSI)
  
} #end of the function