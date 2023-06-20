library(arrangements) #permutations
#library(raster)
library(stringr)
library(doSNOW)
library(foreach)
library(parallel)
library(raster)
library(spatialEco)

projFolders <- "C:/Rarefaction/500m"
outputFolder <- "C:/Rarefaction/outputs/outputs_20220605"
projects <- list.files(projFolders)

for(r in 1:length(projects)){
  projName = projects[r]
  currentdate <- paste(format(Sys.time(), "%Y"),format(Sys.time(), "%m"),format(Sys.time(), "%d"), sep = "")
  
  #enter the input folder where seasonal UDs are stored
  seasonalBBs_fldr = paste(projFolders, projName, "UDsWinter", sep = "/")
  UDs_fldr = paste(projFolders, projName, "UDs_popWinter", sep = "/")
  
  fls <- list.files(UDs_fldr, ".asc$")
  fls_season <- list.files(seasonalBBs_fldr, ".asc$")
  #simulation_count <- 5
  simulation_count <- 100
  ids_unique <- unique(str_split_fixed(fls, "_", 3)[,1])
  
  ## generate random orderings of unique animals
  ## matrix rows = simulation number, columns = different animals (max number of individuals)
  sim_matrix <- data.frame(permutations(x = ids_unique, k = length(ids_unique), nsample = simulation_count))
  
  sim_matrix <- unique(sim_matrix) # get rid of duplicates in sim_matrix
  
  # *********************
  # Parallel process ----
  # *********************
  
  Cores <- detectCores()
  cl <- makeCluster(Cores - 4)
  #cl <- makeCluster(2) 
  registerDoSNOW(cl)
  
  Start.time <- Sys.time()
  
  ## enables to store multiple outputs
  rarefaction <- list()
  
  for(s in 1:nrow(sim_matrix)){
  #pp <- foreach(s = 1:nrow(sim_matrix), .combine = rbind,
  #              .packages = c("raster", "stringr", "spatialEco")) %dopar% {
  #                # create empty data frame for each simulation
                  simulation.df <- data.frame(nAnimals = 1:ncol(sim_matrix)
                                              , simulation = s
                                              , area_50p = NA
                                              , area_75p = NA
                                              , area_90p = NA
                                              , area_95p = NA
                                              , area_99p = NA
                                              , allDuration = NA
                                              , allMigs = NA
                                              , allAnimalYrs = NA
                                              , anID = NA
                  )
                  
                  ## create extent of population-level footprint for when popF > 1
                  
                  for(n in 1:nrow(simulation.df)) {
                    simulation.df$anID[n] <- sim_matrix[s,n] # new animal id
                    
                    if(n == 1){
                      raster1 <- raster(paste(UDs_fldr, fls[ids_unique == sim_matrix[s,1]], sep = "/"))
                      raster_stack <- stack(raster1)
                      
                    }
                    if(n > 1){
                      raster2 <- raster(paste(UDs_fldr, fls[ids_unique == sim_matrix[s,n]], sep = "/")) # stack rasters
                      raster_stack <- stack(raster_saved, raster2)
                      
                    }
                    
                    # sum and mean come up with the same vals, sum allows adding more individuals without making a big stack (faster)
                    raster_saved <- sum(raster_stack) #save the sum of everybody
                    raster_sum <- raster_saved/sum(values(raster_saved)) # make sure values add up to 1 before volume contour
                    plot(raster_sum, main = paste("n = ", n))
                    
                    #VOLUME CONTOUR at 99, 95, 90, 75, 50% isopleths
                    vol99.rasterized <- raster.vol(raster_sum, p = 0.99)
                    simulation.df$area_99p[n] <- cellStats(vol99.rasterized, sum)
                    #plot(vol99.rasterized, main = paste("99% vol contour | n = ", n))
                    
                    vol95.rasterized <- raster.vol(raster_sum, p = 0.95)
                    simulation.df$area_95p[n] <- cellStats(vol95.rasterized, sum)
                    #plot(vol95.rasterized, main = paste("95% vol contour | n = ", n))
                    
                    vol90.rasterized <- raster.vol(raster_sum, p = 0.90)
                    simulation.df$area_90p[n] <- cellStats(vol90.rasterized, sum)
                    #plot(vol90.rasterized, main = paste("90% vol contour | n = ", n))
                    
                    vol75.rasterized <- raster.vol(raster_sum, p = 0.75)
                    simulation.df$area_75p[n] <- cellStats(vol75.rasterized, sum)
                    #plot(vol75.rasterized, main = paste("75% vol contour | n = ", n))
                    
                    vol50.rasterized <- raster.vol(raster_sum, p = 0.50)
                    simulation.df$area_50p[n] <- cellStats(vol50.rasterized, sum)
                    #plot(vol50.rasterized, main = paste("50% vol contour | n = ", n))
                    
                    # add duration, number of an_yrs, and number of migrations for all individuals
                    fls_season_sub <- fls_season[str_split_fixed(fls_season, "_", 3)[,1] %in% sim_matrix[s,1:n]]
                    simulation.df$allDuration[n] <- length(unique(substr((str_split_fixed(fls_season_sub, "_",3)[,2]), 3,4)))
                    simulation.df$allMigs[n] <- NA #na for winter range
                    animalYrs <- paste(paste(str_split_fixed(fls_season_sub, "_",3)[,1],
                                             substr((str_split_fixed(fls_season_sub, "_",3)[,2]), 3,4),
                                             sep = "_")) # paste animalID_year
                    simulation.df$allAnimalYrs[n] <- length(unique(animalYrs))
                  }
                  rarefaction[[s]] <- simulation.df
                }
  
  stopCluster(cl)
  End.time <- Sys.time()
  print(paste(projName, "Start: ",Start.time,"| End: ",End.time, "| Minutes:", 
              round(difftime(End.time, Start.time, units = "mins"), 2)))
  
  individual.df <- data.frame(anID = ids_unique
                              , newAnimalDuration = NA
                              , newAnimalMigs = NA
                              )
  
  for(i in 1:nrow(individual.df)){
    #add number of years and number of migrations for individual
    fls_season_newAn_sub <- fls_season[str_split_fixed(fls_season, "_", 3)[,1] == individual.df$anID[i]]  
    individual.df$newAnimalDuration[i] <- length(unique(substr((str_split_fixed(fls_season_newAn_sub, "_",3)[,2]), 3,4)))
    individual.df$newAnimalMigs[i] <- NA #na for winter range
  }
  
  
  # **************************
  # Merge, add info, save ----
  # **************************
  
  ## Merge individual centroid x and y, individual number of years, and individual number of migs to larger dataframe
  rare_merged <- merge(pp, individual.df, all.x = T)
  
  ## add column for project name
  rare_merged$projName <- projName
  
  ## add column to designate Winter Range
  rare_merged$UD <- "winter"
  
  ## SAVE
  save(rare_merged, file = paste0(outputFolder, "/", "rarefWR_", projName, "_", currentdate,".rda"))
  write.csv(rare_merged, file = paste0(outputFolder, "/", "rarefWR_", projName, "_", currentdate,".csv"), row.names = F)
  save.image(file = paste0(outputFolder, "/", "rarefWR_", projName, "_", currentdate,".RData"))
}


