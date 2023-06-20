library(dplyr)
library(reshape2)
library(stringr)

# *****************************
# Compile rarefaction data ----
# *****************************

outputFolder <- paste0(getwd(), "/sim/outputs")

rares <- list.files(path = outputFolder, pattern = ".rda", full.name = T)

rare.df <- data.frame()
for(r in 1:length(rares)){
  load(rares[r])
  
  rare_merged$population <- str_split_fixed(rares[r], "_", 4)[,2]
  
  #add nAnimals = 0 for each simulation
  df <- setNames(data.frame(matrix(ncol = ncol(rare_merged), 
                                   nrow = length(unique(rare_merged$simulation)))), 
                 c(names(rare_merged)))
  df$simulation <- 1:nrow(df)
  df[,c("nAnimals", "area_50p", "area_75p", "area_90p", "area_95p", "area_99p",
        "allDuration", "allMigs", "allAnimalYrs", "newAnimalDuration", "newAnimalMigs")] <- 0
  df$population <- rare_merged$population[1]
  df$projName <- rare_merged$projName[1]
  df$UD <- rare_merged$UD[1]
  
  rare_wZeros <- rbind(rare_merged, df)
  rare.df <- rbind(rare.df, rare_wZeros)
}

# ****************************
# Rate of increase, sq km ----
# ****************************

# area 50 is number of pixels
# convert to square kilometers
# sheep has 100m pixel resolution
# elk, mule deer, ph 500m pixels

rare.df$sqkm_50p <- (rare.df$area_50p*(500*500))/1000000
rare.df$sqkm_75p <- (rare.df$area_75p*(500*500))/1000000
rare.df$sqkm_90p <- (rare.df$area_90p*(500*500))/1000000
rare.df$sqkm_95p <- (rare.df$area_95p*(500*500))/1000000
rare.df$sqkm_99p <- (rare.df$area_99p*(500*500))/1000000

rare.df <- rare.df[with(rare.df, order(projName, UD, simulation, nAnimals)),]

# ****************
# Clean, save ----
# ****************

## make winter and summer range values 0.01 if they were zero

rare.df$sqkm_50p[rare.df$sqkm_50p == 0] <- .001
rare.df$sqkm_75p[rare.df$sqkm_75p == 0] <- .001
rare.df$sqkm_90p[rare.df$sqkm_90p == 0] <- .001
rare.df$sqkm_95p[rare.df$sqkm_95p == 0] <- .001
rare.df$sqkm_99p[rare.df$sqkm_99p == 0] <- .001

# *****************
# Long Version ----
# *****************

rare_long <- reshape2::melt(rare.df[,c("anID", "nAnimals", "simulation", "allDuration", "allMigs", "allAnimalYrs", "newAnimalDuration", "newAnimalMigs", "projName",
                              "UD", "population", 
                              "sqkm_50p", "sqkm_75p", "sqkm_90p", "sqkm_95p", "sqkm_99p"
                            )]
                   ,id.vars = c("anID", "nAnimals", "simulation", 
                               "allDuration", "allMigs", "allAnimalYrs", "newAnimalDuration", "newAnimalMigs", 
                               "projName", "UD", "population"
                               ))
names(rare_long)[names(rare_long) == "variable"] <- "Contour"
names(rare_long)[names(rare_long) == "value"] <- "area_sqkm"
rare_long$Contour <- gsub("sqkm_", "", rare_long$Contour) 

#take off "p" in contour and make numeric
unique(rare_long$Contour)
rare_long$Contour <- gsub("p", "", rare_long$Contour)
rare_long$Contour <- as.numeric(rare_long$Contour)

save(rare_long, file = paste0(getwd(), "/sim/RarefactionLONGWCovs.rda"))
