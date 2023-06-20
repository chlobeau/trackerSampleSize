library(dplyr)
library(reshape2)

# *****************************
# Compile rarefaction data ----
# *****************************

rares <- list.files(path = "E:/Thesis/RarefactionOutputs", pattern = ".rda", full.name = T)

rare.df <- data.frame()
for(r in 1:length(rares)){
  load(rares[r])
  
  #add nAnimals = 0 for each simulation
  df <- setNames(data.frame(matrix(ncol = ncol(rare_merged), 
                                   nrow = length(unique(rare_merged$simulation)))), 
                 c(names(rare_merged)))
  df$simulation <- 1:nrow(df)
  df[,c("nAnimals", "area_50p", "area_75p", "area_90p", "area_95p", "area_99p",
        "allDuration", "allMigs", "allAnimalYrs", "newAnimalDuration", "newAnimalMigs")] <- 0
  df$projName <- rare_merged$projName[1]
  df$UD <- rare_merged$UD[1]
  
  rare_wZeros <- rbind(rare_merged, df)
  rare.df <- rbind(rare.df, rare_wZeros)
}

# *****************************
# Add mgmt unit covariates ----
# *****************************

## Population Estimates ----
popEsts <- read.csv("E:/Thesis/DataAnalysis/DAU_covariates/PopEstimates/ProjectsAnnualDemographics_wide_20220602.csv")

popEsts$popMean <- rowMeans(popEsts[,12:20], na.rm = T)

popEsts$UNIT %in% unique(rare.df$projName)

#MATCH
#df1$value <- df2$value[match(df1$name,df2$name)]

unique(rare.df$projName) %in% popEsts$UNIT
unique(rare.df$projName)[50]
# need to change RBS08 to RBS8
rare.df$projName[rare.df$projName == "RBS08"] <- "RBS8"

rare.df$spp <- popEsts$Species[match(rare.df$projName, popEsts$UNIT)]
rare.df$popEsts <- popEsts$popMean[match(rare.df$projName, popEsts$UNIT)]

## Winter range ----
wr <- read.csv("E:/Thesis/DataAnalysis/DAU_covariates/WinterRange/ProjectAnnualWinterRange_20220603.csv")
wr$ProjName <- gsub("-", "", wr$ProjName)
wr$ProjName <- gsub(" ", "", wr$ProjName)

unique(rare.df$projName) %in% wr$ProjName

rare.df$projArea_sqkm <- wr$area_sqkm[match(rare.df$projName, wr$ProjName)]
rare.df$projWR_prop <- wr$WR_prop[match(rare.df$projName, wr$ProjName)]

## Canopy cover ----
cc <- read.csv("E:/Thesis/DataAnalysis/DAU_covariates/CanopyCover/ProjectsAnnualCanopyCover_20220603.csv")
unique(rare.df$projName) %in% cc$ProjName

cc$ccMean <- rowMeans(cc[,2:8], na.rm = T)

rare.df$ccMean <- cc$ccMean[match(rare.df$projName, cc$ProjName)]

## Elevation ----
elev <- read.csv("E:/Thesis/DataAnalysis/DAU_covariates/Elev/ProjectsElevMetrics_20220603.csv")
unique(rare.df$projName) %in% elev$ProjName

rare.df$elev_range <- elev$elev_range[match(rare.df$projName, elev$ProjName)]
rare.df$elev_mean <- elev$elev_mean[match(rare.df$projName, elev$ProjName)]
rare.df$elev_95p <- elev$elev_95p[match(rare.df$projName, elev$ProjName)]

## Ecoregions ----
eco <- read.csv("E:/Thesis/DataAnalysis/DAU_covariates/Ecoregions/ProjectsEcoregions_20220603.csv")
unique(rare.df$projName) %in% eco$ProjName

eco$ProjName <- gsub("-", "", eco$ProjName)
eco$ProjName <- gsub(" ", "", eco$ProjName)

rare.df$L3_eco <- eco$L3_eco[match(rare.df$projName, eco$ProjName)]
rare.df$L2_eco <- eco$L2_eco[match(rare.df$projName, eco$ProjName)]

summary(rare.df)


# ADD SEX OF COLLARED INDS AND RATIOS for each simulation?

# ****************************
# Rate of increase, sq km ----
# ****************************

# area 50 is number of pixels
# convert to square kilometers
# sheep has 100m pixel resolution
# elk, mule deer, ph 500m pixels

rare.df$sqkm_50p <- NA
rare.df$sqkm_50p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"] <- (rare.df$area_50p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"]*(100*100))/1000000
rare.df$sqkm_75p <- NA
rare.df$sqkm_75p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"] <- (rare.df$area_75p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"]*(100*100))/1000000
rare.df$sqkm_90p <- NA
rare.df$sqkm_90p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"] <- (rare.df$area_90p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"]*(100*100))/1000000
rare.df$sqkm_95p <- NA
rare.df$sqkm_95p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"] <- (rare.df$area_95p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"]*(100*100))/1000000
rare.df$sqkm_99p <- NA
rare.df$sqkm_99p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"] <- (rare.df$area_99p[rare.df$spp == "RMBS" | rare.df$spp == "DBS"]*(100*100))/1000000

rare.df$sqkm_50p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"] <- (rare.df$area_50p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"]*(500*500))/1000000
rare.df$sqkm_75p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"] <- (rare.df$area_75p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"]*(500*500))/1000000
rare.df$sqkm_90p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"] <- (rare.df$area_90p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"]*(500*500))/1000000
rare.df$sqkm_95p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"] <- (rare.df$area_95p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"]*(500*500))/1000000
rare.df$sqkm_99p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"] <- (rare.df$area_99p[rare.df$spp == "MuleDeer" | rare.df$spp == "Elk" | rare.df$spp == "Pronghorn"]*(500*500))/1000000

rare.df <- rare.df[with(rare.df, order(projName, UD, simulation, nAnimals)),]

# ****************
# Clean, save ----
# ****************

# get rid of bighorn migration rarefaction where < 3 individuals migrate (S16, RBS24, S7)

migs <- rare.df[rare.df$UD == "mig",]
migs <- migs[migs$projName != "S16" & migs$projName != "S7" & migs$projName != "RBS24"
             & migs$projName != "S49" & migs$projName != "RBS24" & migs$projName != "S36",]
unique(migs$projName)

rare.df <- rare.df[rare.df$UD != "mig",]

rare.df <- rbind(rare.df, migs)

## make winter and summer range values 0.01 if they were zero

rare.df$sqkm_50p[rare.df$sqkm_50p == 0] <- .001
rare.df$sqkm_75p[rare.df$sqkm_75p == 0] <- .001
rare.df$sqkm_90p[rare.df$sqkm_90p == 0] <- .001
rare.df$sqkm_95p[rare.df$sqkm_95p == 0] <- .001
rare.df$sqkm_99p[rare.df$sqkm_99p == 0] <- .001

save(rare.df, file = "E:/Thesis/DataAnalysis/Rarefaction/RarefactionDataFrameWCovs_20220629.rda")

# *****************
# Long Version ----
# *****************

load(file = "E:/Thesis/DataAnalysis/Rarefaction/RarefactionDataFrameWCovs_20220629.rda")

rare_long <- melt(rare.df[,c("anID", "nAnimals", "simulation", "allDuration", "allMigs", "allAnimalYrs", "newAnimalDuration", "newAnimalMigs", "projName",
                              "UD", "spp", "popEsts", "projArea_sqkm", "projWR_prop", "ccMean", "elev_range", "elev_mean", "elev_95p", "L2_eco", "L3_eco",
                              "sqkm_50p", "sqkm_75p", "sqkm_90p", "sqkm_95p", "sqkm_99p")]
                   ,id.vars = c("anID", "nAnimals", "simulation", 
                               "allDuration", "allMigs", "allAnimalYrs", "newAnimalDuration", "newAnimalMigs", 
                               "projName", "UD", "spp", "popEsts", "projArea_sqkm", "projWR_prop", 
                               "ccMean", "elev_range", "elev_mean", "elev_95p",  "L2_eco", "L3_eco"))
names(rare_long)[names(rare_long) == "variable"] <- "Contour"
names(rare_long)[names(rare_long) == "value"] <- "area_sqkm"
rare_long$Contour <- gsub("sqkm_", "", rare_long$Contour) 

# make all bighorn one species
unique(rare_long$spp)

rare_long$spp2 <- rare_long$spp
rare_long$spp2[rare_long$spp == "DBS" | rare_long$spp == "RMBS"] <- "Bighorn"
unique(rare_long$spp2)

#take off "p" in contour and make numeric
unique(rare_long$Contour)
rare_long$Contour <- gsub("p", "", rare_long$Contour)
rare_long$Contour <- as.numeric(rare_long$Contour)

#add popDensity
rare_long$popDensity <- rare_long$popEsts/rare_long$projArea_sqkm

# make as.factor() covariates factors ahead of time
# factor(spp) factor(Contour) factor(L3_eco)
rare_long$spp_fact <- as.factor(rare_long$spp2)
rare_long$L3eco_fact <- as.factor(rare_long$L3_eco)
save(rare_long, file = "E:/Thesis/DataAnalysis/Rarefaction/RarefactionLONGWCovs_20220629.rda")
