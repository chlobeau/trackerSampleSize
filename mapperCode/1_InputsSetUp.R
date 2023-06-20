# adapted from WMI Migration Mapper v2_3
# Wyoming Migration Initiative. 2017. Migration Mapper. University of Wyoming, Laramie, Wyoming.

library(raster)
library(stringr)
library(foreign)
library(sf)
library(rdist) #for calculating mig_dists
library(parallel)
library(doSNOW)
library(foreach)
library(rgdal)
library(sp)

#Step 1 Create nonWMI directory folders
#Step 2 Crop out study area
# 3a Corridors
# 3b Winter Range
# 3c Summer Range

# *************************
# INPUTS for all steps ----
# *************************

#create new project folder, store it in a place that makes sense (migMap fldr) and label with DAU and pixel size
NewProjectFolder <- "E:/migMap/250m/A8_250m"

# old input tab6 for one project
oldTabSixFldr <- "E:/migMap/TimingProjFldrs/A8/tabSixOutputs"

shpfl_fldr=paste0(NewProjectFolder,"/tabSixOutputs")
shpfl_name="pointsOut"
migtbl_name="migtime.csv"

# Raster cell size
cell.size=250 #50m, 100m, 250m, 500m

## BB Settings ----
# For when the project has collars with multiple fix intervals less than or greater than the FMV threshold
ci <- read.csv("E:/CollarMgmt/GPSCollarProcessing/Colorado_collarinfo_20220601.csv")
BMvar=1000 # 1400 elk, 1000 md
fixThreshold=5 # 3 elk, 5 md/bh

location.error=20
max.lag=27
time.step=5
mult4buff=0.2
contour=99
cores=2
mindays=30
out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

#Brownian bridge function... FOR FMV
BrownianBridgeCustom<-function (x, y, time.lag, location.error, area.grid = NULL, cell.size = NULL,
                                time.step = 10, max.lag = NULL, BMvar=BMvar)
{
  if (is.null(x) | is.null(y) | (length(x) != length(y))) {
    stop("data is missing or unequal number of x and y coordinates")
  }
  if (is.null(location.error))
    stop("must specify 'location.error'")
  if (is.null(area.grid) & is.null(cell.size)) {
    stop("'area.grid' or 'cell.size' must be specified")
  }
  if (!is.null(area.grid) & is.null(cell.size)) {
    cell.size <- abs(area.grid[1, 1] - area.grid[2, 1])
  }
  if (is.null(area.grid) & !is.null(cell.size)) {
    range.x <- range(x)
    range.y <- range(y)
    min.grid.x <- round(range.x[1] - 1 * sd(x))
    max.grid.x <- round(range.x[2] + 1 * sd(x))
    min.grid.y <- round(range.y[1] - 1 * sd(y))
    max.grid.y <- round(range.y[2] + 1 * sd(y))
    x. <- seq(min.grid.x, max.grid.x, cell.size)
    y. <- seq(min.grid.y, max.grid.y, cell.size)
    area.grid <- merge(x., y.)
  }
  if (is.null(max.lag)) {
    max.lag = max(time.lag) + 1
  }
  if (length(location.error) == 1) {
    location.error <- rep(location.error, length(x))
  }
  n.locs <- length(x)
  # BMvar <- brownian.motion.variance(n.locs, time.lag, location.error,
  #                                   x, y, max.lag)
  BMvar <- rep(BMvar, times = length(x))
  if (is.null(time.step))
    time.step <- 10
  grid.size <- nrow(area.grid)
  probability <- rep(0, grid.size)
  T.Total <- sum(time.lag)
  bbmm <- vector("list", 4)
  names(bbmm) <- c("Brownian motion variance", "x", "y", "probability")
  class(bbmm) <- "bbmm"
  probability <- NULL
  int <- 0
  for (i in 1:(n.locs - 1)) {
    if (time.lag[i] <= max.lag) {
      theta <- NULL
      tm <- 0
      while (tm <= time.lag[i]) {
        alpha <- tm/time.lag[i]
        mu.x <- x[i] + alpha * (x[i + 1] - x[i])
        mu.y <- y[i] + alpha * (y[i + 1] - y[i])
        sigma.2 <- time.lag[i] * alpha * (1 - alpha) *
          BMvar[i] + ((1 - alpha)^2) * (location.error[i]^2) +
          (alpha^2) * (location.error[i + 1]^2)
        ZTZ <- (area.grid[, 1] - mu.x)^2 + (area.grid[,
                                                      2] - mu.y)^2
        theta <- (1/(2 * pi * sigma.2)) * exp(-ZTZ/(2 *
                                                      sigma.2))
        int <- int + theta
        tm <- tm + time.step
      }
    }
  }
  probability <- int/T.Total
  probability <- probability/sum(probability)
  bbmm[[4]] <- probability
  bbmm[[1]] <- BMvar[1]
  bbmm[[2]] <- area.grid[, 1]
  bbmm[[3]] <- area.grid[, 2]
  return(bbmm)
}
#----------- END OF CUSTOM BB FUNCTION (for FMV)

# ************
# FOLDERS ----
# ************

ifelse(!dir.exists(paste0(NewProjectFolder,"/","metadata")), dir.create(paste0(NewProjectFolder,"/","metadata")), "Folder exists already")

ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDs")), dir.create(paste0(NewProjectFolder,"/","UDs")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDsSummer")), dir.create(paste0(NewProjectFolder,"/","UDsSummer")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDsWinter")), dir.create(paste0(NewProjectFolder,"/","UDsWinter")), "Folder exists already")

ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDs_pop")), dir.create(paste0(NewProjectFolder,"/","UDs_pop")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDs_popSummer")), dir.create(paste0(NewProjectFolder,"/","UDs_popSummer")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDs_popWinter")), dir.create(paste0(NewProjectFolder,"/","UDs_popWinter")), "Folder exists already")

ifelse(!dir.exists(paste0(NewProjectFolder,"/","sequences")), dir.create(paste0(NewProjectFolder,"/","sequences")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","sequencesSummer")), dir.create(paste0(NewProjectFolder,"/","sequencesSummer")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","sequencesWinter")), dir.create(paste0(NewProjectFolder,"/","sequencesWinter")), "Folder exists already")

ifelse(!dir.exists(paste0(NewProjectFolder,"/","footprints_pop")), dir.create(paste0(NewProjectFolder,"/","footprints_pop")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","popRasters")), dir.create(paste0(NewProjectFolder,"/","popRasters")), "Folder exists already")

ifelse(!dir.exists(paste0(NewProjectFolder,"/","tabSixOutputs")), dir.create(paste0(NewProjectFolder,"/","tabSixOutputs")), "Folder exists already")

# put projection info in metadata folder
write(out_proj, file = paste0(NewProjectFolder,"/","metadata/out_projection.txt"))

# ***********************
# COPY TAB 6 OUTPUTS ----
# ***********************

## ONE input project ----

# Read in migtime table from old project folder
mt <- read.csv(paste0(oldTabSixFldr, "/migtime.csv"))

#IF MIGRATION LASTED ENTIRE YEAR (bug in app)#
# change migration occured to zero #
# change start/end dates to NA  #
## spring mig
mt$springMig[as.Date(mt$endSpring) - as.Date(mt$startSpring) > 300] <- 0
mt[mt$springMig == 0 , c("startSpring", "endSpring")] <- NA

## fall mig
mt$fallMig[as.Date(mt$endFall) - as.Date(mt$startFall) > 300] <- 0
mt[mt$fallMig == 0 , c("startFall", "endFall")] <- NA

# save migtime table to new project
write.csv(mt, paste0(shpfl_fldr, "/migtime.csv"), row.names = F)

#copy all 4 shapefile files
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.dbf"),
          to = paste0(shpfl_fldr, "/pointsOut.dbf"))
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.prj"),
          to = paste0(shpfl_fldr, "/pointsOut.prj"))
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.shp"),
          to = paste0(shpfl_fldr, "/pointsOut.shp"))
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.shx"),
          to = paste0(shpfl_fldr, "/pointsOut.shx"))
#might have these is shp was edited?
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.sbn"),
          to = paste0(shpfl_fldr, "/pointsOut.sbn"))
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.cpg"),
          to = paste0(shpfl_fldr, "/pointsOut.shx"))
file.copy(from = paste0(oldTabSixFldr, "/pointsOut.sbn"),
          to = paste0(shpfl_fldr, "/pointsOut.sbn"))

## MULTIPLE input projects ----

#MERGE SHAPEFILES FROM MULTIPLE PROJECTS#

shpfl_fldr1="E:/migMap/TimingProjFldrs/D2/tabSixOutputs"
shpfl_fldr2="E:/migMap/TimingProjFldrs/D2_20220406/tabSixOutputs"
shpfl_name="pointsOut"
d1 <- st_read(shpfl_fldr1, shpfl_name)
d2 <- st_read(shpfl_fldr2, shpfl_name)

#bind rows (features) of sf projects
d <- rbind(d1, d2)
d <- as(d, "Spatial")
writeOGR(d, dsn = paste0(shpfl_fldr), layer = paste0(shpfl_name), driver = 'ESRI Shapefile', overwrite = TRUE)

#MERGE MIGRATION TIMING TABLES FROM MULTIPLE PROJECTS#
mt1 <- read.csv(paste(shpfl_fldr1, migtbl_name,sep="/"))
mt2 <- read.csv(paste(shpfl_fldr2, migtbl_name,sep="/"))
mt <- rbind(mt1, mt2)

#IF MIGRATION LASTED ENTIRE YEAR (bug in app)#
# change migration occured to zero #
# change start/end dates to NA  #
## spring mig
mt$springMig[as.Date(mt$endSpring) - as.Date(mt$startSpring) > 300] <- 0
mt[mt$springMig == 0 , c("startSpring", "endSpring")] <- NA

## fall mig
mt$fallMig[as.Date(mt$endFall) - as.Date(mt$startFall) > 300] <- 0
mt[mt$fallMig == 0 , c("startFall", "endFall")] <- NA

write.csv(mt, paste0(shpfl_fldr, "/migtime.csv"), row.names = F)
