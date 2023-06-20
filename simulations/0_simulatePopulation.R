#https://www.danaseidel.com/MovEco-R-Workshop/Materials/Day7/Simulating_Movement/
require(plyr)
require(sp)
require(raster)

# Simulate Data ----
movement <- function(xy, step, heading) {
  
  #First we need to define pi
  pi = 3.141593
  
  #Then we split the starting point into an x_init and y_init
  x_init <- xy[1,1]
  y_init <- xy[1,2]
  
  #Here we translate the negative pi values into positive values
  #The headings now range from 0 to 2*pi
  if (heading < 0) {
    heading <- abs(heading) + pi
  }
  
  #Using this heading, and our sin function, we can solve for the change in y
  #Then we want to create a new value where we alter the starting y value
  y_change <- sin(heading)*step
  y_new <- y_init + y_change
  
  #The use cosine to determine the movement in x direction
  x_change <- cos(heading)*step
  x_new <- x_init + x_change
  
  #Finally, we create a data frame and save our new coordinates
  move.temp <- as.data.frame(matrix(0,1,4))
  move.temp[1,1] <- x_new
  move.temp[1,2] <- y_new
  move.temp[1,3] <- step
  move.temp[1,4] <- heading
  
  return(move.temp)
}

multi.move2 <- function(N, x, start.pts) {
  all.paths <- list()
  
  for (j in 1:N) {
    steps.df <- data.frame(matrix(0,100,4))
    steps.df[1,1:2] <- start.pts[j,]
    colnames(steps.df) <- c("x", "y", "step.length", "turn.angle")
    
    for (i in 2:x[j]) {
      step <- rnorm(n=1, mean=500, sd=100)
      heading <- runif(n=1, min=-pi, max=pi)
      next.pt <- movement(steps.df[(i-1),1:2], step, heading)
      steps.df[i,] <- next.pt
    }
    
    all.paths[[j]] <- steps.df
    all.paths[[j]]$date <- seq(from = as.POSIXct("2022-04-17 00:00:00", tz = "GMT"), by = 13*60*60, length.out = 100)
  }
  return(all.paths)
}

# Population size ----

## a) n = 250 ----
N = 250

## b) n = 500 ----
N = 500

## c) n = 1000? ----
N = 1000

start.pts <- data.frame(matrix(0,N,2))
colnames(start.pts) <- c("x", "y")

## gunnison utms are easting = 331060, northing = 4268386
gunniUTMX = 331060
gunniUTMY = 4268386

## 1. uniform starting distribution ----
simmy = 1

start.pts$x <- runif(n=N, min=gunniUTMX, max=gunniUTMX + 1e+05)
start.pts$y <- runif(n=N, min=gunniUTMY, max=gunniUTMY + 1e+05)
hist(start.pts$x)

multi.paths2 <- multi.move2(N, rep(100, N), start.pts)
names(multi.paths2) <- paste0(as.character(1:N))
df <- ldply(multi.paths2, rbind)

plot(df$x, df$y, col = df$.id)

## 2. bigger spread ----
simmy = 2

start.pts$x <- rnorm(n=N, gunniUTMX, sd = 5000)
start.pts$y <- rnorm(n=N, gunniUTMY, sd = 8000)
hist(start.pts$x)
multi.paths2 <- multi.move2(N, rep(100, N), start.pts)

names(multi.paths2) <- as.character(1:N)
df <- ldply(multi.paths2, rbind)
plot(df$x, df$y, col = df$.id)

## 3. concentrated ----
simmy = 3

start.pts$x <- rnorm(n=N, gunniUTMX, sd = 5)
start.pts$y <- rnorm(n=N, gunniUTMY, sd = 5)
hist(start.pts$x, main = paste("pop = ", N))
multi.paths2 <- multi.move2(N, rep(100, N), start.pts)

names(multi.paths2) <- as.character(1:N)
df <- ldply(multi.paths2, rbind)
plot(df$x, df$y, col = df$.id)

# *******************
# BBMM homerange ----
# *******************

## Inputs ----

# load whatever
cell.size = 500
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

### folders ----

ifelse(!dir.exists(paste0("D:/Thesis/Journal/MEE/sim/projFolders/pop",N)), 
       dir.create(paste0("D:/Thesis/Journal/MEE/sim/projFolders/pop",N)), "Folder exists already")

NewProjectFolder <- paste0("D:/Thesis/Journal/MEE/sim/projFolders/pop",N,"/sim",simmy)

ifelse(!dir.exists(NewProjectFolder), dir.create(NewProjectFolder), "Folder exists already")

ifelse(!dir.exists(paste0(NewProjectFolder,"/","sequencesWinter")), dir.create(paste0(NewProjectFolder,"/","sequencesWinter")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDsWinter")), dir.create(paste0(NewProjectFolder,"/","UDsWinter")), "Folder exists already")
ifelse(!dir.exists(paste0(NewProjectFolder,"/","UDs_popWinter")), dir.create(paste0(NewProjectFolder,"/","UDs_popWinter")), "Folder exists already")

seqs_fldr = paste0(NewProjectFolder, "/sequencesWinter")

BBs_fldr = paste0(NewProjectFolder, "/UDsWinter")
pop_BBs_out_fldr = paste0(NewProjectFolder, "/UDs_popWinter")
pop_Rasters_out_fldr = paste0(NewProjectFolder, "/popRasters")

## Make sequences ----

names(df)[1] <- "id"

#fix the dates  (it is OK to specify GMT, since all dates will be in GMT!)
df$date <- as.POSIXct(strptime(df$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
df$id <- as.character(df$id)
df$year <- as.numeric(strftime(df$date, format = "%Y", tz = "GMT"))
head(df)

#make data spatial
xy <- df[,c("x", "y")]
d <- SpatialPointsDataFrame(coords = xy, data = df, 
                            proj4string = CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 + towgs84=0,0,0"))
#transform to aea
d <- spTransform(d, CRS(out_proj))

#load study area grd that aligns with larger raster extent

ext <- raster::extent(d)
multiplyers <- c((ext[2]-ext[1])*mult4buff, (ext[4]-ext[3])*mult4buff)   # add about 20% around the edges of your extent (you can adjust this if necessary)
ext <- raster::extend(ext, multiplyers)
grd <- raster(ext)
res(grd) <- cell.size
projection(grd) <- proj4string(d)


#reduce the dataset to columns of interest
d <- d[,c("id", "date", "year")]

#build a database of the potential winter periods
ids <- unique(d$id)

for(i in 1:length(unique(ids))){
  tmp <- d[d$id == ids[i],]
  tmp <- cbind(as.data.frame(tmp)[,c("id","date")],coordinates(tmp))    # get it out of sp object
  tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
  names(tmp) <- c("id","date","x","y")   # rename columns
  if(nrow(tmp)>0){
    foreign::write.dbf(tmp, file = paste0(seqs_fldr,"/",paste0(ids[i],"_wi22"),".dbf"))    #write as dbf
  }
  else{
    next
  }
}

## Individual WR ----
#create database for files
fls <- list.files(seqs_fldr, ".dbf$")
print(paste0("You have ", length(fls), " sequences."))
d <- do.call(rbind, lapply(1:length(fls), function(i){
  db <- foreign::read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
  db$wint <- sub(".dbf","",fls[i])
  return(db)
}))

#loop through all individuals create asciis
u <- unique(d$wint)

require(parallel)
require(doSNOW)
regBB <- list()
Cores <- detectCores()
cl <- makeCluster(Cores - 2)
#cl <- makeCluster(2)
registerDoSNOW(cl)
start.time <- Sys.time()

indBB <- foreach(i = 1:length(u), .combine = rbind,
                 .packages = c("BBMM", "raster")) %dopar%{
                   #for(i in 1:length(u)){                   
                   temp <- d[d$wint == u[i],]
                   temp <- temp[order(temp$date),]
                   temp$date <- as.POSIXct(temp$date, "%Y-%m-%d %H:%M:%S", tz = "GMT")
                   
                   jul <- as.numeric(strftime(temp$date, format = "%j", tz = "GMT"))
                   
                       #prepare only the cells to run BB over
                       ext2 <- raster::extent(temp)
                       multiplyers <- c((ext2[2]-ext2[1])*mult4buff, (ext2[4]-ext2[3])*mult4buff)   # add about mult4buff around the edges of your extent (you can adjust this if necessary)
                       ext2 <- raster::extend(ext2, multiplyers)
                       cels <- raster::cellsFromExtent(grd, ext2)
                       
                       bb <- R.utils::withTimeout({
                           try(BrownianBridgeCustom(x=temp$x,
                                                    y=temp$y,
                                                    time.lag=diff(as.numeric(temp$date)/60),
                                                    area.grid=coordinates(grd)[cels,],
                                                    max.lag=max.lag*60,
                                                    time.step=time.step,
                                                    BMvar=BMvar,
                                                    location.error=location.error), #this is the location error of your collars
                               silent=TRUE)
                         }, envir=environment(), timeout = 10800, onTimeout = "warning")
                       
                       #set to 0 any values that are outside of the < 0.9999 contour
                       cutoff <- sort(bb$probability, decreasing=TRUE)
                       vlscsum <- cumsum(cutoff)
                       cutoff <- cutoff[vlscsum > .9999][1]
                       bb$probability[bb$probability < cutoff] <- 0
                       
                       #rescale probabilities so they equal 1
                       bb$probability <- bb$probability/sum(bb$probability)
                       
                       #output ASCII file
                       m <- data.frame(x=coordinates(grd)[,1], y=coordinates(grd)[,2],z=0)
                       m$z[cels] <- bb$probability
                       m <- SpatialPixelsDataFrame(points = m[c("x", "y")], data=m)
                       m <- as(m, "SpatialGridDataFrame")
                       write.asciigrid(m, paste(BBs_fldr,"/",u[i],"_ASCII.asc",sep=""), attr=3)
                       
                       
                       
                   }

stopCluster(cl)

## Population WR ----
fls <- list.files(BBs_fldr, ".asc$")

bb <- stack(paste(BBs_fldr, "/", fls, sep=""))
if(ncell(bb)*nlayers(bb) < 10000000){
    bb <- mean(bb)
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }else{
    beginCluster(n=cores)
    bb <- clusterR(bb, fun=calc, list(fun=sum))
    endCluster()
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }
  #output averaged individual ASCII file
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste(pop_BBs_out_fldr,"/","averageUD_winter.asc",sep=""), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste(pop_BBs_out_fldr,"/","averageUD_winter.img",sep=""), format="HFA") #change this fldr location to population loc for SR, WR, migs
  print(paste0("End time: ", Sys.time()))

plot(bb)
