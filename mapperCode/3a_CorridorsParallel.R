# adapted from WMI Migration Mapper v2_3
# Wyoming Migration Initiative. 2017. Migration Mapper. University of Wyoming, Laramie, Wyoming.

# ***********
# Inputs ----
# ***********

# NewProjectFolder from step 1

seqs_fldr= paste0(NewProjectFolder, "/sequences")
metadata_fldr=paste0(NewProjectFolder, "/metadata")

lns_out_fldr = paste0(NewProjectFolder, "/migrationLines")

BBs_fldr = paste0(NewProjectFolder, "/UDs")
pop_BBs_out_fldr= paste0(NewProjectFolder, "/UDs_pop")
pop_Rasters_out_fldr = paste0(NewProjectFolder, "/popRasters")
pop_footprint_out_fldr = paste0(NewProjectFolder, "/footprints_pop")

#load study area grd that aligns with larger raster extent
load(paste0(NewProjectFolder, "/studyarea_grd.rda"))
cell.size #from step 1 make sure this matches studyarea_grd resolution
res(grd.studyarea)


# ***************************
# 1. Migration Sequences ----
# ***************************
## based on C:\MAPP_2.3\app\scriptsv2\create.seqs.R
# This function takes the outputs from Migration Mapper tab 6 and creates
# dbf files of each individual animal's migration sequences

# read in migration timing

mt <- read.csv(paste(shpfl_fldr, migtbl_name,sep="/"))

# read in the shapefile

d <- st_read(shpfl_fldr, shpfl_name)

d <- as(d, "Spatial")

print(paste0("Your shapefiles has ", nrow(d), " rows."))

#reproject to new projection, if necessary
proj <- proj4string(d)
if(proj != out_proj){
  d <- spTransform(d, CRS(out_proj))
}

#reduce the dataset to columns of interest
d <- d[,c("newUid","nwMstrD")]
names(d) <- c("id","date")

#fix the dates  (it is OK to specify GMT, since all dates will be in GMT!)
d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
mt$startSpring <- as.POSIXct(strptime(paste(mt$startSpring, "00:00:01"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
mt$endSpring <- as.POSIXct(strptime(paste(mt$endSpring, "23:59:59"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
mt$startFall <- as.POSIXct(strptime(paste(mt$startFall, "00:00:01"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
mt$endFall <- as.POSIXct(strptime(paste(mt$endFall, "23:59:59"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")

# loop through spring migrations
for(i in 1:nrow(mt)){
  if(mt$springMig[i]==1 & is.na(mt$startSpring[i])==FALSE){ ##jgage
    tmp <- d[d$id == mt$newUid[i] & d$date >= mt$startSpring[i] & d$date <= mt$endSpring[i],]
    tmp <- as.data.frame(tmp)    # get it out of sp object
    tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
    names(tmp) <- c("id","date","x","y")   # rename columns
    if(nrow(tmp)>0){ ##jgage
      write.dbf(tmp, file = paste0(seqs_fldr,"/",mt$newUid[i], "_sp",substr(mt$nsdYear[i],3,4),".dbf"))    #write as dbf        ##jgage
    } ##jgage
  }else{
    next
  }
}

# loop through fall migrations
for(i in 1:nrow(mt)){
  if(mt$fallMig[i]==1 & is.na(mt$startFall[i])==FALSE){ ##jgage
    tmp <- d[d$id == mt$newUid[i] & d$date >= mt$startFall[i] & d$date <= mt$endFall[i],]
    tmp <- as.data.frame(tmp)    # get it out of sp object
    tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
    names(tmp) <- c("id","date","x","y")   # rename columns
    if(nrow(tmp)>0){ ##jgage
      write.dbf(tmp, file = paste0(seqs_fldr,"/",mt$newUid[i], "_fa",substr(mt$nsdYear[i],3,4),".dbf"))    #write as dbf    #write as dbf        ##jgage
    } ##jgage
  }else{
    next
  }
}

# ******************************
# 2. Individual Migrations ----
# ******************************

## based on C:\MAPP_2.3\app\scriptsv2\create.BBs.R
# This function takes a folder of migration sequences (in .dbf format)
# and creates regular brownian bridges from them
# location.error is in meters, cell size is in meters. max.lag is in hours
# contour is the % contour around UD for the footprint. Usually 99

#load up the data into a single database
fls <- list.files(seqs_fldr, ".dbf$")
print(paste0("You have ", length(fls), " sequences."))
d <- do.call(rbind, lapply(1:length(fls), function(i){
  db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
  db$mig <- sub(".dbf","",fls[i])
  return(db)
}))
#check and make sure the columns are correct.
all(c("date","x","y") %in% names(d))

d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")

#if there are less than 4 rows of d, stop - not enough seqs to calculate BBs
nrow(d) < 4

u <- unique(d$mig)

# need to create a vector if id_yrs for linking to migtimes table ####jgage
# this will then be a new column in the metadata table
idsSeasonsYears<-strsplit(u,'_') ### jgage
uids<-lapply(idsSeasonsYears, `[[`, 1) ### jgage
uids<-unlist(uids) ### jgage
years<-sapply(idsSeasonsYears, tail, 1) ### jgage
substrRight <- function(x, n){ ### jgage
  substr(x, nchar(x)-n+1, nchar(x)) ### jgage
} ### jgage
years<-paste0('20',substrRight(years, 2)) ### jgage
id_yrs<-paste0(uids,'_',years) ### jgage

#designate BMvar for individuals based on collar fix interval
uBMvar <- uids
for(i in 1:length(uBMvar)){
  #fixInt <- ci$FixInterval_hrs[uids[i]]
  fixInt <- ci$FixInterval_hrs[ci$Animal_ID == uids[i]]
  
  if(fixInt > fixThreshold | is.na(fixInt)){
    uBMvar[i] <- BMvar
  }
  else{
    uBMvar[i] <- 0
  }
}
uBMvar #look at bmvars

###USES GRD.STUDYAREA INSTEAD OF GRD###

regBB <- list()
cl <- makeCluster(cores)
registerDoSNOW(cl)
start.time <- Sys.time()

indBB <- foreach(i = 1:length(u), .combine = rbind,
                 .packages = c("BBMM", "raster")) %dopar%{
#for(i in 1:length(u)){
  temp <- d[d$mig==u[i],]
  temp <- temp[order(temp$date),]
  
  jul <- as.numeric(strftime(temp$date, format = "%j", tz = "GMT"))
  if(nrow(temp) <= length(unique(jul))){
    BB_meta <- data.frame(input.file=u[i],
                        brownian.motion.variance=NA,
                        grid.size=NA,
                        grid.cell.size=NA,
                        date.created=Sys.time(),
                        execution_time=paste0(round(difftime(Sys.time(), start.time, units="min"),2)," minutes"),
                        num.locs=nrow(temp),
                        Start.Date=min(temp$date),
                        End.Date=max(temp$date),
                        num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                        errors="Only 1 point per day, on average.",
                        notes=mt[which(mt$id_yr==id_yrs[i]),'notes'])
  }  #filter out animals with only 1 pt per day on average
  else{
    if(nrow(temp) < 4){
      BB_meta <- data.frame(input.file=u[i],
                          brownian.motion.variance=NA,
                          grid.size=NA,
                          grid.cell.size=NA,
                          date.created=Sys.time(),
                          execution_time=paste0(round(difftime(Sys.time(), start.time, units="min"),2)," minutes"),
                          num.locs=nrow(temp),
                          Start.Date=min(temp$date),
                          End.Date=max(temp$date),
                          num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                          errors="Less than 4 points.",
                          notes=mt[which(mt$id_yr==id_yrs[i]),'notes'])
    } #also filter out animals w <4 pts
    else{
      #prepare only the cells to run BB over
      ext2 <- raster::extent(temp)
      multiplyers <- c((ext2[2]-ext2[1])*mult4buff, (ext2[4]-ext2[3])*mult4buff)   # add about mult4buff around the edges of your extent (you can adjust this if necessary)
      ext2 <- raster::extend(ext2, multiplyers)
      cels <- cellsFromExtent(grd.studyarea, ext2)
      
      # this is the function to calculate the regular BB
      
      if(uBMvar[i] == 0){
        bb <- R.utils::withTimeout({
          try(BBMM::brownian.bridge(x=temp$x,
                                    y=temp$y,
                                    time.lag=diff(as.numeric(temp$date)/60),
                                    area.grid=coordinates(grd.studyarea)[cels,],
                                    max.lag=max.lag*60,
                                    time.step=time.step,
                                    location.error=location.error), #this is the location error of your collars
              silent=TRUE)
        }, envir=environment(), timeout = 10800, onTimeout = "warning")
      }else{ #THIS IS FMV code
        bb <- R.utils::withTimeout({
          try(BrownianBridgeCustom(x=temp$x,
                                   y=temp$y,
                                   time.lag=diff(as.numeric(temp$date)/60),
                                   area.grid=coordinates(grd.studyarea)[cels,],
                                   max.lag=max.lag*60,
                                   time.step=time.step,
                                   BMvar=BMvar,
                                   location.error=location.error), #this is the location error of your collars
              silent=TRUE)
        }, envir=environment(), timeout = 10800, onTimeout = "warning")
      }
      #set to 0 any values that are outside of the < 0.9999 contour
      cutoff <- sort(bb$probability, decreasing=TRUE)
      vlscsum <- cumsum(cutoff)
      cutoff <- cutoff[vlscsum > .9999][1]
      bb$probability[bb$probability < cutoff] <- 0
      
      #rescale probabilities so they equal 1
      bb$probability <- bb$probability/sum(bb$probability)
      
      #output ASCII file
      m <- data.frame(x=coordinates(grd.studyarea)[,1], y=coordinates(grd.studyarea)[,2],z=0)
      m$z[cels] <- bb$probability
      m <- SpatialPixelsDataFrame(points = m[c("x", "y")], data=m)
      m <- as(m, "SpatialGridDataFrame")
      write.asciigrid(m, paste0(BBs_fldr,"/",u[i],"_ASCII.asc"), attr=3)
      
      BB_meta <- data.frame(input.file=u[i],
                            brownian.motion.variance=round(bb[[1]],2),
                            grid.size=length(bb$x),
                            grid.cell.size=abs(bb$x[1]-bb$x[2]),
                            date.created=Sys.time(),
                            execution_time=paste0(round(difftime(Sys.time(), start.time, units="min"),2)," minutes"),
                            num.locs=nrow(temp),
                            Start.Date=min(temp$date),
                            End.Date=max(temp$date),
                            num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                            errors="None",
                            notes=mt[which(mt$id_yr==id_yrs[i]),'notes'])
      
    }
  }
  regBB[[i]] <- BB_meta
}

stopCluster(cl)
write.csv(indBB, file=paste0(metadata_fldr,"/metadata_migration.csv"), row.names=FALSE)

# ****************************************
# 3. Population and ind. average migs ----
# ****************************************
## based on C:\MAPP_2.3\app\scriptsv2\create.BBs.avgs.R
# This function takes in individual BBs and creates averaged individual BBs
# and population UDs

fls <- list.files(BBs_fldr, ".asc$")

print(paste("You have", length(fls), "sequences with successful BBs.", sep=" "))

ids <- str_split_fixed(fls, "_",3)[,1]
ids_unique <- unique(ids)
print(paste("You have", length(ids_unique), "unique individuals with successful BBs.", sep=" "))

#average corridors for individuals w multiple years
for(i in 1:length(ids_unique)){
#for(i in 94:length(ids_unique)){  
  if(i == round(length(ids_unique)/2,0))
    print("You are 1/4 done.")
  if(length(ids[ids == ids_unique[i]])==1){
    bb <- raster(paste0(BBs_fldr, "/", fls[ids == ids_unique[i]]))   # when there is just 1 individual
    cutoff <- sort(values(bb), decreasing=TRUE)
    vlscsum <- cumsum(cutoff)
    cutoff <- cutoff[vlscsum > contour/100][1]
    # bbtemp <- list("Brownian motion variance" = 0, "x" = coordinates(bb)[,1], "y" = coordinates(bb)[,2], "probability" = values(bb))
    # qtl <- bbmm.contour(bbtemp, levels = contour, plot = FALSE)
    bb <- reclassify(bb, rcl=matrix(c(-1,cutoff,0),2,3, byrow=T))
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }else{
    fls2 <- fls[ids == ids_unique[i]]
    bb <- raster(paste0(BBs_fldr, "/", fls2[1]))
    cutoff <- sort(values(bb), decreasing=TRUE)
    vlscsum <- cumsum(cutoff)
    cutoff <- cutoff[vlscsum > contour/100][1]
    # bbtemp <- list("Brownian motion variance" = 0, "x" = coordinates(bb)[,1], "y" = coordinates(bb)[,2], "probability" = values(bb))
    # qtl <- bbmm.contour(bbtemp, levels = contour, plot = FALSE)
    bb <- reclassify(bb, rcl=matrix(c(-1,cutoff,0),2,3, byrow=T))
    bb <- bb/sum(values(bb))   #verify that they add up to 1
    for(e in 2:length(fls2)){
      bb2 <- raster(paste0(BBs_fldr, "/", fls2[e]))
      cutoff <- sort(values(bb2), decreasing=TRUE)
      vlscsum <- cumsum(cutoff)
      cutoff <- cutoff[vlscsum > contour/100][1]
      # bbtemp <- list("Brownian motion variance" = 0, "x" = coordinates(bb2)[,1], "y" = coordinates(bb2)[,2], "probability" = values(bb2))
      # qtl <- bbmm.contour(bbtemp, levels = contour, plot = FALSE)
      bb2 <- reclassify(bb2, rcl=matrix(c(-1,cutoff,0),2,3, byrow=T))
      bb2 <- bb2/sum(values(bb2))   #verify that they add up to 1
      bb <- addLayer(bb, bb2)
    }
    if(nlayers(bb) != length(fls2))
      stop("You have a problem. See error 1.")
    bb <- mean(bb)
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }
  #output averaged individual ASCII file
  m <- as(bb, "SpatialGridDataFrame")
  
  # make sure the cell size is teh same in x and y direction.
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste0(pop_BBs_out_fldr,"/",ids_unique[i],"_ASCII.asc"), attr=1)
  
}

#create overall averages
fls <- list.files(pop_BBs_out_fldr, ".asc$")

if(length(fls)==1){
  bb <- raster(paste0(pop_BBs_out_fldr, "/", fls[1]))
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste0(pop_Rasters_out_fldr,"/","averageUD.asc"), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste0(pop_Rasters_out_fldr,"/","averageUD.img"), format="HFA")
}else{  #if there is at least 1 individual
  bb <- stack(paste0(pop_BBs_out_fldr, "/", fls))
  if(nlayers(bb) != length(fls))
    stop("You have a problem. See error 3.")
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
  write.asciigrid(m, paste0(pop_Rasters_out_fldr,"/","averageUD_migration.asc"), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste0(pop_Rasters_out_fldr,"/","averageUD_migration.img"), format="HFA")
  
}

plot(bb)

# ******************
# 4. Footprints ----
# ******************

fls <- list.files(BBs_fldr, ".asc$")

ids <- str_split_fixed(fls, "_",3)[,1]
ids_unique <- unique(ids)

#this is now for the footprints (without removing the 99% prior to averaging individuals)
for(i in 1:length(ids_unique)){
  if(i == round(length(ids_unique)/2,0))
    print("You are 3/4 done.")
  if(length(ids[ids == ids_unique[i]])==1){
    bb <- raster(paste0(BBs_fldr, "/", fls[ids == ids_unique[i]]))   # when there is just 1 individual
  }else{
    fls2 <- fls[ids == ids_unique[i]]
    bb <- raster(paste0(BBs_fldr, "/", fls2[1]))
    for(e in 2:length(fls2)){
      bb2 <- raster(paste0(BBs_fldr, "/", fls2[e]))
      bb <- addLayer(bb, bb2)
    }
    if(nlayers(bb) != length(fls2))
      stop("You have a problem. See error 1.")
    bb <- mean(bb)
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }
  #99% contours
  cutoff <- sort(values(bb), decreasing=TRUE)
  vlscsum <- cumsum(cutoff)
  cutoff <- cutoff[vlscsum > contour/100][1]
  # bb <- list('Brownian motion variance'=0,x=coordinates(bb)[,1],y=coordinates(bb)[,2],probability=values(bb))
  # contours <- bbmm.contour(bb, levels=contour, plot=F)
  # Create data.frame indicating cells within the each contour and export as Ascii Grid
  contour.99 <- data.frame(x = coordinates(bb)[,1], y = coordinates(bb)[,2], probability = values(bb))
  # contour.99 <- contour.99[contour.99$probability >= contours$Z[1],]
  contour.99$in.out <- ifelse(contour.99$probability >= cutoff, 1, 0)
  #write out footprint for individual
  m <- SpatialPixelsDataFrame(points = contour.99[c("x", "y")], data=contour.99)
  m <- as(m, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste0(pop_footprint_out_fldr,"/",ids_unique[i],"_99pct_contour.asc"), attr=ncol(m))
}

fls <- list.files(pop_BBs_out_fldr, ".asc$")

#create overall population footprint
if(length(fls)==1){
  fls <- list.files(pop_footprint_out_fldr, ".asc$")
  bb <- raster(paste0(pop_footprint_out_fldr, "/", fls[1]))
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste0(pop_footprint_out_fldr,"/","popFootprint.asc"), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste0(pop_footprint_out_fldr,"/","popFootprint.img"), format="HFA")
  print(paste0("End time: ", Sys.time()))
  return("Done. Check your folders. Note you only had 1 ID so there is nothing to average.")
}else{  #if there is at least 1 individual
  bb <- stack(paste0(pop_BBs_out_fldr, "/", fls))
  if(nlayers(bb) != length(fls))
    stop("You have a problem. See error 3.")
  if(ncell(bb)*nlayers(bb) < 10000000){
    bb <- mean(bb)
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }else{
    beginCluster(n=cores)
    bb <- clusterR(bb, fun=calc, list(fun=sum))
    endCluster()
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }
  #create overall population footprint
  fls <- list.files(pop_footprint_out_fldr, ".asc$")
  bb <- stack(paste0(pop_footprint_out_fldr, "/", fls))
  if(nlayers(bb) != length(fls))
    stop("You have a problem. See error 4.")
  if(ncell(bb)*nlayers(bb) < 10000000){
    bb <- sum(bb)
  }else{
    beginCluster(n=cores)
    bb <- clusterR(bb, fun=calc, list(fun=sum))
    endCluster()
  }
  #output averaged individual ASCII file
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste0(pop_footprint_out_fldr,"/","popFootprint.asc"), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste0(pop_footprint_out_fldr,"/","popFootprint.img"), format="HFA")
  print(paste0("End time: ", Sys.time()))
  print("Done. Check your folders.")
}

# ******************************
#5. Metadata ----
# ******************************

## Annual collars deployed ----
# Create metadata for number of collars deployed each nsdYear
yearInds <- data.frame(Year = unique(mt$nsdYear),
                       ActiveCollars = NA)

for(y in 1:nrow(yearInds)){
  yearInds$ActiveCollars[y] <- length(unique(mt$newUid[mt$nsdYear == yearInds$Year[y]])) 
}

write.csv(yearInds, paste(metadata_fldr, "metadata_annualCollars.csv",sep="/"), row.names = F) #save to metadata fldr

## Migration distances ----
mig_dists <- do.call(rbind, lapply(1:length(u), function(i){
  return(data.frame(mig=u[i], max_dist=max(rdist(d[d$mig==u[i],c("x","y")]))/1000))
}))

head(mig_dists)
hist(mig_dists$max_dist)
mig_dists$mean_max_dist <- mean(mig_dists$max_dist)
mig_dists$sd_max_dist <- sd(mig_dists$max_dist)
mig_dists$min_max_dist <- min(mig_dists$max_dist)
mig_dists$max_max_dist <- max(mig_dists$max_dist)
write.csv(mig_dists, file=paste0(metadata_fldr,"/migration_distance_info.csv"), row.names=FALSE)

# ******************************
# 6. Create lines -----
# ******************************

#code to create a spatial lines data frame from a sequences folder
# adapted from WMI Jerod Merkle, 25 Feb 2019

#check the new directories
if(dir.exists(lns_out_fldr)==FALSE){
  dir.create(lns_out_fldr)
}


#load up teh data into a single database
fls <- dir(seqs_fldr)  
print(paste0("You have ", length(fls), " sequences."))
d <- do.call(rbind, lapply(1:length(fls), function(i){
  db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
  db$mig <- sub(".dbf","",fls[i])
  print(nrow(db))
  return(db)
}))


#check and make sure the columns are correct.
all(c("date","x","y") %in% names(d))

d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
d <- d[order(d$mig, d$date),]
coordinates(d) <- c("x","y")
proj4string(d) <- out_proj

u <- unique(d$mig)
lns <- lapply(1:length(u), function(e){
  return(Lines(list(Line(coordinates(d[d$mig==u[e],]))),ID=u[e]))
})
lns <- SpatialLines(lns)
df <- data.frame(mig=u,
                 firstdate=do.call(c, lapply(u, function(e){min(d$date[d$mig==e], na.rm=TRUE)})),
                 lastdate=do.call(c, lapply(u, function(e){max(d$date[d$mig==e], na.rm=TRUE)})))
rownames(df) <- df$mig
df$id <- str_split_fixed(df$mig, "_", 2)[,1]
df$season <- substr(str_split_fixed(df$mig, "_", 2)[,2],1,2)
df$year <- as.numeric(paste0("20",substr(str_split_fixed(df$mig, "_", 2)[,2],3,4)))
lns <- SpatialLinesDataFrame(lns, df)
proj4string(lns) <- out_proj
writeOGR(lns, lns_out_fldr, "migration_lines", driver="ESRI Shapefile")
