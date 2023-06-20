# adapted from WMI Migration Mapper v2_3
# Wyoming Migration Initiative. 2017. Migration Mapper. University of Wyoming, Laramie, Wyoming.

# ***********
# Inputs ----
# ***********

# NewProjectFolder and BB settings from step 1

mig.metadata.file=paste0(NewProjectFolder, "/metadata/metadata_migration.csv") #need to run migration bbs first

seqs_fldr = paste0(NewProjectFolder, "/sequencesWinter")
metadata_fldr=paste0(NewProjectFolder, "/metadata")

BBs_fldr = paste0(NewProjectFolder, "/UDsWinter")
pop_BBs_out_fldr = paste0(NewProjectFolder, "/UDs_popWinter")
pop_Rasters_out_fldr = paste0(NewProjectFolder, "/popRasters")

#load study area grd that aligns with larger raster extent
load(paste0(NewProjectFolder, "/studyarea_grd.rda"))
cell.size # from step 1 make sure this matches studyarea_grd resolution
res(grd.studyarea)

# *************************
# 1. Winter Sequences ----
# *************************
## based on C:\MAPP_2.3\app\scriptsv2\create.seqs.W.R
# This function takes the outputs from Migration Mapper tab 6 and creates
# dbf files of each individual animal's migration sequences

fl <- read.csv(mig.metadata.file)

tempYr<-fl$input.file ##jgage
tempYr<-as.character(tempYr) ##jgage
tempYr<-substr(tempYr,nchar(tempYr)-1,nchar(tempYr)) ##jgage
tempSeas<-fl$input.file ##jgage
tempSeas<-as.character(tempSeas) ##jgage
tempSeas<-substr(tempSeas,nchar(tempSeas)-3,nchar(tempSeas)-2) ##jgage

fl$yr<-tempYr ##jgage
fl$seas<-tempSeas ##jgage
fl$yr <- as.numeric(ifelse(as.numeric(fl$yr) > 80, paste0("19",fl$yr), paste0("20",fl$yr)))
fl$Start.Date <- as.POSIXct(strptime(fl$Start.Date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
fl$End.Date <- as.POSIXct(strptime(fl$End.Date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")

fl$id <- str_split_fixed(fl$input.file, "_",2)[,1]

# read in the shapefile

d <- st_read(shpfl_fldr, shpfl_name)

d <- as(d, "Spatial")

#reduce the dataset to columns of interest
d <- d[,c("newUid", "nwMstrD")]
names(d) <- c("id","date")

#fix the dates  (it is OK to specify GMT, since all dates will be in GMT!)
d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
d$id <- as.character(d$id)
d$year <- as.numeric(strftime(d$date, format = "%Y", tz = "GMT"))

#build a database of the potential winter periods

ids <- unique(d$id)
dts <- do.call(rbind, lapply(1:length(ids), function(i){
  yr_rng <- sort(unique(d$year[d$id == ids[i]]))
  tmp2 <- do.call(rbind, lapply(0:length(yr_rng), function(e){
    if(e == 0){
      wint.end <- fl$Start.Date[fl$id==ids[i] & fl$seas == "sp" & fl$yr == yr_rng[e+1]]
      if(length(wint.end)==0){
        wint.end <- NA
      }
      return(data.frame(id=ids[i], wintr=paste(yr_rng[e+1]-1, yr_rng[e+1], sep="_"),
                        id_yr=paste(ids[i],substr(yr_rng[e+1],3,4),sep="_"), wint.start=NA,
                        wint.end=as.character(wint.end-86400)))
    }else{
      if(e == length(yr_rng)){
        wint.start <- fl$End.Date[fl$id==ids[i] & fl$seas == "fa" & fl$yr == yr_rng[e]]
        if(length(wint.start)==0){
          wint.start <- NA
        }
        return(data.frame(id=ids[i], wintr=paste(yr_rng[e], yr_rng[e]+1, sep="_"),
                          id_yr=paste(ids[i],substr(yr_rng[e]+1,3,4),sep="_"),
                          wint.start=as.character(wint.start+86400),
                          wint.end=NA))
      }else{
        wint.start <- fl$End.Date[fl$id==ids[i] & fl$seas == "fa" & fl$yr == yr_rng[e]]
        wint.end <- fl$Start.Date[fl$id==ids[i] & fl$seas == "sp" & fl$yr == yr_rng[e+1]]
        if(length(wint.end)==0){
          wint.end <- NA
        }
        if(length(wint.start)==0){
          wint.start <- NA
        }
        return(data.frame(id=ids[i], wintr=paste(yr_rng[e], yr_rng[e+1], sep="_"),
                          id_yr=paste(ids[i],substr(yr_rng[e+1],3,4),sep="_"),
                          wint.start=as.character(wint.start+86400),
                          wint.end=as.character(wint.end-86400)))
      }
    }
  }))
  return(tmp2)
}))

dts$wint.start <- as.character(dts$wint.start) #fix the dates  (it is OK to specify GMT, since all dates will be in GMT!)
dts$wint.end <- as.character(dts$wint.end)

mnths <- as.numeric(substr(na.omit(dts$wint.start), 6,7))
mn_strt <- mean(as.Date(strptime(paste0(ifelse(mnths > 8, "2017","2018"),"-",substr(na.omit(as.character(dts$wint.start)),6,10)),format = "%Y-%m-%d")), na.rm = T)

## BIGHORN - Statewide mean winter start ----
# Oct 30
#mn_strt <- as.Date("2017-10-30")

if(as.numeric(substr(mn_strt,6,7)) > 8){
  yrs <- substr(dts$wintr[is.na(dts$wint.start)==TRUE],1,4)
}else{
  yrs <- substr(dts$wintr[is.na(dts$wint.start)==TRUE],6,9)
}

dts$wint.start[is.na(dts$wint.start)==TRUE] <- paste0(yrs, "-", substr(mn_strt,6,10), " 23:00:00")

mnths <- as.numeric(substr(na.omit(dts$wint.end), 6,7))
mn_end <- mean(as.Date(strptime(paste0(ifelse(mnths > 8, "2017","2018"),"-",substr(na.omit(as.character(dts$wint.end)),6,10)),format = "%Y-%m-%d")), na.rm = T) #leap year in 2016, rendered this all na, need to add na.rm = T
#mn_end <- mean(as.Date(strptime(paste0(ifelse(mnths > 8, "2017","2018"),"-",substr(na.omit(as.character(dts$wint.end)),6,10)),format = "%Y-%m-%d"))) #leap year in 2016, rendered this all na, need to add na.rm = T

## BIGHORN - Statewide mean winter end ----
# May 27
#mn_end <- as.Date("2018-05-27")

if(as.numeric(substr(mn_end,6,7)) > 8){
  yrs <- substr(dts$wintr[is.na(dts$wint.end)==TRUE],1,4)
}else{
  yrs <- substr(dts$wintr[is.na(dts$wint.end)==TRUE],6,9)
}

dts$wint.end[is.na(dts$wint.end)==TRUE] <- paste0(yrs, "-", substr(mn_end,6,10), " 01:00:00")

dts$wint.start <- as.POSIXct(strptime(dts$wint.start,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
dts$wint.end <- as.POSIXct(strptime(dts$wint.end,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")

d$id <- as.character(d$id)
dts$id <- as.character(dts$id)
dts$numb_of_pts <- NA

for(i in 1:nrow(dts)){
  tmp <- d[d$date >= dts$wint.start[i] & d$date <= dts$wint.end[i] & d$id == dts$id[i],]
  tmp <- cbind(as.data.frame(tmp)[,c("id","date")],coordinates(tmp))    # get it out of sp object
  tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
  names(tmp) <- c("id","date","x","y")   # rename columns
  dts$numb_of_pts[i] <- nrow(tmp)
  if(nrow(tmp)>0){
    write.dbf(tmp, file = paste0(seqs_fldr,"/",sub("_","_wi",dts$id_yr[i]),".dbf"))    #write as dbf
  }
  else{
    next
  }
}

# *********************
# 2. Individual WR ----
# *********************
## based on C:\MAPP_2.3\app\scriptsv2\create.BBs.W.R
# This function takes a folder of migration sequences (in .dbf format)
# and creates regular brownian bridges from them

#create database for files
fls <- list.files(seqs_fldr, ".dbf$")
print(paste0("You have ", length(fls), " sequences."))
d <- do.call(rbind, lapply(1:length(fls), function(i){
  db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
  db$wint <- sub(".dbf","",fls[i])
  return(db)
}))

#loop through all individuals create asciis
u <- unique(d$wint)

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
#migtime<-migtime[,c('id_yr','newUid','nsdYear','notes')] ##jgage

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

regBB <- list()
Cores <- detectCores()
cl <- makeCluster(Cores - 4)
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
                   if(length(unique(jul)) < mindays){ #filters out animals with less than mindays of data
                     BB_meta <- data.frame(input.file=u[i],
                                           brownian.motion.variance=NA,
                                           grid.size=NA,
                                           grid.cell.size=NA,
                                           date.created=Sys.time(),
                                           execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                                           num.locs=nrow(temp),
                                           Start.Date=min(temp$date),
                                           End.Date=max(temp$date),
                                           num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                                           errors=paste0("Less than ",mindays, " days of data."),
                                           note1="")
                   } 
                   
                   else{ #also filter out animals with only 1 pt per day on average
                     if(nrow(temp) <= length(unique(jul))){
                       BB_meta <- data.frame(input.file=u[i],
                                             brownian.motion.variance=NA,
                                             grid.size=NA,
                                             grid.cell.size=NA,
                                             date.created=Sys.time(),
                                             execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                                             num.locs=nrow(temp),
                                             Start.Date=min(temp$date),
                                             End.Date=max(temp$date),
                                             num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                                             errors="Only 1 point per day, on average.",
                                             note1="")#jgage
                     } 
                     if(nrow(temp) < 4){ #also filter out animals w <4 pts
                       BB_meta <- data.frame(input.file=u[i],
                                             brownian.motion.variance=NA,
                                             grid.size=NA,
                                             grid.cell.size=NA,
                                             date.created=Sys.time(),
                                             execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                                             num.locs=nrow(temp),
                                             Start.Date=min(temp$date),
                                             End.Date=max(temp$date),
                                             num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                                             errors="Less than 4 points.",
                                             note1="")
                     } 
                     else{
                       #prepare only the cells to run BB over
                       ext2 <- raster::extent(temp)
                       multiplyers <- c((ext2[2]-ext2[1])*mult4buff, (ext2[4]-ext2[3])*mult4buff)   # add about mult4buff around the edges of your extent (you can adjust this if necessary)
                       ext2 <- raster::extend(ext2, multiplyers)
                       cels <- raster::cellsFromExtent(grd.studyarea, ext2)
                       
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
                       write.asciigrid(m, paste(BBs_fldr,"/",u[i],"_ASCII.asc",sep=""), attr=3)
                       
                       BB_meta <- data.frame(input.file=u[i],
                                             brownian.motion.variance=round(bb[[1]],2),
                                             grid.size=length(bb$x),
                                             grid.cell.size=abs(bb$x[1]-bb$x[2]),
                                             date.created=Sys.time(),
                                             execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                                             num.locs=nrow(temp),
                                             Start.Date=min(temp$date),
                                             End.Date=max(temp$date),
                                             num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                                             errors="None",
                                             note1="")
                     }
                   }
                   regBB[[i]] <- BB_meta
                 }

stopCluster(cl)
write.csv(indBB, file=paste(metadata_fldr,"/metadata_winter.csv",sep=""), row.names=FALSE)

# *********************
# 3. Population WR ----
# *********************
## based on C:\MAPP_2.3\app\scriptsv2\create.BB.avgs.W.R
# This function takes in individual BBs and creates averaged individual BBs
# and population UDs

fls <- list.files(BBs_fldr, ".asc$")

print(paste("You have ", length(fls), " sequences with successful Winter Range BBs.", sep=""))

ids <- str_split_fixed(fls, "_",3)[,1]

ids_unique <- unique(ids)
print(paste("You have ", length(ids_unique), " unique individuals with successful Winter Range BBs.", sep=""))

#average individual winter ranges over years
for(i in 1:length(ids_unique)){
  if(i == round(length(ids_unique)/2,0))
    print("You are half done.")
  if(length(ids[ids == ids_unique[i]])==1){
    bb <- raster(paste(BBs_fldr, "/", fls[ids == ids_unique[i]], sep=""))   # when there is just 1 individual
  }else{
    fls2 <- fls[ids == ids_unique[i]]
    bb <- raster(paste(BBs_fldr, "/", fls2[1], sep=""))
    for(e in 2:length(fls2)){
      bb <- addLayer(bb, raster(paste(BBs_fldr, "/", fls2[e], sep="")))
    }
    if(nlayers(bb) != length(fls2))
      stop("You have a problem. See error 1.")
    bb <- mean(bb)
    bb <- bb/sum(values(bb))   #verify that they add up to 1
  }
  
  #output averaged individual ASCII file
  m <- as(bb, "SpatialGridDataFrame")
  
  # make sure the cell size is the same in x and y direction.
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste(pop_BBs_out_fldr,"/",ids_unique[i],"_ASCII.asc",sep=""), attr=1)
  
}

#create overall averages
fls <- list.files(pop_BBs_out_fldr, ".asc$")

if(length(fls)==1){
  bb <- raster(paste(pop_BBs_out_fldr, "/", fls[1], sep=""))
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize")
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2)
  write.asciigrid(m, paste(pop_Rasters_out_fldr,"/","averageUD_winter.asc",sep=""), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste(pop_Rasters_out_fldr,"/","averageUD_winter.img",sep=""), format="HFA")
  print(paste0("End time: ", Sys.time()))
  print("Done. Check your folders. Note that you only have 1 individual, so no averaging occured!")
}else{                                                                      # THIS IS WHAT I'LL WANT IN RAREFACTION CURVE
  bb <- stack(paste(pop_BBs_out_fldr, "/", fls, sep=""))
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
  write.asciigrid(m, paste(pop_Rasters_out_fldr,"/","averageUD_winter.asc",sep=""), attr=1)
  projection(bb) <- out_proj
  writeRaster(bb, filename = paste(pop_Rasters_out_fldr,"/","averageUD_winter.img",sep=""), format="HFA") #change this fldr location to population loc for SR, WR, migs
  print(paste0("End time: ", Sys.time()))
}

plot(bb)
