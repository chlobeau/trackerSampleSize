#Load in saved blank rasters (made w 0_EmptyPixelStateExtent.R)
# when running WMI BBMM scripts, outside the application use this workflow to create a grid to calculate BBs over
# then find/replace grd in the WMI scripts with grd.studyarea generated here

# ***********
# Inputs ----
# ***********

blnkrstr_fldr <- "E:/migMap/mapperCode/BlankRasters"

cell.size #from step 1, be sure it matches what you want
##load blank raster below based on cell.size
load(paste0(blnkrstr_fldr,"/", "BlankRasterExtent_CO_UT_50m_20220304.rda"))
#load(paste0(blnkrstr_fldr,"/", "BlankRasterExtent_CO_UT_100m_20220304.rda"))
#load(paste0(blnkrstr_fldr,"/", "BlankRasterExtent_CO_UT_250m_20220304.rda"))
#load(paste0(blnkrstr_fldr,"/", "BlankRasterExtent_CO_UT_500m_20220304.rda"))
#rename blank raster loaded above to grd.bigempty
grd.bigempty <- grd.agg50m

#Create smaller extent of particular study area to crop from larger grid
#read in data, re-project 

d <- st_read(shpfl_fldr, shpfl_name)
d <- st_transform(d, out_proj)

st_crs(d)==st_crs(grd.bigempty) #check that crs match

#extend shapefile extent, then crop from statewide raster
ext.studyarea <- raster::extent(d)
multiplyers.studyarea <- c((ext.studyarea[2]-ext.studyarea[1])*mult4buff, (ext.studyarea[4]-ext.studyarea[3])*mult4buff)   # add about 20% around the edges of your extent (you can adjust this if necessary)
ext.studyarea <- raster::extend(ext.studyarea, multiplyers.studyarea)
grd.studyarea <- crop(grd.bigempty, raster::extent(ext.studyarea))

origin(grd.bigempty)
origin(grd.studyarea)
extent(d)
extent(grd.studyarea)
mosaic(grd.bigempty, grd.studyarea, fun = sum) #check that they play nicely together

save(grd.studyarea, file = paste0(NewProjectFolder,"/","studyarea_grd.rda"))
