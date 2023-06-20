# This code creates a large empty pixel grid (greater than the state of CO or greater than the state of CO and UT)
# to be able to crop out smaller study area extents while maintaining a common origin to resample/combine data statewide
# because of the large extent we need to use projection crs NAD83 / Conus Albers 
# lat long wouldn't convert nicely to a grid

library(sf)
library(maps) 
library(raster)

#CO <- map("state", "Colorado")
CO_UT <- map("state", c("Colorado", "Utah"))

# input settings, should match WMI settings 
mult4buff=.2
cell.size=50 #50m, 250m, 500m
out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

#Create a blank grid that extends beyond state borders to clip out smaller study areas
CO_UT <- st_as_sf(CO_UT)
CO_UT <- st_transform(CO_UT, out_proj)
ext.bigempty <- raster::extent(CO_UT)
multiplyers.bigempty <- c((ext.bigempty[2]-ext.bigempty[1])*mult4buff, (ext.bigempty[4]-ext.bigempty[3])*mult4buff)   # add about 20% around the edges of your extent (you can adjust this if necessary)
ext.bigempty <- raster::extend(ext.bigempty, multiplyers.bigempty)
grd.bigempty <- raster(ext=ext.bigempty,
                       res=cell.size,
                       crs = out_proj)

#extent of 500m pixel raster was slightly different which probably isn't a big deal but I'm going to be a perfectionist about it
#use aggregate/disaggregate function to resize pixels to ensure extent is the same
grd.agg500m <- aggregate(grd.bigempty, fact=10)
origin(grd.agg500m)
extent(grd.agg500m)
grd.agg250m <- disaggregate(grd.agg, fact=2)
origin(grd.agg250m)
extent(grd.agg250m)
grd.agg50m <- disaggregate(grd.agg, fact=10)
extent(grd.agg50m)
origin(grd.agg50m)

save(grd.agg50m, file = "BlankRasterExtent_CO_UT_50m_20220304.rda")
save(grd.agg250m, file = "BlankRasterExtent_CO_UT_250m_20220304.rda")
save(grd.agg500m, file = "BlankRasterExtent_CO_UT_500m_20220304.rda")

