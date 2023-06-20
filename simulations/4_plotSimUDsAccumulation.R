require(mapview)
require(raster)
require(spatialEco)
require(ggplot2)
require(cowplot)

# Cumulative curves ----
load(file = paste0(getwd(), "/sim/RarefactionLONGWCovs.rda")) #load simulated rarefaction data
load(file = paste0(getwd(), "/simPredictedData.rda")) #load new data

#make contour a factor and add "%"
allNewData$Contour <- paste0(allNewData$Contour, "%")
allNewData$Contour <- factor(allNewData$Contour, levels = c("99%", "95%", "90%", "75%", "50%"))

rare_long$Contour <- paste0(rare_long$Contour, "%")
rare_long$Contour <- factor(rare_long$Contour, levels = c("99%", "95%", "90%", "75%", "50%"))

rare_long$projName <- factor(rare_long$projName, levels = c("sim1", "sim2", "sim3"))
allNewData$projName <- factor(allNewData$projName, levels = c("sim1", "sim2", "sim3"))

unique(rare_long$projName)
unique(allNewData$projName)

#for allNewData make column in data frame for whether nanimals are > actual

allNewData$sampleSize <- 100
allNewData$predicted <- 0
allNewData$predicted[allNewData$nAnimals > allNewData$sampleSize] <- 1

rare_long2 <- rare_long[rare_long$Contour == "50%" |
                          rare_long$Contour == "75%" |
                          rare_long$Contour == "95%",]
allNewData2 <- allNewData[allNewData$Contour == "50%" |
                            allNewData$Contour == "75%" |
                            allNewData$Contour == "95%",]

## Plots ----

plotCurves <- function(pop, sim){
  
  plottedCurve <- ggplot(rare_long2[rare_long2$population == pop & rare_long2$projName == sim,], 
                         aes(x = nAnimals, y = area_sqkm, group = Contour, color = Contour))+
    geom_line(aes(x = nAnimals, y = area_sqkm, group = interaction(simulation, Contour)), color = "red", alpha = .2)+
    geom_line(data = allNewData2[allNewData2$population == pop & allNewData2$projName == sim & allNewData2$predicted == 1,], aes(linetype='Extrapolated'))+
    geom_line(data = allNewData2[allNewData2$population == pop & allNewData2$projName == sim & allNewData2$predicted == 0,], aes(linetype='Empirical'), size = 1)+
    
    scale_color_manual(values = c("grey", "grey40", "black"))+
    scale_linetype_manual('', values=c(1,3))+   
    guides(color = guide_legend(order = 1),
           size = guide_legend(order = 2))+
    theme_bw()+
    labs(title = "",
         x = "Number of animals collared (n)",
         y = expression(Area~(km^2)))+
    scale_x_continuous(expand = c(0, 0), limits = c(0, 150)
                       #                    ,breaks = c(25, 43, 47, 62, 75, 100) # adjust
    ) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  
  
  return(plottedCurve)
}

wrP100S1 <- plotCurves(pop = "pop100", sim = "sim1")
wrP100S2 <- plotCurves(pop = "pop100", sim = "sim2")
wrP100S3 <- plotCurves(pop = "pop100", sim = "sim3")

wrP250S1 <- plotCurves(pop = "pop250", sim = "sim1")
wrP250S2 <- plotCurves(pop = "pop250", sim = "sim2")
wrP250S3 <- plotCurves(pop = "pop250", sim = "sim3")

wrP500S1 <- plotCurves(pop = "pop500", sim = "sim1")
wrP500S2 <- plotCurves(pop = "pop500", sim = "sim2")
wrP500S3 <- plotCurves(pop = "pop500", sim = "sim3")

wrP750S1 <- plotCurves(pop = "pop750", sim = "sim1")
wrP750S2 <- plotCurves(pop = "pop750", sim = "sim2")
wrP750S3 <- plotCurves(pop = "pop750", sim = "sim3")

wrP1000S1 <- plotCurves(pop = "pop1000", sim = "sim1")
wrP1000S2 <- plotCurves(pop = "pop1000", sim = "sim2")
wrP1000S3 <- plotCurves(pop = "pop1000", sim = "sim3")


# UDs @ varying contours ----

projFolders <- paste0(getwd(),"/sim/projFolders")
out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

mapSimUDs <- function(pop, sim){
  s1_fldr <- paste(projFolders, pop, sim, "UDsWinter", sep = "/")
  fls <- list.files(s1_fldr, ".asc$")
  
  s1 <- stack(paste(s1_fldr, fls, sep = "/"))
  
  s1 <- sum(s1)
  s1 <- s1/sum(values(s1))
  
  projection(s1) <- out_proj
  s1 <- terra::rast(s1)
  
  s195.rasterized <- raster.vol(s1, p = 0.95)
  s195.rasterized <- raster(s195.rasterized)
  polygonized.s195 <- rasterToPolygons(s195.rasterized, dissolve = TRUE)#convert raster 
  poly.s195 <- polygonized.s195[polygonized.s195$layer == 1,] #subset by just the home range polygons and ignore the outer perimter "box" polygon
  
  s175.rasterized <- raster.vol(s1, p = 0.75)
  s175.rasterized <- raster(s175.rasterized)
  polygonized.s175 <- rasterToPolygons(s175.rasterized, dissolve = TRUE)#convert raster 
  poly.s175 <- polygonized.s175[polygonized.s175$layer == 1,] #subset by just the home range polygons and ignore the outer perimter "box" polygon
  
  s150.rasterized <- raster.vol(s1, p = 0.50)
  s150.rasterized <- raster(s150.rasterized)
  polygonized.s150 <- rasterToPolygons(s150.rasterized, dissolve = TRUE)#convert raster 
  poly.s150 <- polygonized.s150[polygonized.s150$layer == 1,] #subset by just the home range polygons and ignore the outer perimter "box" polygon
  
  ### plot
  s1[s1 == 0] <- NA
  
  s1.spdf <- as(raster(s1), "SpatialPixelsDataFrame")
  s1.df <- as.data.frame(s1.spdf)
  colnames(s1.df) <- c("value", "x", "y")
  
  Map <- ggplot() +  
    geom_tile(data=s1.df, aes(x=x, y=y, fill=value), alpha=0.8)+ 
    geom_polygon(data=poly.s150, aes(x=long, y=lat, group=group), 
                 fill=NA, color="black", size=0.25)+
    geom_polygon(data=poly.s175, aes(x=long, y=lat, group=group), 
                 fill=NA, color="gray35", size=0.25)+
    geom_polygon(data=poly.s195, aes(x=long, y=lat, group=group),
                 fill=NA, color="gray", size=0.25)+
    scale_fill_gradient(low = "#F5F1ED", high = "black")+
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position = "none")+
    labs(title = "", x = "", y = "")
  
  return(Map)
}

p100s1Map <- mapSimUDs(pop = "pop100", sim = "sim1")
p100s2Map <- mapSimUDs(pop = "pop100", sim = "sim2")
p100s3Map <- mapSimUDs(pop = "pop100", sim = "sim3")

p250s1Map <- mapSimUDs(pop = "pop250", sim = "sim1")
p250s2Map <- mapSimUDs(pop = "pop250", sim = "sim2")
p250s3Map <- mapSimUDs(pop = "pop250", sim = "sim3")

p500s1Map <- mapSimUDs(pop = "pop500", sim = "sim1")
p500s2Map <- mapSimUDs(pop = "pop500", sim = "sim2")
p500s3Map <- mapSimUDs(pop = "pop500", sim = "sim3")

p750s1Map <- mapSimUDs(pop = "pop750", sim = "sim1")
p750s2Map <- mapSimUDs(pop = "pop750", sim = "sim2")
p750s3Map <- mapSimUDs(pop = "pop750", sim = "sim3")

p1000s1Map <- mapSimUDs(pop = "pop1000", sim = "sim1")
p1000s2Map <- mapSimUDs(pop = "pop1000", sim = "sim2")
p1000s3Map <- mapSimUDs(pop = "pop1000", sim = "sim3")

# Sequences ----

plotSeqs <- function(pop, sim){
  seqs_fldr <- paste(projFolders, pop, sim, "sequencesWinter", sep = "/")
  fls <- list.files(seqs_fldr, ".dbf$")
  
  d <- do.call(rbind, lapply(1:length(fls), function(i){
    db <- foreign::read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
    db$wint <- sub(".dbf","",fls[i])
    return(db)
  }))
  
  simNumber <- substr(sim, nchar(sim), nchar(sim))
  popSize <- gsub("pop", "", pop)
  
  plottedSeqs <- ggplot() +  
    geom_point(data=d, aes(x=x, y=y, color = id), size = .5, alpha = .5)+
    scale_color_manual(values = rainbow(popSize))+
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position = "none",
          plot.title = element_text(size = 10))+
    labs(title = paste("Population", popSize, "Simulation", simNumber), 
         x = "", y = "")
  
  return(plottedSeqs)
}

p100s1Seqs <- plotSeqs(pop = "pop100", sim = "sim1")
p100s2Seqs <- plotSeqs(pop = "pop100", sim = "sim2")
p100s3Seqs <- plotSeqs(pop = "pop100", sim = "sim3")

p250s1Seqs <- plotSeqs(pop = "pop250", sim = "sim1")
p250s2Seqs <- plotSeqs(pop = "pop250", sim = "sim2")
p250s3Seqs <- plotSeqs(pop = "pop250", sim = "sim3")

p500s1Seqs <- plotSeqs(pop = "pop500", sim = "sim1")
p500s2Seqs <- plotSeqs(pop = "pop500", sim = "sim2")
p500s3Seqs <- plotSeqs(pop = "pop500", sim = "sim3")

p750s1Seqs <- plotSeqs(pop = "pop750", sim = "sim1")
p750s2Seqs <- plotSeqs(pop = "pop750", sim = "sim2")
p750s3Seqs <- plotSeqs(pop = "pop750", sim = "sim3")

p1000s1Seqs <- plotSeqs(pop = "pop1000", sim = "sim1")
p1000s2Seqs <- plotSeqs(pop = "pop1000", sim = "sim2")
p1000s3Seqs <- plotSeqs(pop = "pop1000", sim = "sim3")

# Cowplot grid ----
pdf("simulationsP100b.pdf", 
    width = 9, height= 12)
plot_grid(p100s1Seqs, p100s1Map, wrP100S1,
          p100s2Seqs, p100s2Map, wrP100S2,
          p100s3Seqs, p100s3Map, wrP100S3, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()

pdf("simulationsP250.pdf", 
    width = 9, height= 12)
plot_grid(p250s1Seqs, p250s1Map, wrP250S1,
          p250s2Seqs, p250s2Map, wrP250S2,
          p250s3Seqs, p250s3Map, wrP250S3, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()

pdf("simulationsP500.pdf", 
    width = 9, height= 12)
plot_grid(p500s1Seqs, p500s1Map, wrP500S1,
          p500s2Seqs, p500s2Map, wrP500S2,
          p500s3Seqs, p500s3Map, wrP500S3, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()

pdf("simulationsP750.pdf", 
    width = 9, height= 12)
plot_grid(p750s1Seqs, p750s1Map, wrP750S1,
          p750s2Seqs, p750s2Map, wrP750S2,
          p750s3Seqs, p750s3Map, wrP750S3, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()

pdf("simulationsP1000.pdf", 
    width = 9, height= 12)
plot_grid(p1000s1Seqs, p1000s1Map, wrP1000S1,
          p1000s2Seqs, p1000s2Map, wrP1000S2,
          p1000s3Seqs, p1000s3Map, wrP1000S3, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()
