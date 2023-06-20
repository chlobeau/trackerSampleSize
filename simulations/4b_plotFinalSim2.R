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

allNewData$sampleSize <- as.numeric(stringr::str_remove(allNewData$population, "pop"))

rare_long2 <- rare_long[rare_long$Contour == "50%" |
                          rare_long$Contour == "75%" |
                          rare_long$Contour == "95%",]
allNewData2 <- allNewData[allNewData$Contour == "50%" |
                            allNewData$Contour == "75%" |
                            allNewData$Contour == "95%",]

## Plots ----

plotCurves <- function(pop, sim){
  
  pred.subset <- pred.df[pred.df$population == pop & pred.df$projName == sim,]
  
  plottedCurve <- ggplot(rare_long2[rare_long2$population == pop & rare_long2$projName == sim,], 
                         aes(x = nAnimals, y = area_sqkm, group = Contour, color = Contour))+
    geom_line(aes(x = nAnimals, y = area_sqkm, group = interaction(simulation, Contour)), color = "red", alpha = .2)+
    geom_line(data = allNewData2[allNewData2$population == pop & allNewData2$projName == sim,], size = 1)+
    
    geom_segment(aes(x = pred.subset$nPredicted[pred.subset$Contour == 50], y = 0, 
                     xend = pred.subset$nPredicted[pred.subset$Contour == 50], 
                     yend = pred.subset$plateau[pred.subset$Contour == 50]), color = "black") +
    geom_segment(aes(x = pred.subset$nPredicted[pred.subset$Contour == 75], y = 0, 
                     xend = pred.subset$nPredicted[pred.subset$Contour == 75], 
                     yend = pred.subset$plateau[pred.subset$Contour == 75]), color = "grey40") +
    geom_segment(aes(x = pred.subset$nPredicted[pred.subset$Contour == 95], y = 0, 
                     xend = pred.subset$nPredicted[pred.subset$Contour == 95], 
                     yend = pred.subset$plateau[pred.subset$Contour == 95]), color = "grey") +
    
    scale_color_manual(values = c("grey", "grey40", "black"))+
    theme_bw()+
    labs(title = "",
         x = "Number of animals collared (n)",
         y = expression(Area~(km^2)))+
    scale_x_continuous(expand = c(0, 0), limits = c(0, as.numeric(stringr::str_remove(pop, "pop")))
                       #                    ,breaks = c(25, 43, 47, 62, 75, 100) # adjust
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 850))
  
  
  return(plottedCurve)
}

wrP250S2 <- plotCurves(pop = "pop250", sim = "sim2")
wrP500S2 <- plotCurves(pop = "pop500", sim = "sim2")
wrP1000S2 <- plotCurves(pop = "pop1000", sim = "sim2")


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

p250s2Map <- mapSimUDs(pop = "pop250", sim = "sim2")
p500s2Map <- mapSimUDs(pop = "pop500", sim = "sim2")
p1000s2Map <- mapSimUDs(pop = "pop1000", sim = "sim2")

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
    labs(title = paste("Population N =", popSize), 
         x = "", y = "")
  
  return(plottedSeqs)
}

p250s2Seqs <- plotSeqs(pop = "pop250", sim = "sim2")
p500s2Seqs <- plotSeqs(pop = "pop500", sim = "sim2")
p1000s2Seqs <- plotSeqs(pop = "pop1000", sim = "sim2")

# Cowplot grid ----
pdf("simulations.pdf", 
    width = 9, height= 10)
plot_grid(p100s2Seqs, p100s2Map, wrP100S2,
          p500s2Seqs, p500s2Map, wrP500S2,
          p1000s2Seqs, p1000s2Map, wrP1000S2, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()

pdf("simulations2.pdf", 
    width = 9, height= 10)
plot_grid(p250s2Seqs, p250s2Map, wrP250S2,
          p500s2Seqs, p500s2Map, wrP500S2,
          p1000s2Seqs, p1000s2Map, wrP1000S2, 
          ncol = 3, rel_widths = c(1,1,2))
dev.off()

load(file = "./sim/predictProjSummaries.rda")
pred.df$population <- factor(pred.df$population, levels = c("pop100", "pop250", "pop500", "pop750", "pop1000"))

require(dplyr)
pred <- pred.df %>% dplyr::arrange(projName, Contour, population) %>%
  subset(population != "pop250" & population != "pop750" &
           Contour != 90 & Contour != 99 &
           projName == "sim2") %>%
  mutate(population = as.numeric(stringr::str_remove(population, "pop"))) %>%
  dplyr::select(-c(projName, UD, nAnimalsActual))

pred$nPredicted[pred$nPredicted > pred$population] <- pred$population[pred$nPredicted > pred$population]

write.csv(pred, "supplementarySimTable.csv", row.names = F)

pred2 <- pred.df %>% dplyr::arrange(projName, Contour, population) %>%
  subset(population != "pop100" & population != "pop750" &
           Contour != 90 & Contour != 99 &
           projName == "sim2") %>%
  mutate(population = as.numeric(stringr::str_remove(population, "pop"))) %>%
  dplyr::select(-c(projName, UD, nAnimalsActual))

pred2$nPredicted[pred2$nPredicted > pred2$population] <- pred2$population[pred2$nPredicted > pred2$population]

# look at all new data area estimates at variance for those nPredicted values

pred2$estMean <- NA
pred2$estSD <- NA
pred2$areaActual <- NA

for(i in 1:nrow(pred2)){
  sub <- rare_long[rare_long$projName == "sim2" &
                     rare_long$population == paste0("pop", pred2$population[i]) &
                     rare_long$Contour == paste0(pred2$Contour[i], "%") &
                     rare_long$nAnimals == pred2$nPredicted[i],]
  pred2$estMean[i] <- mean(sub$area_sqkm)
  pred2$estSD[i] <- sd(sub$area_sqkm)
  pred2$areaActual[i] <- max(sub$area_sqkm)
}

# summary <- rare_long[rare_long$projName == "sim2" &
#                        rare_long$population == "pop250" &
#                        rare_long$Contour == "95%",] %>% 
#   group_by(nAnimals) %>%
#   summarise(areaMean = mean(area_sqkm),
#             areaSD = sd(area_sqkm))

write.csv(pred2, "supplementarySimTable2.csv", row.names = F)

