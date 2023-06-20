library(nlme)
library(aomisc)
library(ggplot2)
library(drc)

load(file = paste0(getwd(), "/sim/RarefactionLONGWCovs.rda"))

# *********
# LOOP ----
# *********

ContourNew <- unique(rare_long$Contour)

pred.df <- data.frame()
asymptote.df <- data.frame(unique(rare_long[,c("projName", "UD")]))

for(c in 1:length(ContourNew)){
  asymptote.df$Contour  <-  ContourNew[c]
  pred.df <- rbind(pred.df, asymptote.df)
}
pred.df$nPredicted <- NA
pred.df$plateau <- NA
pred.df$nAnimalsActual <- NA

seasons <- unique(pred.df$UD)
allList <- list()
for(s in 1:length(seasons)){
  UD_sub <- rare_long[rare_long$UD == seasons[s],]
  UD_projects <- unique(UD_sub$projName)
  seasonList <- list()
  for(p in 1:length(UD_projects)){
    projUD_sub <- UD_sub[UD_sub$projName == UD_projects[p],]
    projectList <- list()
    for(c in 1:length(ContourNew)){
      #group the data by simulation
      sppUDgrp <- groupedData(area_sqkm ~ nAnimals | simulation, 
                              data=projUD_sub[projUD_sub$Contour == ContourNew[c],
                                              c("nAnimals", "area_sqkm", "simulation", "Contour")])
      #cycle through all groups
      
      mod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                 data = projUD_sub[projUD_sub$Contour == ContourNew[c],])
      
      newdata <- expand.grid(nAnimals = 1:4000, Contour = ContourNew[c])
      newdata$area_sqkm <- predict(mod, newdata = newdata)
      
      #save new data to list to plot
      projectList[[c]] <- newdata
      
      # pull out sample sizes
      ypredict <- as.numeric(coef(mod)[names(coef(mod)) == "plateau:(Intercept)"])
      #ypredict <- as.numeric(coef(mod.nlme)[names(coef(mod.nlme)) == "plateau"][1,]) #for nlme
      nPredicted <- min(newdata$nAnimals[newdata$area_sqkm >= (ypredict - (ypredict * .01))])
      
      #add predictions to dataframe
      pred.df$nPredicted[pred.df$UD == seasons[s] &
                           pred.df$projName == UD_projects[p] &
                           pred.df$Contour == ContourNew[c]
                           ] <- nPredicted
      pred.df$plateau[pred.df$UD == seasons[s] &
                        pred.df$projName == UD_projects[p] &
                        pred.df$Contour == ContourNew[c]
                      ] <- ypredict
      pred.df$nAnimalsActual[pred.df$UD == seasons[s] &
                               pred.df$projName == UD_projects[p] &
                               pred.df$Contour == ContourNew[c]
                             ] <- max(projUD_sub$nAnimals[projUD_sub$Contour == ContourNew[c]])
      #add spp to dataframe
      pred.df$spp[pred.df$UD == seasons[s] &
                    pred.df$projName == UD_projects[p] &
                    pred.df$Contour == ContourNew[c]] <- projUD_sub$spp2[1]
    }
    newdata2 <- do.call("rbind", projectList)
    newdata2$projName <- UD_projects[p]
    newdata2$spp <- projUD_sub$spp2[1]
    
    ##plot each project UD
    pl  <-  ggplot(projUD_sub, aes(nAnimals, area_sqkm, color = as.factor(Contour))) +
      geom_point()+
      geom_line(data = newdata2, aes(group = as.factor(Contour)), color = "black")+
      xlim(c(0,200))+
      labs(title = paste(UD_projects[p], seasons[s]))
    plot(pl)
    seasonList[[p]] <- newdata2
  }
  seas.df <- do.call("rbind", seasonList)
  seas.df$UD <- seasons[s]
  allList[[s]] <- seas.df
}
allNewData <- do.call("rbind", allList)

save(allNewData, file = "simPredictedData.rda")
# 
pred.df <- pred.df[with(pred.df, order(projName, UD, Contour)),]
save(pred.df, file = "predictProjSummaries.rda")