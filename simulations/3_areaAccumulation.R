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
asymptote.df <- data.frame(unique(rare_long[,c("population", "projName", "UD")]))

for(c in 1:length(ContourNew)){
  asymptote.df$Contour  <-  ContourNew[c]
  pred.df <- rbind(pred.df, asymptote.df)
}
pred.df$nPredicted <- NA
pred.df$plateau <- NA
pred.df$nAnimalsActual <- NA


allList <- list()

# populations <- unique(pred.df$population)
# for(p in 1:length(populations)){
#   pop_sub <- rare_long[rare_long$population == populations[p],]
#   pop_projects <- unique(pop_sub$projName)
#   
# }

populations <- unique(pred.df$population)
for(s in 1:length(populations)){
  pop_sub <- rare_long[rare_long$population == populations[s],]
  pop_projects <- unique(pop_sub$projName)
  seasonList <- list()
  for(p in 1:length(pop_projects)){
    projUD_sub <- pop_sub[pop_sub$projName == pop_projects[p],]
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
      pred.df$nPredicted[pred.df$UD == "winter" &
                           pred.df$population == populations[s] &
                           pred.df$projName == pop_projects[p] &
                           pred.df$Contour == ContourNew[c]
                           ] <- nPredicted
      pred.df$plateau[pred.df$UD == "winter" &
                        pred.df$population == populations[s] &
                        pred.df$projName == pop_projects[p] &
                        pred.df$Contour == ContourNew[c]
                      ] <- ypredict
      pred.df$nAnimalsActual[pred.df$UD == "winter" &
                               pred.df$projName == pop_projects[p] &
                               pred.df$Contour == ContourNew[c]
                             ] <- max(projUD_sub$nAnimals[projUD_sub$Contour == ContourNew[c]])
      #add spp to dataframe
      pred.df$spp[pred.df$UD == "winter" &
                    pred.df$projName == pop_projects[p] &
                    pred.df$Contour == ContourNew[c]] <- projUD_sub$spp2[1]
    }
    newdata2 <- do.call("rbind", projectList)
    newdata2$projName <- pop_projects[p]
    newdata2$spp <- projUD_sub$spp2[1]
    
    # ##plot each project UD
    # pl  <-  ggplot(projUD_sub, aes(nAnimals, area_sqkm, color = as.factor(Contour))) +
    #   geom_point()+
    #   geom_line(data = newdata2, aes(group = as.factor(Contour)), color = "black")+
    #   xlim(c(0,200))+
    #   labs(title = paste(pop_projects[p], populations[s]))
    # plot(pl)
     seasonList[[p]] <- newdata2
  }
  seas.df <- do.call("rbind", seasonList)
  seas.df$UD <- "winter"
  seas.df$population <- populations[s]
  allList[[s]] <- seas.df
}
allNewData <- do.call("rbind", allList)

save(allNewData, file = "simPredictedData.rda")
# 
pred.df <- pred.df[with(pred.df, order(projName, UD, Contour)),]
save(pred.df, file = "predictProjSummaries.rda")

# load(file = "predictProjSummaries.rda")
# pred.df$population <- factor(pred.df$population, levels = c("pop100", "pop250", "pop500", "pop750", "pop1000"))
# require(dplyr)
# 
# pred.df <- pred.df %>% dplyr::arrange(projName, Contour, population) %>%
#   dplyr::select(-c(UD, nAnimalsActual))
# write.csv(pred.df, "predictProjSummaries.csv", row.names = F)
