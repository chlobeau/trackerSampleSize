library(ggplot2)
library(nlme)
library(aomisc)

load(file = "D:/Thesis/cleanArchive/RarefactionLONGwProp_20220831.rda")

# **************
# **************
# MIGRATION ----
# **************
# **************


# *************
# Elk Migs ----
# *************

elkMigsub <- rare_long[rare_long$spp2 == "Elk" & rare_long$UD == "mig",]
elkMigmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                 data = elkMigsub)

#group the data by simulation (groupedData nlme)
elkMigGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                         data=elkMigsub[,c("nAnimals", "areaProp", 
                                           "simulation", "Contour", "popEsts", "projArea_sqkm",
                                           "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
elkMigfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                data = elkMigGrp)

# look at variation between groups
plot(intervals(elkMigfL), layout = c(3,1))

## ***************************************
## elkMig1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
elkMig0 <- nlme(elkMigfL,
                 random = init + m ~ 1,
                 start = coef(elkMigmod),
                 control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(elkMig0) # maybe some high correlation between m and init
plot(elkMig0)

elkMig1 <- update(elkMig0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                  start = c(fixef(elkMig0)[1],0,fixef(elkMig0)[2],0,fixef(elkMig0)[3]))
plot(elkMig1)
summary(elkMig1)

## *************************
## elkMig2 - popEsts ----
## *************************

elkMig2 <- update(elkMig1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                  start = c(fixef(elkMig1)[1:2],0,fixef(elkMig1)[3:4],0,fixef(elkMig1)[5]))
plot(elkMig2)
summary(elkMig2)

## **************************
## elkMig3 - allDuration ----
## **************************
elkMig3 <- update(elkMig1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                  start = c(fixef(elkMig1)[1:2],0,fixef(elkMig1)[3:4],0,fixef(elkMig1)[5]))
plot(elkMig3)
summary(elkMig3) #iteration 1 didn't converge

## **************************
## elkMig4 - projWR_prop ----
## **************************
elkMig4 <- update(elkMig1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                  start = c(fixef(elkMig1)[1:2],0,fixef(elkMig1)[3:4],0,fixef(elkMig1)[5]))
plot(elkMig4)
summary(elkMig4)

## ***************************
## elkMig5 - Canopy Cover ----
## ***************************
elkMig5 <- update(elkMig1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                  start = c(fixef(elkMig1)[1:2],0,fixef(elkMig1)[3:4],0,fixef(elkMig1)[5]))
plot(elkMig5)
summary(elkMig5)

## ********************
## elkMigModSel table ----
## ********************
AIC(elkMig1, elkMig2, elkMig3, elkMig4, elkMig5)

elkMigModList <- list(elkMig1, elkMig2, elkMig3, elkMig4, elkMig5)

elkMigModSel <- data.frame(modelName = c("elkMig1", "elkMig2", "elkMig3", "elkMig4", "elkMig5"),
                         df = NA,
                         AIC = NA,
                         logLik = NA,
                         initInt = NA,
                         mInt = NA,
                         plateauInt = NA)

for(l in 1:length(elkMigModList)){
  elkMigModSel$AIC[l] <- AIC(elkMigModList[[l]])
  elkMigModSel$logLik[l] <- logLik(elkMigModList[[l]])
  elkMigModSel$initInt[l] <- fixef(elkMigModList[[l]])[names(fixef(elkMigModList[[l]])) == "init.(Intercept)"]
  elkMigModSel$mInt[l] <- fixef(elkMigModList[[l]])[names(fixef(elkMigModList[[l]])) == "m.(Intercept)"]
  
  elkMigModSel$plateauInt[l] <- fixef(elkMigModList[[l]])[names(fixef(elkMigModList[[l]])) == "plateau"]
  elkMigModSel$df[l] <- length(fixef(elkMigModList[[l]])) + 4 # I just figured this out based on looking at anova table
}


# **************
# Deer Migs ----
# **************

deerMigsub <- rare_long[rare_long$spp2 == "MuleDeer" & rare_long$UD == "mig",]
deerMigmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                 data = deerMigsub)

#group the data by simulation (groupedData nlme)
deerMigGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                         data=deerMigsub[,c("nAnimals", "areaProp", 
                                           "simulation", "Contour", "popEsts", "projArea_sqkm",
                                           "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
deerMigfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                    data = deerMigGrp)

# look at variation between groups
plot(intervals(deerMigfL), layout = c(3,1))

## ***************************************
## deerMig1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
deerMig0 <- nlme(deerMigfL,
                random = init + m ~ 1,
                start = coef(deerMigmod),
                control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(deerMig0) # maybe some high correlation between m and init
plot(deerMig0)

deerMig1 <- update(deerMig0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                  start = c(fixef(deerMig0)[1],0,fixef(deerMig0)[2],0,fixef(deerMig0)[3]))
plot(deerMig1)
summary(deerMig1)

## *************************
## deerMig2 - popEsts ----
## *************************

deerMig2 <- update(deerMig1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                  start = c(fixef(deerMig1)[1:2],0,fixef(deerMig1)[3:4],0,fixef(deerMig1)[5]))
plot(deerMig2)
summary(deerMig2)

## **************************
## deerMig3 - allDuration ----
## **************************
deerMig3 <- update(deerMig1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                  start = c(fixef(deerMig1)[1:2],0,fixef(deerMig1)[3:4],0,fixef(deerMig1)[5]))
plot(deerMig3)
summary(deerMig3) #iteration 1 didn't converge

## **************************
## deerMig4 - projWR_prop ----
## **************************
deerMig4 <- update(deerMig1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                  start = c(fixef(deerMig1)[1:2],0,fixef(deerMig1)[3:4],0,fixef(deerMig1)[5]))
plot(deerMig4)
summary(deerMig4)

## ***************************
## deerMig5 - Canopy Cover ----
## ***************************
deerMig5 <- update(deerMig1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                  start = c(fixef(deerMig1)[1:2],0,fixef(deerMig1)[3:4],0,fixef(deerMig1)[5]))
plot(deerMig5)
summary(deerMig5)

## ********************
## deerMigModSel table ----
## ********************
AIC(deerMig1, deerMig2, deerMig3, deerMig4, deerMig5)
AICc(deerMig1, deerMig2, deerMig3, deerMig4, deerMig5)

deerMigModList <- list(deerMig1, deerMig2, deerMig3, deerMig4, deerMig5)

deerMigModSel <- data.frame(modelName = c("deerMig1", "deerMig2", "deerMig3", "deerMig4", "deerMig5"),
                           df = NA,
                           AIC = NA,
                           logLik = NA,
                           initInt = NA,
                           mInt = NA,
                           plateauInt = NA)

for(l in 1:length(deerMigModList)){
  deerMigModSel$AIC[l] <- AIC(deerMigModList[[l]])
  deerMigModSel$logLik[l] <- logLik(deerMigModList[[l]])
  deerMigModSel$initInt[l] <- fixef(deerMigModList[[l]])[names(fixef(deerMigModList[[l]])) == "init.(Intercept)"]
  deerMigModSel$mInt[l] <- fixef(deerMigModList[[l]])[names(fixef(deerMigModList[[l]])) == "m.(Intercept)"]
  
  deerMigModSel$plateauInt[l] <- fixef(deerMigModList[[l]])[names(fixef(deerMigModList[[l]])) == "plateau"]
  deerMigModSel$df[l] <- length(fixef(deerMigModList[[l]])) + 4 # I just figured this out based on looking at anova table
}



# *****************
# Bighorn Migs ----
# *****************

bhMigsub <- rare_long[rare_long$spp2 == "Bighorn" & rare_long$UD == "mig",]
bhMigmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                  data = bhMigsub)

#group the data by simulation (groupedData nlme)
bhMigGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                          data=bhMigsub[,c("nAnimals", "areaProp", 
                                             "simulation", "Contour", "popEsts", "projArea_sqkm",
                                             "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
bhMigfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                     data = bhMigGrp)

# look at variation between groups
plot(intervals(bhMigfL), layout = c(3,1))

## ***************************************
## bhMig1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
bhMig0 <- nlme(bhMigfL,
                 random = init + m ~ 1,
                 start = coef(bhMigmod),
                 control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(bhMig0) # maybe some high correlation between m and init
plot(bhMig0)

bhMig1 <- update(bhMig0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                   start = c(fixef(bhMig0)[1],0,fixef(bhMig0)[2],0,fixef(bhMig0)[3]))
plot(bhMig1)
summary(bhMig1)

## *************************
## bhMig2 - popEsts ----
## *************************

bhMig2 <- update(bhMig1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                   start = c(fixef(bhMig1)[1:2],0,fixef(bhMig1)[3:4],0,fixef(bhMig1)[5]))
plot(bhMig2)
summary(bhMig2)

## **************************
## bhMig3 - allDuration ----
## **************************
bhMig3 <- update(bhMig1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                   start = c(fixef(bhMig1)[1:2],0,fixef(bhMig1)[3:4],0,fixef(bhMig1)[5]))
plot(bhMig3)
summary(bhMig3) #iteration 1 didn't converge

## **************************
## bhMig4 - projWR_prop ----
## **************************
bhMig4 <- update(bhMig1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                   start = c(fixef(bhMig1)[1:2],0,fixef(bhMig1)[3:4],0,fixef(bhMig1)[5]))
plot(bhMig4)
summary(bhMig4)

## ***************************
## bhMig5 - Canopy Cover ----
## ***************************
bhMig5 <- update(bhMig1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                   start = c(fixef(bhMig1)[1:2],0,fixef(bhMig1)[3:4],0,fixef(bhMig1)[5]))
plot(bhMig5)
summary(bhMig5)

## ********************
## bhMigModSel table ----
## ********************
AIC(bhMig1, bhMig2, bhMig3, bhMig4, bhMig5)

bhMigModList <- list(bhMig1, bhMig2, bhMig3, bhMig4, bhMig5)

bhMigModSel <- data.frame(modelName = c("bhMig1", "bhMig2", "bhMig3", "bhMig4", "bhMig5"),
                            df = NA,
                            AIC = NA,
                            logLik = NA,
                            initInt = NA,
                            mInt = NA,
                            plateauInt = NA)

for(l in 1:length(bhMigModList)){
  bhMigModSel$AIC[l] <- AIC(bhMigModList[[l]])
  bhMigModSel$logLik[l] <- logLik(bhMigModList[[l]])
  bhMigModSel$initInt[l] <- fixef(bhMigModList[[l]])[names(fixef(bhMigModList[[l]])) == "init.(Intercept)"]
  bhMigModSel$mInt[l] <- fixef(bhMigModList[[l]])[names(fixef(bhMigModList[[l]])) == "m.(Intercept)"]
  
  bhMigModSel$plateauInt[l] <- fixef(bhMigModList[[l]])[names(fixef(bhMigModList[[l]])) == "plateau"]
  bhMigModSel$df[l] <- length(fixef(bhMigModList[[l]])) + 4 # I just figured this out based on looking at anova table
}



# *****************
# *****************
# WINTER RANGE ----
# *****************
# *****************


# *************
# Elk WR ----
# *************

elkWRsub <- rare_long[rare_long$spp2 == "Elk" & rare_long$UD == "winter",]
elkWRmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                 data = elkWRsub)

#group the data by simulation (groupedData nlme)
elkWRGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                         data=elkWRsub[,c("nAnimals", "areaProp", 
                                           "simulation", "Contour", "popEsts", "projArea_sqkm",
                                           "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
elkWRfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                    data = elkWRGrp)

# look at variation between groups
plot(intervals(elkWRfL), layout = c(3,1))

## ***************************************
## elkWR1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
elkWR0 <- nlme(elkWRfL,
                random = init + m ~ 1,
                start = coef(elkWRmod),
                control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(elkWR0) # maybe some high correlation between m and init
plot(elkWR0)

elkWR1 <- update(elkWR0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                  start = c(fixef(elkWR0)[1],0,fixef(elkWR0)[2],0,fixef(elkWR0)[3]))
plot(elkWR1)
summary(elkWR1)

## *************************
## elkWR2 - popEsts ----
## *************************

elkWR2 <- update(elkWR1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                  start = c(fixef(elkWR1)[1:2],0,fixef(elkWR1)[3:4],0,fixef(elkWR1)[5]))
plot(elkWR2)
summary(elkWR2)

## **************************
## elkWR3 - allDuration ----
## **************************
elkWR3 <- update(elkWR1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                  start = c(fixef(elkWR1)[1:2],0,fixef(elkWR1)[3:4],0,fixef(elkWR1)[5]))
plot(elkWR3)
summary(elkWR3) #iteration 1 didn't converge

## **************************
## elkWR4 - projWR_prop ----
## **************************
elkWR4 <- update(elkWR1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                  start = c(fixef(elkWR1)[1:2],0,fixef(elkWR1)[3:4],0,fixef(elkWR1)[5]))
plot(elkWR4)
summary(elkWR4)

## ***************************
## elkWR5 - Canopy Cover ----
## ***************************
elkWR5 <- update(elkWR1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                  start = c(fixef(elkWR1)[1:2],0,fixef(elkWR1)[3:4],0,fixef(elkWR1)[5]))
plot(elkWR5)
summary(elkWR5)

## ********************
## elkWRModSel table ----
## ********************
AIC(elkWR1, elkWR2, elkWR3, elkWR4, elkWR5)

elkWRModList <- list(elkWR1, elkWR2, elkWR3, elkWR4, elkWR5)

elkWRModSel <- data.frame(modelName = c("elkWR1", "elkWR2", "elkWR3", "elkWR4", "elkWR5"),
                           df = NA,
                           AIC = NA,
                           logLik = NA,
                           initInt = NA,
                           mInt = NA,
                           plateauInt = NA)

for(l in 1:length(elkWRModList)){
  elkWRModSel$AIC[l] <- AIC(elkWRModList[[l]])
  elkWRModSel$logLik[l] <- logLik(elkWRModList[[l]])
  elkWRModSel$initInt[l] <- fixef(elkWRModList[[l]])[names(fixef(elkWRModList[[l]])) == "init.(Intercept)"]
  elkWRModSel$mInt[l] <- fixef(elkWRModList[[l]])[names(fixef(elkWRModList[[l]])) == "m.(Intercept)"]
  
  elkWRModSel$plateauInt[l] <- fixef(elkWRModList[[l]])[names(fixef(elkWRModList[[l]])) == "plateau"]
  elkWRModSel$df[l] <- length(fixef(elkWRModList[[l]])) + 4 # I just figured this out based on looking at anova table
}

## ******************* 
## predict elkWR ----
## *******************

summary(elkWRGrp$popEsts)

elkWRNdat <- expand.grid(nAnimals = 1:200, Contour = c(50,75,90,95,99) 
                          , popEsts = c(min(elkWRGrp$popEsts), mean(elkWRGrp$popEsts), quantile(elkWRGrp$popEsts, .75))
                          #,projWR_prop = c(min(elkWRGrp$projWR_prop), median(elkWRGrp$projWR_prop), max(elkWRGrp$projWR_prop))
                          #, allDuration = c(1,2,3,5)
)
elkWRNdat$simulation <- 1
elkWRNdat$areaProp <- predict(elkWR2, newdata = elkWRNdat)
elkWRNdat$spp2 <- "Elk"
elkWRNdat$UD <- "winter"

elkWRNdat$popEstsNumber[elkWRNdat$popEsts == min(elkWRGrp$popEsts)] <- 1
elkWRNdat$popEstsNumber[elkWRNdat$popEsts == mean(elkWRGrp$popEsts)] <- 2
elkWRNdat$popEstsNumber[elkWRNdat$popEsts == quantile(elkWRGrp$popEsts, .75)] <- 3

ggplot(elkWRNdat, aes(nAnimals, areaProp, color = as.factor(popEsts))) + 
  geom_line() +
  facet_grid(~ Contour) +
  ylim(0,1)+
  xlim(0,100)+
  labs(title = "Elk Winter range")




# **************
# Deer WR ----
# **************

deerWRsub <- rare_long[rare_long$spp2 == "MuleDeer" & rare_long$UD == "winter",]
deerWRmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                  data = deerWRsub)

#group the data by simulation (groupedData nlme)
deerWRGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                          data=deerWRsub[,c("nAnimals", "areaProp", 
                                             "simulation", "Contour", "popEsts", "projArea_sqkm",
                                             "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
deerWRfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                     data = deerWRGrp)

# look at variation between groups
plot(intervals(deerWRfL), layout = c(3,1))

## ***************************************
## deerWR1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
deerWR0 <- nlme(deerWRfL,
                 random = init + m ~ 1,
                 start = coef(deerWRmod),
                 control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(deerWR0) # maybe some high correlation between m and init
plot(deerWR0)

deerWR1 <- update(deerWR0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                   start = c(fixef(deerWR0)[1],0,fixef(deerWR0)[2],0,fixef(deerWR0)[3]))
plot(deerWR1)
summary(deerWR1)

## *************************
## deerWR2 - popEsts ----
## *************************

deerWR2 <- update(deerWR1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                   start = c(fixef(deerWR1)[1:2],0,fixef(deerWR1)[3:4],0,fixef(deerWR1)[5]))
plot(deerWR2)
summary(deerWR2)

## **************************
## deerWR3 - allDuration ----
## **************************
deerWR3 <- update(deerWR1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                   start = c(fixef(deerWR1)[1:2],0,fixef(deerWR1)[3:4],0,fixef(deerWR1)[5]))
plot(deerWR3)
summary(deerWR3) #iteration 1 didn't converge

## **************************
## deerWR4 - projWR_prop ----
## **************************
deerWR4 <- update(deerWR1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                   start = c(fixef(deerWR1)[1:2],0,fixef(deerWR1)[3:4],0,fixef(deerWR1)[5]))
plot(deerWR4)
summary(deerWR4)

## ***************************
## deerWR5 - Canopy Cover ----
## ***************************
deerWR5 <- update(deerWR1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                   start = c(fixef(deerWR1)[1:2],0,fixef(deerWR1)[3:4],0,fixef(deerWR1)[5]))
plot(deerWR5)
summary(deerWR5)

## ********************
## deerWRModSel table ----
## ********************
AIC(deerWR1, deerWR2, deerWR3, deerWR4, deerWR5)

deerWRModList <- list(deerWR1, deerWR2, deerWR3, deerWR4, deerWR5)

deerWRModSel <- data.frame(modelName = c("deerWR1", "deerWR2", "deerWR3", "deerWR4", "deerWR5"),
                            df = NA,
                            AIC = NA,
                            logLik = NA,
                            initInt = NA,
                            mInt = NA,
                            plateauInt = NA)

for(l in 1:length(deerWRModList)){
  deerWRModSel$AIC[l] <- AIC(deerWRModList[[l]])
  deerWRModSel$logLik[l] <- logLik(deerWRModList[[l]])
  deerWRModSel$initInt[l] <- fixef(deerWRModList[[l]])[names(fixef(deerWRModList[[l]])) == "init.(Intercept)"]
  deerWRModSel$mInt[l] <- fixef(deerWRModList[[l]])[names(fixef(deerWRModList[[l]])) == "m.(Intercept)"]
  
  deerWRModSel$plateauInt[l] <- fixef(deerWRModList[[l]])[names(fixef(deerWRModList[[l]])) == "plateau"]
  deerWRModSel$df[l] <- length(fixef(deerWRModList[[l]])) + 4 # I just figured this out based on looking at anova table
}

## ******************* 
## predict deerWR ----
## *******************

summary(deerWRGrp$popEsts)

deerWRNdat <- expand.grid(nAnimals = 1:200, Contour = c(50,75,90,95,99) 
                           , popEsts = c(min(deerWRGrp$popEsts), mean(deerWRGrp$popEsts), quantile(deerWRGrp$popEsts, .75))
                           #,projWR_prop = c(min(deerWRGrp$projWR_prop), median(deerWRGrp$projWR_prop), max(deerWRGrp$projWR_prop))
                           #, allDuration = c(1,2,3,5)
)
deerWRNdat$simulation <- 1
deerWRNdat$areaProp <- predict(deerWR2, newdata = deerWRNdat)
deerWRNdat$spp2 <- "MuleDeer"
deerWRNdat$UD <- "winter"

deerWRNdat$popEstsNumber[deerWRNdat$popEsts == min(deerWRGrp$popEsts)] <- 1
deerWRNdat$popEstsNumber[deerWRNdat$popEsts == mean(deerWRGrp$popEsts)] <- 2
deerWRNdat$popEstsNumber[deerWRNdat$popEsts == quantile(deerWRGrp$popEsts, .75)] <- 3

ggplot(deerWRNdat, aes(nAnimals, areaProp, color = as.factor(popEsts))) + 
  geom_line() +
  facet_grid(~ Contour) +
  ylim(0,1)+
  xlim(0,100)+
  labs(title = "Deer Winter range")


# *****************
# Bighorn WR ----
# *****************

bhWRsub <- rare_long[rare_long$spp2 == "Bighorn" & rare_long$UD == "winter",]
bhWRmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                data = bhWRsub)

#group the data by simulation (groupedData nlme)
bhWRGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                        data=bhWRsub[,c("nAnimals", "areaProp", 
                                         "simulation", "Contour", "popEsts", "projArea_sqkm",
                                         "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
bhWRfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                   data = bhWRGrp)

# look at variation between groups
plot(intervals(bhWRfL), layout = c(3,1))

## ***************************************
## bhWR1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
bhWR0 <- nlme(bhWRfL,
               random = init + m ~ 1,
               start = coef(bhWRmod),
               control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(bhWR0) # maybe some high correlation between m and init
plot(bhWR0)

bhWR1 <- update(bhWR0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                 start = c(fixef(bhWR0)[1],0,fixef(bhWR0)[2],0,fixef(bhWR0)[3]))
plot(bhWR1)
summary(bhWR1)

## *************************
## bhWR2 - popEsts ----
## *************************

bhWR2 <- update(bhWR1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                 start = c(fixef(bhWR1)[1:2],0,fixef(bhWR1)[3:4],0,fixef(bhWR1)[5]))
plot(bhWR2)
summary(bhWR2)

## **************************
## bhWR3 - allDuration ----
## **************************
bhWR3 <- update(bhWR1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                 start = c(fixef(bhWR1)[1:2],0,fixef(bhWR1)[3:4],0,fixef(bhWR1)[5]))
plot(bhWR3)
summary(bhWR3) #iteration 1 didn't converge

## **************************
## bhWR4 - projWR_prop ----
## **************************
bhWR4 <- update(bhWR1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                 start = c(fixef(bhWR1)[1:2],0,fixef(bhWR1)[3:4],0,fixef(bhWR1)[5]))
plot(bhWR4)
summary(bhWR4)

## ***************************
## bhWR5 - Canopy Cover ----
## ***************************
bhWR5 <- update(bhWR1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                 start = c(fixef(bhWR1)[1:2],0,fixef(bhWR1)[3:4],0,fixef(bhWR1)[5]))
plot(bhWR5)
summary(bhWR5)

## ********************
## bhWRModSel table ----
## ********************
AIC(bhWR1, bhWR2, bhWR3, bhWR4, bhWR5)

bhWRModList <- list(bhWR1, bhWR2, bhWR3, bhWR4, bhWR5)

bhWRModSel <- data.frame(modelName = c("bhWR1", "bhWR2", "bhWR3", "bhWR4", "bhWR5"),
                          df = NA,
                          AIC = NA,
                          logLik = NA,
                          initInt = NA,
                          mInt = NA,
                          plateauInt = NA)

for(l in 1:length(bhWRModList)){
  bhWRModSel$AIC[l] <- AIC(bhWRModList[[l]])
  bhWRModSel$logLik[l] <- logLik(bhWRModList[[l]])
  bhWRModSel$initInt[l] <- fixef(bhWRModList[[l]])[names(fixef(bhWRModList[[l]])) == "init.(Intercept)"]
  bhWRModSel$mInt[l] <- fixef(bhWRModList[[l]])[names(fixef(bhWRModList[[l]])) == "m.(Intercept)"]
  
  bhWRModSel$plateauInt[l] <- fixef(bhWRModList[[l]])[names(fixef(bhWRModList[[l]])) == "plateau"]
  bhWRModSel$df[l] <- length(fixef(bhWRModList[[l]])) + 4 # I just figured this out based on looking at anova table
}

## ******************* 
## predict bhWR ----
## *******************

summary(bhWRGrp$popEsts)

bhWRNdat <- expand.grid(nAnimals = 1:200, Contour = c(50,75,90,95,99) 
                         , popEsts = c(min(bhWRGrp$popEsts), median(bhWRGrp$popEsts), quantile(bhWRGrp$popEsts, .75))
                         #,projWR_prop = c(min(bhWRGrp$projWR_prop), median(bhWRGrp$projWR_prop), max(bhWRGrp$projWR_prop))
                         #, allDuration = c(1,2,3,5)
)
bhWRNdat$simulation <- 1
bhWRNdat$areaProp <- predict(bhWR2, newdata = bhWRNdat)
bhWRNdat$spp2 <- "Bighorn"
bhWRNdat$UD <- "winter"

bhWRNdat$popEstsNumber[bhWRNdat$popEsts == min(bhWRGrp$popEsts)] <- 1
bhWRNdat$popEstsNumber[bhWRNdat$popEsts == median(bhWRGrp$popEsts)] <- 2
bhWRNdat$popEstsNumber[bhWRNdat$popEsts == quantile(bhWRGrp$popEsts, .75)] <- 3

ggplot(bhWRNdat, aes(nAnimals, areaProp, color = as.factor(popEsts))) + 
  geom_line() +
  facet_grid(~ Contour) +
  ylim(0,1)+
  xlim(0,20)+
  labs(title = "Bighorn Winter range")



# *****************
# *****************
# SUMMER RANGE ----
# *****************
# *****************

# *************
# Elk SR ----
# *************

elkSRsub <- rare_long[rare_long$spp2 == "Elk" & rare_long$UD == "summer",]
elkSRmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                 data = elkSRsub)

#group the data by simulation (groupedData nlme)
elkSRGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                         data=elkSRsub[,c("nAnimals", "areaProp", 
                                           "simulation", "Contour", "popEsts", "projArea_sqkm",
                                           "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
elkSRfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                    data = elkSRGrp)

# look at variation between groups
plot(intervals(elkSRfL), layout = c(3,1))

## ***************************************
## elkSR1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
elkSR0 <- nlme(elkSRfL,
                random = init + m ~ 1,
                start = coef(elkSRmod),
                control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(elkSR0) # maybe some high correlation between m and init
plot(elkSR0)

elkSR1 <- update(elkSR0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                  start = c(fixef(elkSR0)[1],0,fixef(elkSR0)[2],0,fixef(elkSR0)[3]))
plot(elkSR1)
summary(elkSR1)

## *************************
## elkSR2 - popEsts ----
## *************************

elkSR2 <- update(elkSR1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                  start = c(fixef(elkSR1)[1:2],0,fixef(elkSR1)[3:4],0,fixef(elkSR1)[5]))
plot(elkSR2)
summary(elkSR2)

## **************************
## elkSR3 - allDuration ----
## **************************
elkSR3 <- update(elkSR1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                  start = c(fixef(elkSR1)[1:2],0,fixef(elkSR1)[3:4],0,fixef(elkSR1)[5]))
plot(elkSR3)
summary(elkSR3) #iteration 1 didn't converge

## **************************
## elkSR4 - projWR_prop ----
## **************************
elkSR4 <- update(elkSR1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                  start = c(fixef(elkSR1)[1:2],0,fixef(elkSR1)[3:4],0,fixef(elkSR1)[5]))
plot(elkSR4)
summary(elkSR4)

## ***************************
## elkSR5 - Canopy Cover ----
## ***************************
elkSR5 <- update(elkSR1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                  start = c(fixef(elkSR1)[1:2],0,fixef(elkSR1)[3:4],0,fixef(elkSR1)[5]))
plot(elkSR5)
summary(elkSR5)

## ********************
## elkSRModSel table ----
## ********************
AIC(elkSR1, elkSR2, elkSR3, elkSR4, elkSR5)

elkSRModList <- list(elkSR1, elkSR2, elkSR3, elkSR4, elkSR5)

elkSRModSel <- data.frame(modelName = c("elkSR1", "elkSR2", "elkSR3", "elkSR4", "elkSR5"),
                           df = NA,
                           AIC = NA,
                           logLik = NA,
                           initInt = NA,
                           mInt = NA,
                           plateauInt = NA)

for(l in 1:length(elkSRModList)){
  elkSRModSel$AIC[l] <- AIC(elkSRModList[[l]])
  elkSRModSel$logLik[l] <- logLik(elkSRModList[[l]])
  elkSRModSel$initInt[l] <- fixef(elkSRModList[[l]])[names(fixef(elkSRModList[[l]])) == "init.(Intercept)"]
  elkSRModSel$mInt[l] <- fixef(elkSRModList[[l]])[names(fixef(elkSRModList[[l]])) == "m.(Intercept)"]
  
  elkSRModSel$plateauInt[l] <- fixef(elkSRModList[[l]])[names(fixef(elkSRModList[[l]])) == "plateau"]
  elkSRModSel$df[l] <- length(fixef(elkSRModList[[l]])) + 4 # I just figured this out based on looking at anova table
}

## ******************* 
## predict elkSR ----
## *******************

summary(elkSRGrp$popEsts)

elkSRNdat <- expand.grid(nAnimals = 1:200, Contour = c(50,75,90,95,99) 
                          , popEsts = c(min(elkSRGrp$popEsts), mean(elkSRGrp$popEsts), quantile(elkSRGrp$popEsts, .75))
                          #,projWR_prop = c(min(elkSRGrp$projWR_prop), median(elkSRGrp$projWR_prop), max(elkSRGrp$projWR_prop))
                          #, allDuration = c(1,2,3,5)
)
elkSRNdat$simulation <- 1
elkSRNdat$areaProp <- predict(elkSR2, newdata = elkSRNdat)
elkSRNdat$spp2 <- "Elk"
elkSRNdat$UD <- "summer"

elkSRNdat$popEstsNumber[elkSRNdat$popEsts == min(elkSRGrp$popEsts)] <- 1
elkSRNdat$popEstsNumber[elkSRNdat$popEsts == mean(elkSRGrp$popEsts)] <- 2
elkSRNdat$popEstsNumber[elkSRNdat$popEsts == quantile(elkSRGrp$popEsts, .75)] <- 3

ggplot(elkSRNdat, aes(nAnimals, areaProp, color = as.factor(popEsts))) + 
  geom_line() +
  facet_grid(~ Contour) +
  ylim(0,1)+
  xlim(0,100)+
  labs(title = "Elk summer")

# **************
# Deer SR ----
# **************

deerSRsub <- rare_long[rare_long$spp2 == "MuleDeer" & rare_long$UD == "summer",]
deerSRmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                  data = deerSRsub)

#group the data by simulation (groupedData nlme)
deerSRGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                          data=deerSRsub[,c("nAnimals", "areaProp", 
                                             "simulation", "Contour", "popEsts", "projArea_sqkm",
                                             "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
deerSRfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                     data = deerSRGrp)

# look at variation between groups
plot(intervals(deerSRfL), layout = c(3,1))

## ***************************************
## deerSR1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
deerSR0 <- nlme(deerSRfL,
                 random = init + m ~ 1,
                 start = coef(deerSRmod),
                 control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(deerSR0) # maybe some high correlation between m and init
plot(deerSR0)

deerSR1 <- update(deerSR0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                   start = c(fixef(deerSR0)[1],0,fixef(deerSR0)[2],0,fixef(deerSR0)[3]))
plot(deerSR1)
summary(deerSR1)

## *************************
## deerSR2 - popEsts ----
## *************************

deerSR2 <- update(deerSR1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                   start = c(fixef(deerSR1)[1:2],0,fixef(deerSR1)[3:4],0,fixef(deerSR1)[5]))
plot(deerSR2)
summary(deerSR2)

## **************************
## deerSR3 - allDuration ----
## **************************
deerSR3 <- update(deerSR1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                   start = c(fixef(deerSR1)[1:2],0,fixef(deerSR1)[3:4],0,fixef(deerSR1)[5]))
plot(deerSR3)
summary(deerSR3) #iteration 1 didn't converge

## **************************
## deerSR4 - projWR_prop ----
## **************************
deerSR4 <- update(deerSR1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                   start = c(fixef(deerSR2)[1:2],0,fixef(deerSR2)[4:5],0,fixef(deerSR2)[7]))
plot(deerSR4)
summary(deerSR4) #didn't run

## ***************************
## deerSR5 - Canopy Cover ----
## ***************************
deerSR5 <- update(deerSR1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                   start = c(fixef(deerSR2)[1:2],0,fixef(deerSR2)[4:5],0,fixef(deerSR2)[7]))
plot(deerSR5)
summary(deerSR5) #didn't run

## ********************
## deerSRModSel table ----
## ********************
AIC(deerSR1, deerSR2, deerSR3)
summary(deerSR2)

deerSRModList <- list(deerSR1, deerSR2, deerSR3)

deerSRModSel <- data.frame(modelName = c("deerSR1", "deerSR2", "deerSR3"),
                            df = NA,
                            AIC = NA,
                            logLik = NA,
                            initInt = NA,
                            mInt = NA,
                            plateauInt = NA)

for(l in 1:length(deerSRModList)){
  deerSRModSel$AIC[l] <- AIC(deerSRModList[[l]])
  deerSRModSel$logLik[l] <- logLik(deerSRModList[[l]])
  deerSRModSel$initInt[l] <- fixef(deerSRModList[[l]])[names(fixef(deerSRModList[[l]])) == "init.(Intercept)"]
  deerSRModSel$mInt[l] <- fixef(deerSRModList[[l]])[names(fixef(deerSRModList[[l]])) == "m.(Intercept)"]
  
  deerSRModSel$plateauInt[l] <- fixef(deerSRModList[[l]])[names(fixef(deerSRModList[[l]])) == "plateau"]
  deerSRModSel$df[l] <- length(fixef(deerSRModList[[l]])) + 4 # I just figured this out based on looking at anova table
}

## ******************* 
## predict deerSR ----
## *******************

summary(deerSRGrp$popEsts)

deerSRNdat <- expand.grid(nAnimals = 1:200, Contour = c(50,75,90,95,99) 
                           , popEsts = c(min(deerSRGrp$popEsts), mean(deerSRGrp$popEsts), quantile(deerSRGrp$popEsts, .75))
                           #,projWR_prop = c(min(deerSRGrp$projWR_prop), median(deerSRGrp$projWR_prop), max(deerSRGrp$projWR_prop))
                           #, allDuration = c(1,2,3,5)
)
deerSRNdat$simulation <- 1
deerSRNdat$areaProp <- predict(deerSR2, newdata = deerSRNdat)
deerSRNdat$spp2 <- "MuleDeer"
deerSRNdat$UD <- "summer"

deerSRNdat$popEstsNumber[deerSRNdat$popEsts == min(deerSRGrp$popEsts)] <- 1
deerSRNdat$popEstsNumber[deerSRNdat$popEsts == mean(deerSRGrp$popEsts)] <- 2
deerSRNdat$popEstsNumber[deerSRNdat$popEsts == quantile(deerSRGrp$popEsts, .75)] <- 3

ggplot(deerSRNdat, aes(nAnimals, areaProp, color = as.factor(popEsts))) + 
  geom_line() +
  facet_grid(~ Contour) +
  ylim(0,1)+
  xlim(0,200)+
  labs(title = "Deer summer")


# *****************
# Bighorn SR ----
# *****************

bhSRsub <- rare_long[rare_long$spp2 == "Bighorn" & rare_long$UD == "summer",]
bhSRmod <- drm(area_sqkm ~ nAnimals, fct = DRC.asymReg(),
                data = bhSRsub)

#group the data by simulation (groupedData nlme)
bhSRGrp <- groupedData(areaProp ~ nAnimals | simulation, 
                        data=bhSRsub[,c("nAnimals", "areaProp", 
                                         "simulation", "Contour", "popEsts", "projArea_sqkm",
                                         "allDuration", "projWR_prop", "ccMean")])


#cycle through all groups
bhSRfL <- nlsList(areaProp ~ NLS.asymReg(nAnimals, init, m, plateau),
                   data = bhSRGrp)

# look at variation between groups
plot(intervals(bhSRfL), layout = c(3,1))

## ***************************************
## bhSR1 - Mixed Model with Contour ----
## ***************************************

# fit mixed model to group
# use starting coefficients from naive model
# use random intercepts for m only 
# relax convergence criteria (control)
bhSR0 <- nlme(bhSRfL,
               random = init + m ~ 1,
               start = coef(bhSRmod),
               control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))
summary(bhSR0) # maybe some high correlation between m and init
plot(bhSR0)

bhSR1 <- update(bhSR0, fixed = list(init + m  ~ Contour, plateau ~ 1),
                 start = c(fixef(bhSR0)[1],0,fixef(bhSR0)[2],0,fixef(bhSR0)[3]))
plot(bhSR1)
summary(bhSR1)

## *************************
## bhSR2 - popEsts ----
## *************************

bhSR2 <- update(bhSR1, fixed = list(init + m  ~ Contour + popEsts, plateau ~ 1),
                 start = c(fixef(bhSR1)[1:2],0,fixef(bhSR1)[3:4],0,fixef(bhSR1)[5]))
plot(bhSR2)
summary(bhSR2)

## **************************
## bhSR3 - allDuration ----
## **************************
bhSR3 <- update(bhSR1, fixed = list(init + m  ~ Contour + allDuration, plateau ~ 1),
                 start = c(fixef(bhSR1)[1:2],0,fixef(bhSR1)[3:4],0,fixef(bhSR1)[5]))
plot(bhSR3)
summary(bhSR3) #iteration 1 didn't converge

## **************************
## bhSR4 - projWR_prop ----
## **************************
bhSR4 <- update(bhSR1, fixed = list(init + m  ~ Contour + projWR_prop, plateau ~ 1),
                 start = c(fixef(bhSR1)[1:2],0,fixef(bhSR1)[3:4],0,fixef(bhSR1)[5]))
plot(bhSR4)
summary(bhSR4)

## ***************************
## bhSR5 - Canopy Cover ----
## ***************************
bhSR5 <- update(bhSR1, fixed = list(init + m  ~ Contour + ccMean, plateau ~ 1),
                 start = c(fixef(bhSR1)[1:2],0,fixef(bhSR1)[3:4],0,fixef(bhSR1)[5]))
plot(bhSR5)
summary(bhSR5)

## ********************
## bhSRModSel table ----
## ********************
AIC(bhSR1, bhSR2, bhSR3, bhSR4, bhSR5)

bhSRModList <- list(bhSR1, bhSR2, bhSR3, bhSR4, bhSR5)

bhSRModSel <- data.frame(modelName = c("bhSR1", "bhSR2", "bhSR3", "bhSR4", "bhSR5"),
                          df = NA,
                          AIC = NA,
                          logLik = NA,
                          initInt = NA,
                          mInt = NA,
                          plateauInt = NA)

for(l in 1:length(bhSRModList)){
  bhSRModSel$AIC[l] <- AIC(bhSRModList[[l]])
  bhSRModSel$logLik[l] <- logLik(bhSRModList[[l]])
  bhSRModSel$initInt[l] <- fixef(bhSRModList[[l]])[names(fixef(bhSRModList[[l]])) == "init.(Intercept)"]
  bhSRModSel$mInt[l] <- fixef(bhSRModList[[l]])[names(fixef(bhSRModList[[l]])) == "m.(Intercept)"]
  
  bhSRModSel$plateauInt[l] <- fixef(bhSRModList[[l]])[names(fixef(bhSRModList[[l]])) == "plateau"]
  bhSRModSel$df[l] <- length(fixef(bhSRModList[[l]])) + 4 # I just figured this out based on looking at anova table
}

## ******************* 
## predict bhSR ----
## *******************

summary(bhSRGrp$popEsts)

bhSRNdat <- expand.grid(nAnimals = 1:200, Contour = c(50,75,90,95,99) 
                         , popEsts = c(min(bhSRGrp$popEsts), median(bhSRGrp$popEsts), quantile(bhSRGrp$popEsts, .75))
                         #,projWR_prop = c(min(bhSRGrp$projWR_prop), median(bhSRGrp$projWR_prop), max(bhSRGrp$projWR_prop))
                         #, allDuration = c(1,2,3,5)
)
bhSRNdat$simulation <- 1
bhSRNdat$areaProp <- predict(bhSR2, newdata = bhSRNdat)
bhSRNdat$spp2 <- "Bighorn"
bhSRNdat$UD <- "summer"

bhSRNdat$popEstsNumber[bhSRNdat$popEsts == min(bhSRGrp$popEsts)] <- 1
bhSRNdat$popEstsNumber[bhSRNdat$popEsts == median(bhSRGrp$popEsts)] <- 2
bhSRNdat$popEstsNumber[bhSRNdat$popEsts == quantile(bhSRGrp$popEsts, .75)] <- 3

ggplot(bhSRNdat, aes(nAnimals, areaProp, color = as.factor(popEsts))) + 
  geom_line() +
  facet_grid(~ Contour) +
  ylim(0,1)+
  xlim(0,20)+
  labs(title = "Bighorn summer")

# *******************************
# Summarize model parameters ----
# *******************************
library(dplyr)
library(tidyr)

migParams <- data.frame()

elkMig2p <- data.frame(coef(summary(elkMig2)))
elkMig2p$sppUD <- "elkMig"
elkMig2p$param <- row.names(elkMig2p)
migParams <- rbind(elkMig2p, migParams)

deerMig2p <- data.frame(coef(summary(deerMig2)))
deerMig2p$sppUD <- "deerMig"
deerMig2p$param <- row.names(deerMig2p)
migParams <- rbind(deerMig2p, migParams)

bhMig2p <- data.frame(coef(summary(bhMig2)))
bhMig2p$sppUD <- "bhMig"
bhMig2p$param <- row.names(bhMig2p)
migParams <- rbind(bhMig2p, migParams)

migParams_wide <- migParams %>%
  select(sppUD, param, Value, Std.Error) %>%
  pivot_wider(names_from = sppUD,
              values_from = c(Value, Std.Error))

wrParams <- data.frame()
elkWR2p <- data.frame(coef(summary(elkWR2)))
elkWR2p$sppUD <- "elkWR"
elkWR2p$param <- row.names(elkWR2p)
wrParams <- rbind(elkWR2p, wrParams)

deerWR2p <- data.frame(coef(summary(deerWR2)))
deerWR2p$sppUD <- "deerWR"
deerWR2p$param <- row.names(deerWR2p)
wrParams <- rbind(deerWR2p, wrParams)

bhWR2p <- data.frame(coef(summary(bhWR2)))
bhWR2p$sppUD <- "bhWR"
bhWR2p$param <- row.names(bhWR2p)
wrParams <- rbind(bhWR2p, wrParams)


wrParams_wide <- wrParams %>%
  select(sppUD, param, Value, Std.Error) %>%
  pivot_wider(names_from = sppUD,
              values_from = c(Value, Std.Error))

elkSR2p <- data.frame(coef(summary(elkSR2)))
elkSR2p$sppUD <- "elkSR"
elkSR2p$param <- row.names(elkSR2p)
srParams <- elkSR2p

deerSR2p <- data.frame(coef(summary(deerSR2)))
deerSR2p$sppUD <- "deerSR"
deerSR2p$param <- row.names(deerSR2p)
srParams <- rbind(deerSR2p, srParams)

bhSR2p <- data.frame(coef(summary(bhSR2)))
bhSR2p$sppUD <- "bhSR"
bhSR2p$param <- row.names(bhSR2p)
srParams <- rbind(bhSR2p, srParams)

srParams_wide <- srParams %>%
  select(sppUD, param, Value, Std.Error) %>%
  pivot_wider(names_from = sppUD,
              values_from = c(Value, Std.Error))

write.csv(migParams_wide, "migParams_20220727.csv", row.names = F)
write.csv(wrParams_wide, "wrParams_20220727.csv", row.names = F)
write.csv(srParams_wide, "srParams_20220727.csv", row.names = F)

allModSel <- rbind(elkMigModSel, deerMigModSel, bhMigModSel,
                   elkWRModSel, deerWRModSel, bhWRModSel,
                   elkSRModSel, deerSRModSel, bhSRModSel)

write.csv(allModSel, "allModSel_20220727.csv", row.names = F)

#save.image(file = "F:/Thesis/DataAnalysis/ModelSel/areaProp_modelSel_20220727.RData")

## ******************* 
## predict Box 1 ----
## *******************

summary(deerMigGrp$popEsts)

deerMigNdat <- expand.grid(nAnimals = 1:200
                           #, Contour = c(50,75,90,95,99) 
                           , Contour = 50
                         #, popEsts = c(min(deerMigGrp$popEsts), mean(deerMigGrp$popEsts), quantile(deerMigGrp$popEsts, .75))
                         , popEsts = 10000
)
deerMigNdat$simulation <- 1
deerMigNdat$areaProp <- predict(deerMig2, newdata = deerMigNdat)

.99 * as.numeric(fixef(deerMig2)[names(fixef(deerMig2)) == "plateau"])

min(deerMigNdat$nAnimals[deerMigNdat$areaProp >= 0.941808])
