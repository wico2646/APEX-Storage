###Placeholder for cleaned up version of code for the beta early thaw project

## bring in the data
BetaVegCombo <- read_csv("Data/BetaVegCombo.csv") #vegetation data
BetaEnvNew <- read_csv("Data/BetaEnvNew.csv") #environmentals and flux

## grab relevant packages
library(vegan) #for veg community analysis
library(MASS) #for several things
library(rcompanion) #for transformations and misc utility functions
library(lme4) #for mixed effects models
library(nlme) #for mixed effects models
library(ggplot2) #for graphics
library(dplyr) #for organization
library(pvclust) #for cluster analysis
library(indicspecies) #for indicator species analysis

## mixed effects models
#active layer thickness
ALTmod <- lme(MaxThaw ~ Thaw.Stage * Year, random = ~1|PlotName, data = BetaEnvNew)
ALTmod
anova.lme(ALTmod)
summary(ALTmod)

#methane flux
BetaEnvNew$JulyPeakMeth.umol.m2.min <- BetaEnvNew$`JulyPeakMeth.umol/m2/min` #rename
CH4mod <- lme(JulyPeakMeth.umol.m2.min ~ Thaw.Stage * Year, random = ~1|PlotName, data = BetaEnvNew)
CH4mod
anova.lme(CH4mod)
summary(CH4mod)

#NEE
NEEmod <- lme(JulyNEE ~ Thaw.Stage * Year, random = ~1|PlotName, data = BetaEnvNew)
NEEmod
anova.lme(NEEmod)
summary(NEEmod)

#ER
ERmod <- lme(JulyER ~ Thaw.Stage * Year, random = ~1|PlotName, data = BetaEnvNew)
ERmod
anova.lme(ERmod)
summary(ERmod)

#GPP
GPPmod <- lme(GPP ~ Thaw.Stage * Year, random = ~1|PlotName, data = BetaEnvNew)
GPPmod
anova.lme(GPPmod)
summary(GPPmod)


## veg analysis
#anosim
anosimThaw<-anosim(BetaVegCombo, BetaEnvNew$Thaw.Stage, distance="bray")
summary(anosimThaw)
plot(anosimThaw)