###Placeholder for cleaned up version of code for the beta early thaw project

## bring in the data
BetaVegCombo <- read_csv("Data/BetaVegCombo.csv") #vegetation data
#rownames(BetaVegCombo) <- BetaVegCombo[,1] #deprecated
BetaVegCombo = BetaVegCombo[,-1] #ignore the first row for analyses
BetaEnvNew <- read_csv("Data/BetaEnvNew.csv") #environmentals and flux
BetaEnvNew$Year <- as.factor(BetaEnvNew$Year) #change year to a factor
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

#simper
betaSIMPER<-simper(BetaVegCombo, BetaEnvNew$Thaw.Stage)
summary(betaSIMPER) #breakdown by species
lapply(betaSIMPER, FUN=function(x){x$overall}) #summarize percent dissimilarity by group

#ward's cluster analysis

#NMDS fig

##indicator species analyses
multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage)
indval <- multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage, control = how(nperm=999))
summary.multipatt(indval, indvalcomp = TRUE)
indval$sign

indvalori <- multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage, duleg = TRUE,             
                      control = how(nperm=999))
summary.multipatt(indvalori, indvalcomp = TRUE)

indvaloriward <- multipatt(BetaVegCombo, BetaEnvNew$WardCluster, duleg = TRUE,             
                          control = how(nperm=999))
summary(indvaloriward)

#including species pairs
BetaComb <- combinespecies(BetaVegCombo, max.order = 2)$XC
dim(BetaComb)
pairs <- multipatt(BetaComb, BetaEnvNew$Thaw.Stage, duleg = TRUE, 
                 control = how(nperm=999))
summary.multipatt(pairs, indvalcomp = TRUE)

indvalspcomb <- multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage, duleg = TRUE, 
                         control = how(nperm=999))
summary.multipatt(indvalspcomb, indvalcomp = TRUE)

#indicators function for pairs and trios
#stable
scstable <- indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=0,           
                     max.order = 3, verbose=TRUE,             
                     At=0.5, Bt=0.2, enableFixed = TRUE)
print(scstable, sqrtIVt = 0.6)
sc2stable <- pruneindicators(scstable, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2stable)

#early
scearly <- indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=1,           
                    max.order = 3, verbose=TRUE,             
                    At=0.5, Bt=0.2, enableFixed = TRUE)
print(scearly, sqrtIVt = 0.6)
sc2early <- pruneindicators(scearly, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2early)

#intermediate
scint <- indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=2,           
                  max.order = 3, verbose=TRUE,             
                  At=0.5, Bt=0.2, enableFixed = TRUE)
print(scint, sqrtIVt = 0.6)
sc2int <- pruneindicators(scint, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2int)

#advanced
scadv <- indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=3,           
                  max.order = 3, verbose=TRUE,             
                  At=0.5, Bt=0.2, enableFixed = TRUE)
print(scadv, sqrtIVt = 0.6)
sc2adv <- pruneindicators(scadv, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2adv)


### figs
ALTfig <- ggplot(data = BetaEnvNew, aes(x = Thaw.Stage, y = -1*MaxThaw, fill = Year))+
    geom_boxplot() +
    geom_hline(yintercept = 0) +  #surface
    geom_hline(yintercept=-46) +  #mean historical ALT
    geom_hline(yintercept=-46+9.01, linetype="dashed") + #Sd+
    geom_hline(yintercept=-46-9.01, linetype="dashed") + #Sd-
    geom_hline(yintercept = -120) #probe length
ALTfig

