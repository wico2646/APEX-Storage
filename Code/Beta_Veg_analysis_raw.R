BetaVeg17 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVeg17.csv", row.names=1)
View(BetaVeg17)
BetaVeg18 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVeg18.csv", row.names=1)
View(BetaVeg18)
BetaVeg19 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVeg19.csv", row.names=1)
View(BetaVeg19)
BetaVegCombo <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVegCombo.csv", row.names=1)
View(BetaVegCombo)
BetaVegBO <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVegBO.csv", row.names=1)
View(BetaVegBO)
BetaVegBryo<-read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVegBryo.csv", row.names=1)
View(BetaVegBryo)
BetaVegNoBryo<-read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaVegNoBryo.csv", row.names=1)
View(BetaVegNoBryo)
install.packages("vegan")
library(vegan)
install.packages("MASS")
library(MASS)
library(phyloseqGraphTest)
library(ggplot2)
library(plyr)
library(dplyr)
library(installr)
updateR()

beta17.dist<-vegdist(BetaVeg17)
beta18.dist<-vegdist(BetaVeg18)

beta17.mds<-metaMDS(BetaVeg17, trace=FALSE)
beta18.mds<-metaMDS(BetaVeg18, trace=FALSE)

beta17.mds
beta18.mds

plot(beta17.mds, type="t")
plot(beta18.mds, type="t")

beta17.pca<-rda(BetaVeg17, scale=TRUE)
beta18.pca<-rda(BetaVeg18, scale=TRUE)
beta19.pca<-rda(BetaVeg19, scale=TRUE)
betacombo.pca<-rda(BetaVegCombo, scale=TRUE)

beta17.pca
beta18.pca
beta19.pca
betacombo.pca

plot(beta17.pca)
plot(beta18.pca)
plot(beta19.pca)
plot(betacombo.pca, scaling=3)
text(betacombo.pca, display = "sites", scaling = 3, cex = 0.8, col = "darkcyan")
library(zoom)
zm(type="zoom in")
zm

biplot(beta17.pca)
biplot(beta18.pca)

mod17<-decorana(BetaVeg17)
mod18<-decorana(BetaVeg18)

plot(mod17)
plot(mod18)

BetaEnv17 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaEnv17.csv", row.names=1)
View(BetaEnv17)

BetaEnv18 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaEnv18.csv", row.names=1)
View(BetaEnv18)

BetaEnvCombo <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaEnvCombo.csv", row.names=1)
View(BetaEnvCombo)

BetaEnvBO <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaEnvBO.csv", row.names=1)
View(BetaEnvBO)

BetaEnvNew <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaEnvNew.csv", row.names=1)
View(BetaEnvNew)
BetaEnvNew$Year<-factor(BetaEnvNew$Year)

plot(beta17.mds, type="t", display="sites")
plot(beta18.mds, type="t", display="sites")
plot(betacombo.pca, display="sites")

beta17.pcscores<-scores(beta17.pca)
beta18.pcscores<-scores(beta18.pca)
betacombo.pcscores<-scores(betacombo.pca)

beta17.pcscores
beta18.pcscores
betacombo.pcscores

icechangemod<-lm(BetaEnvNew$IceChange17to19~betacombo.pcscores$sites)
summary(icechangemod)

thawmod<-lm(BetaEnv18$Sep18.DepthToIce~beta18.pcscores$sites)
summary(thawmod)

thawmod2<-lm(BetaEnv18$Sep18.DepthToIce~BetaEnv18$Location+beta18.pcscores$sites)
summary(thawmod2)

thawmbo<-lm(BetaEnvCombo$Sep18.DepthToIce~BetaEnvCombo$Location+betacombo.pcscores$sites)
summary(thawmbo)

plantonly<-lm(BetaEnvCombo$Sep18.DepthToIce~betacombo.pcscores$sites)
summary(plantonly)

methanemod1<-lm(BetaEnvNew$JulyPeakMeth.umol.m2.min~betacombo.pcscores$sites)
summary(methanemod1)
icemod1<-lm(BetaEnvNew$MaxThaw~betacombo.pcscores$sites)
summary(icemod1)
co2mod<-lm(BetaEnvNew$PeakER~betacombo.pcscores$sites)
summary(co2mod)

methice<-lm(BetaEnvCombo$JulyMethane~BetaEnvCombo$Sep18.DepthToIce)
summary(methice)

methanemod2<-lm(BetaEnvCombo$JulyMethane~betacombo.pcscores$sites+BetaEnvCombo$Sep18.DepthToIce+
                    BetaVegCombo$OXMI+BetaVegCombo$CHCA
                +BetaVegCombo$ERCH+BetaVegCombo$VAVI+BetaVegCombo$VAUL+BetaVegCombo$RUCH+BetaVegCombo$LEDE
                +BetaVegCombo$PLSC+BetaVegCombo$SPRI+BetaVegCombo$SPFU+BetaVegCombo$SPGI+BetaVegCombo$SPWA
                +BetaVegCombo$SPCA)
summary(methanemod2)
stepAIC(methanemod2, direction = "both", trace=FALSE )

methanemod3<-lm(BetaEnvCombo$JulyMethane~BetaEnvCombo$Sep18.DepthToIce+ 
                    BetaVegCombo$ERCH + BetaVegCombo$SPFU + BetaVegCombo$SPWA)
summary(methanemod3)

methanemod4<-lm(BetaEnvCombo$JulyMethane~BetaEnvCombo$Location)
summary(methanemod4)
anova(methanemod4)
plot(methanemod4)
plot(BetaEnvCombo$Location,BetaEnvCombo$JulyMethane)

methanemod5<-lm(BetaEnvCombo$JulyMethane~betacombo.pcscores$sites+BetaEnvCombo$Location+
                    BetaEnvCombo$Sep18.DepthToIce+BetaVegCombo$OXMI+BetaVegCombo$CHCA
                +BetaVegCombo$ERCH+BetaVegCombo$VAVI+BetaVegCombo$VAUL+BetaVegCombo$RUCH+
                    BetaVegCombo$LEDE+BetaVegCombo$PLSC+BetaVegCombo$SPRI+BetaVegCombo$SPFU+
                    BetaVegCombo$SPGI+BetaVegCombo$SPWA+BetaVegCombo$SPCA)
summary(methanemod5)


keyspmod<-lm(BetaEnvCombo$Sep18.DepthToIce~betacombo.pcscores$sites+BetaVegCombo$OXMI+BetaVegCombo$CHCA
             +BetaVegCombo$ERCH+BetaVegCombo$VAVI+BetaVegCombo$VAUL+BetaVegCombo$RUCH+BetaVegCombo$LEDE
             +BetaVegCombo$PLSC+BetaVegCombo$SPRI+BetaVegCombo$SPFU+BetaVegCombo$SPGI+BetaVegCombo$SPWA
             +BetaVegCombo$SPCA)
summary(keyspmod)
stepAIC(keyspmod, direction = "both", trace=FALSE )

keyspmod2<-lm(BetaEnvCombo$Sep18.DepthToIce~betacombo.pcscores$sites+BetaVegCombo$OXMI+BetaVegCombo$CHCA
              +BetaVegCombo$ERCH+BetaVegCombo$VAUL+BetaVegCombo$LEDE
              +BetaVegCombo$SPRI+BetaVegCombo$SPCA)
summary(keyspmod2)

keyspmod3<-lm(BetaEnvCombo$Sep18.DepthToIce~betacombo.pcscores$sites+BetaVegCombo$OXMI+BetaVegCombo$CHCA
              +BetaVegCombo$ERCH+BetaVegCombo$VAVI+BetaVegCombo$VAUL+BetaVegCombo$RUCH+BetaVegCombo$LEDE
              +BetaVegCombo$PLSC+BetaVegCombo$SPRI+BetaVegCombo$SPFU+BetaVegCombo$SPGI+BetaVegCombo$SPWA
              +BetaVegCombo$SPCA+BetaEnvCombo$Location*BetaEnvCombo$Year)
summary(keyspmod3)




nonlinmod1<-lm(BetaEnvCombo$Sep18.DepthToIce~betacombo.pcscores$sites+betacombo.pcscores$sites^2)
summary(nonlinmod1)

plot(nonlinmod1)


plot(beta17.pca)
plot(beta18.pca)
plot(beta18.pca, display="sites")
plot(beta18.pca, display="species")

mod18<-(beta18.pca)
with(BetaEnv18, levels(Location))
scl<-3
colvec<-c("red2", "green4", "mediumblue")
plot(mod18, type="n", scaling=scl)
with(BetaEnv18, points(mod18, display = "sites", col = colvec[Location],
                       scaling = scl, pch = 21, bg = colvec[Location]))
text(mod18, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(BetaEnv18, legend("topright", legend = levels(Location), bty = "n",
                       col = colvec, pch = 21, pt.bg = colvec))



efBeta<-envfit(beta18.pca, BetaEnv18$Sep17to18, permutations = 0, choices=c(1,2), 
               display = "sites", w  = weights(beta18.pca), na.rm = TRUE)
plot(efBeta)

efBetaSep<-envfit(beta18.pca, BetaEnv18$Sep18.DepthToIce, permutations = 0, choices=c(1,2), 
                  display = "sites", w  = weights(beta18.pca), na.rm = TRUE)
plot(efBetaSep)

mod17<-(beta17.pca)
with(BetaEnv18, levels(Location))
scl<-3
colvec<-c("red2", "green4", "mediumblue")
plot(mod17, type="n", scaling=scl)
with(BetaEnv18, points(mod17, display = "sites", col = colvec[Location],
                       scaling = scl, pch = 21, bg = colvec[Location]))
text(mod17, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(BetaEnv18, legend("topright", legend = levels(Location), bty = "n",
                       col = colvec, pch = 21, pt.bg = colvec))
with(BetaEnv18, ordiellipse(beta17.pca, Location, scaling=scl, kind = "ehull", conf = 0.95))



efBeta<-envfit(beta17.pca, BetaEnv18$Sep17to18, permutations = 0, choices=c(1,2), 
               display = "sites", w  = weights(beta18.pca), na.rm = TRUE)
plot(efBeta)

efBetaSep<-envfit(beta17.pca, BetaEnv18$Sep18.DepthToIce, permutations = 0, choices=c(1,2), 
                  display = "sites", w  = weights(beta18.pca), na.rm = TRUE)
plot(efBetaSep)

modcombo<-(betacombo.pca)
with(BetaEnv18, levels(Location))
scl<-3
colvec<-c("red2", "green4", "mediumblue")
plot(mod18, type="n", scaling=scl)
with(BetaEnv18, points(mod18, display = "sites", col = colvec[Location],
                       scaling = scl, pch = 21, bg = colvec[Location]))
text(mod18, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(BetaEnv18, legend("topright", legend = levels(Location), bty = "n",
                       col = colvec, pch = 21, pt.bg = colvec))
with(BetaEnv18, ordiellipse(mod18, Location, scaling=scl, kind = "ehull", conf = 0.95))

mod19<-(beta19.pca)
with(BetaEnv18, levels(Location))
scl<-3
colvec<-c("red2", "green4", "mediumblue")
plot(mod19, type="n", scaling=scl)
with(BetaEnv18, points(mod19, display = "sites", col = colvec[Location],
                       scaling = scl, pch = 21, bg = colvec[Location]))
text(mod19, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(BetaEnv18, legend("topright", legend = levels(Location), bty = "n",
                       col = colvec, pch = 21, pt.bg = colvec))
with(BetaEnv18, ordiellipse(mod19, Location, scaling=scl, kind = "ehull", conf = 0.95))

efBeta<-envfit(beta19.pca, BetaEnv18$Sep18.DepthToIce, permutations = 0, choices=c(1,2), 
               display = "sites", w  = weights(beta19.pca), na.rm = TRUE)
plot(efBeta)
plot(beta19.pca)
plot(beta19.pca, display="species")

fullmodel<-lm(BetaEnv18$Sep17to18~.,data=BetaVeg18)
stepwisemod<-stepAIC(fullmodel, direction="both", trace=FALSE)

modcombo<-(betacombo.pca)
with(BetaEnvCombo, levels(Location))

scl<-3.0
colvec<-c("red2", "green4", "mediumblue")
plot(modcombo, type="n", scaling=scl)
with(BetaEnvCombo, points(betacombo.pca, display = "sites", col = colvec[Location],
                          scaling = scl, pch = 21, bg = colvec[Location]))
text(modcombo, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")

with(BetaEnvCombo, legend("topright", legend = levels(Location), bty = "n",
                          col = colvec, pch = 21, pt.bg = colvec))
with(BetaEnvCombo, ordiellipse(betacombo.pca, Location, scaling=scl, kind = "ehull", conf = 0.95))
efBetaCombo<-envfit(betacombo.pca, BetaEnvCombo$Sep18.DepthToIce, permutations = 0, choices=c(1,2), 
                    display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaCombo)

efBetaMeth<-envfit(betacombo.pca, BetaEnvCombo$JulyMethane, permutations = 0, choices=c(1,2), 
                   display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaMeth)

plot(betacombo.pca, display="sites")
text(modcombo, display="sites")


beta.ca<-cca(BetaVegCombo)
eftest<-envfit(beta.ca, BetaEnvCombo, permutations=999)
eftest

plot(betacombo.pca, display = "sites", type = "p")
with(BetaEnvCombo, ordiellipse(betacombo.pca, Location, scaling=scl, kind = "ehull", conf = 0.95))
with(BetaEnvCombo, ordispider(betacombo.pca, Location, col = "blue", label= FALSE))
with(BetaEnvCombo, ordihull(betacombo.ca, Location, col="blue", lty=2))


anosimZone<-anosim(BetaVegCombo, BetaEnvNew$Thaw.Stage, distance="bray")
summary(anosimZone)
plot(anosimZone)


anosimward<-anosim(BetaVegCombo, BetaEnvNew$WardCluster, distance="bray")
summary(anosimward)
plot(anosimward)

anosimTime<-anosim(BetaVegCombo, BetaEnvNew$Year, distance="bray")
summary(anosimTime)
plot(anosimTime)

anosimZoneTime<-anosim(BetaVegCombo, BetaEnvNew$Stage.Year, distance="bray")
summary(anosimZoneTime)
plot(anosimZoneTime)

anosimBOTime<-anosim(BetaVegBO, BetaEnvBO$Year, distance = "bray")
summary(anosimBOTime)
plot(anosimBOTime)


anosimTreat<-anosim(BetaVegCombo, BetaEnvCombo$Treatment, distance="bray")
summary(anosimTreat)

pcmodel1<-lm(betacombo.pcscores$sites~BetaEnvCombo$Location*BetaEnvCombo$Year)
summary(pcmodel1)

plot(anosimZone)
plot(anosimTime)
plot(anosimTreat)




fullvspmodel<-lm(BetaEnvNew$MaxThaw~BetaVegCombo$BEGL+BetaVegCombo$CACA+BetaVegCombo$CAAQ+
                     BetaVegCombo$CAGY+BetaVegCombo$CHCA+BetaVegCombo$DRRO+BetaVegCombo$EQSC+BetaVegCombo$ERCH
                 +BetaVegCombo$ERVA+BetaVegCombo$LALA+BetaVegCombo$LEDE+BetaVegCombo$LEGR+BetaVegCombo$OXMI+
                     BetaVegCombo$PIMA+BetaVegCombo$RUCH+BetaVegCombo$VAUL+BetaVegCombo$VAVI)
stepAIC(fullvspmodel, direction = "both", trace=FALSE )
summary(fullvspmodel)


fullcspmodel<-lm(BetaEnvNew$MaxThaw~BetaVegCombo$AUPA+BetaVegCombo$CAGI+BetaVegCombo$CAST1+
                     BetaVegCombo$CAST2+BetaVegCombo$CLAM+BetaVegCombo$CLCE+BetaVegCombo$CLCR+BetaVegCombo$CLGR
                 +BetaVegCombo$CLMI+BetaVegCombo$DIUN+BetaVegCombo$HYSP+BetaVegCombo$MYAN+BetaVegCombo$PEAP+
                     BetaVegCombo$PENE+BetaVegCombo$PLSC+BetaVegCombo$POST+BetaVegCombo$PTCR+BetaVegCombo$RHPS
                 +BetaVegCombo$SPAN+BetaVegCombo$SPCA+BetaVegCombo$SPFU+BetaVegCombo$SPGI+BetaVegCombo$SPMA+
                     BetaVegCombo$SPRI+BetaVegCombo$SPSQ+BetaVegCombo$SPWA+BetaVegCombo$TONI)
stepAIC(fullcspmodel, direction = "both", trace=FALSE )
summary(fullcspmodel)

fullpspmodel<-lm(BetaEnv18$Sep18.DepthToIce~BetaVeg18$EQSC+BetaVeg18$LEDE+BetaVeg18$OXMI+
                     BetaVeg18$RUCH+BetaVeg18$VAUL+BetaVeg18$VAVI+BetaVeg18$CLAM+BetaVeg18$PLSC
                 +BetaVeg18$SPCA+BetaVeg18$SPRI+BetaVeg18$SPWA+BetaVeg18$TONI)
step(fullpspmodel, direction = "both", trace=FALSE)
summary(fullpspmodel)

fullmodel<-lm(BetaEnv18$Sep18.DepthToIce~BetaEnv18$Zone+BetaEnv18$Location+BetaEnv18$Treatment+
                  BetaVeg18$EQSC+BetaVeg18$OXMI+BetaVeg18$VAUL+BetaVeg18$VAVI+BetaVeg18$SPCA
              +BetaVeg18$SPRI)
step(fullmodel, direction="both", trace=FALSE)

anovamod<-(aov(BetaEnv18$Sep18.DepthToIce~BetaEnv18$Zone+BetaEnv18$Location+BetaEnv18$Treatment))
summary(anovamod)

anovamod2<-(aov(BetaEnv18$Sep18.DepthToIce~BetaEnv18$Zone+
                    BetaVeg18$OXMI+BetaVeg18$SPRI))
summary(anovamod2)

oxmimod<-lme(OXMI~Location, data=tester,
             random = ~ 1 | Year,
             method="ML", na.action = na.omit)
summary(oxmimod)
anova(oxmimod)

oxmimod<-lm(OXMI~Location*Year, data=tester,  na.action = na.omit)
summary(oxmimod)
anova(oxmimod)

plot(tester$Location, tester$OXMI)

library("nlme")
library("PerformanceAnalytics")
chart.Correlation(BetaVegCombo[1:63], histogram=FALSE, method='spearman')
cor(BetaVegCombo[1:63], method='spearman')





###new analysis
BetaEnvNew <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/BetaEnvNew.csv", row.names=1)
View(BetaEnvNew)
BetaEnvNew$Year<-factor(BetaEnvNew$Year)
BetaEnvNew$WardCluster<-factor(BetaEnvNew$WardCluster)

methmodnew<-lm(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
summary(methmodnew)
anova(methmodnew)
plot(BetaEnvNew$JulyPeakMeth~BetaEnvNew$Location*BetaEnvNew$Year)

methcheck<-lm(BetaEnvNew$JulyPeakMeth~BetaEnvNew$Treatment+BetaEnvNew$Thaw.Stage*BetaEnvNew$Year
              +BetaEnvNew$MaxThaw)
summary(methcheck)
anova(methcheck)
plot(BetaEnvNew$JulyPeakMeth~BetaEnvNew$Treatment+BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
step(methcheck, direction="both", trace=FALSE)
anova(step(methcheck, direction="both", trace=FALSE))


icemodnew<-lm(BetaEnvNew$MaxThaw~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
summary(icemodnew)
anova(icemodnew)
plot(BetaEnvNew$MaxThaw~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
iceanova<-aov(BetaEnvNew$MaxThaw~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
TukeyHSD(iceanova)

icemod2<-lm(BetaEnvNew$SepIceDepth~BetaEnvNew$Stage.Year)
summary(icemod2)
anova(icemod2)

clusticemod<-lm(BetaEnvNew$SepIceDepth~BetaEnvNew$WardCluster*BetaEnvNew$Thaw.Stage)
summary(clusticemod)
anova(clusticemod)


icechangemodnew<-lm(BetaEnvNew$IceChange17to19~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
summary(icechangemodnew)
anova(icechangemodnew)

neemod<-lm(BetaEnvNew$JulyNEE~BetaEnvNew$Location*BetaEnvNew$Year)
summary(neemod)
anova(neemod)
plot(BetaEnvNew$JulyNEE~BetaEnvNew$Location*BetaEnvNew$Year)

ermod<-lm(BetaEnvNew$PeakER~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
summary(ermod)
anova(ermod)
eranova<-aov(BetaEnvNew$PeakER~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
plot(BetaEnvNew$PeakER~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
TukeyHSD(eranova)

richmod<-lm(BetaEnvNew$Richness~BetaEnvNew$Location*BetaEnvNew$Year)
summary(richmod)
anova(richmod)
plot(BetaEnvNew$Richness~BetaEnvNew$Location)

carbonmod<-lm(BetaEnvNew$TotalC~BetaEnvNew$Location*BetaEnvNew$Year)
summary(carbonmod)
anova(carbonmod)
plot(BetaEnvNew$TotalC~BetaEnvNew$Location*BetaEnvNew$Year)

newbetamod<-lm(BetaEnvNew$JulyPeakMeth~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year+BetaEnvNew$Soil.moisture.mV
               +BetaEnvNew$Soil.temp.C+BetaEnvNew$Treatment)
summary(newbetamod)
anova(newbetamod)

methmod<-lm(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
summary(methmod)
anova(methmod)
methanova<-aov(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
plot(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$Thaw.Stage*BetaEnvNew$Year)
TukeyHSD(eranova)

mod1<-lm(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$MaxThaw)
summary(mod1)

###just the BO
View(BetaEnvBO)

BOicemod<-lm(BetaEnvBO$SepIceDepth~BetaEnvBO$Year)
summary(BOicemod)
boxplot(BetaEnvBO$SepIceDepth~BetaEnvBO$Year)

BOmethmod<-lm(BetaEnvBO$JulyPeakMeth~BetaEnvBO$Year)
summary(BOmethmod)
boxplot(BetaEnvBO$JulyPeakMeth~BetaEnvBO$Year)

BOneemod<-lm(BetaEnvBO$JulyNEE~BetaEnvBO$Year)
summary(BOneemod)
boxplot(BetaEnvBO$JulyNEE~BetaEnvBO$Year)

BOermod<-lm(BetaEnvBO$JulyER~BetaEnvBO$Year)
summary(BOermod)
boxplot(BetaEnvBO$JulyER~BetaEnvBO$Year)

### Key species
BOOXMImod<-lm(BetaVegBO$OXMI~BetaEnvBO$Year)
summary(BOOXMImod)

BOCHCAmod<-lm(BetaVegBO$CHCA~BetaEnvBO$Year)
summary(BOCHCAmod)



###PCA new
scl<-3.0
colvec<-c("red2", "green4", "mediumblue", "purple")
colvec2<-c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080")
plot(modcombo, type="n", scaling=scl)
with(BetaEnvNew, points(betacombo.pca, display = "sites", col = colvec3[Stage.Year],
                        scaling = scl, pch = 21, bg = colvec3[Stage.Year]))
text(modcombo, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
text(modcombo, display = "sites", scaling = scl, cex = 0.8, col="darkcyan")

with(BetaEnvNew, legend("topright", legend = levels(Location), bty = "n",
                        col = colvec2, pch = 21, pt.bg = colvec2))
with(BetaEnvNew, ordiellipse(betacombo.pca, Location, scaling=scl, kind = "ehull", conf = 0.95))
efBetaIce<-envfit(betacombo.pca, BetaEnvNew$SepIceDepth, permutations = 0, choices=c(1,2), 
                  display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaIce)

efBetaIceChange<-envfit(betacombo.pca, BetaEnvNew$IceChange17to19, permutations = 0, choices=c(1,2), 
                        display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaIceChange)

efBetaMethPC<-envfit(betacombo.pca, BetaEnvNew$JulyPeakMeth, permutations = 0, choices=c(1,2), 
                     display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaMethPC)

efBetaMoistPC<-envfit(betacombo.pca, BetaEnvNew$Soil.moisture.mV, permutations = 0, choices=c(1,2), 
                      display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaMoistPC)

efBetaTempPC<-envfit(betacombo.pca, BetaEnvNew$Soil.temp.C, permutations = 0, choices=c(1,2), 
                     display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaTempPC)


###subset by year
env<-BetaEnvNew
theme_set(theme_bw())


plot(betacombo.pcscores$sites, col=colvec2[BetaEnvNew$Location], pch=21, 
     bg=colvec2[BetaEnvNew$Location])




pchvec<-c(1,10,19)
colvec2<-c("#96d6e9","#000080","#c016bc","#ff8000")
colvec3<-c("#dcf1f8", "#b9e3f0", "#96d6e9", "#0000ff", "#0000bf", "#000080", 
           "#ec53e8", "#e723e1", "#c016bc", "#ffba75", "#ff9f40", "#ff8000")
plot(modcombo, type="n", scaling=scl)
with(BetaEnvNew, points(betacombo.pca, display = "sites", col = colvec2[Location],
                        scaling = scl, pch = 21, bg = colvec2[Location]))

text(modcombo, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
text(modcombo, display = "sites", scaling = scl, cex = 0.8, col="darkcyan")

with(BetaEnvNew, legend("topright", legend = levels(Location), bty = "n",
                        col = colvec2, pch = 21, pt.bg = colvec2))
with(BetaEnvNew, ordiellipse(betacombo.pca, Location, scaling=scl, kind = "ehull", conf = 0.95))
efBetaIce<-envfit(betacombo.pca, BetaEnvNew$SepIceDepth, permutations = 0, choices=c(1,2), 
                  display = "sites", w  = weights(betacombo.pca), na.rm = TRUE)
plot(efBetaIce)



###NMDS
library(vegan)
library(ggplot2)
# Assuming 'site_name' is the name of the column with site identifiers
rownames(BetaVegCombo) <- BetaVegCombo$Plot.Year

# Remove the site name column from the dataframe
BetaVegCombo$Plot.Year <- NULL

#betanew.dist<-vegdist(BetaVegCombo)

betanew.mds<-metaMDS(BetaVegCombo, trace=FALSE)
#betaVP.mds<-metaMDS(BetaVegNoBryo)
#betaVP.mds
#plot(betaVP.mds)
#text(betaVP.mds, display="species")
#with(BetaEnvNew, ordiellipse(betaVP.mds, Thaw.Stage, scaling=scl,
#kind = "ehull", conf = 0.95))
#with(BetaEnvNew, points(betaVP.mds, display = "sites", col = colvec3[Stage.Year],
#scaling = scl, pch = 21, bg = colvec3[Stage.Year]))

###NMDS fig
View(BetaEnvNew)
betanew.mds
stressplot(betanew.mds)
BetaEnvNew$Thaw.Stage<-as.factor(BetaEnvNew$Thaw.Stage)

plot(betanew.mds)
points(betanew.mds, display="sites", pch=c(15,16,17) [as.factor(BetaEnvNew$Year)],
       col=colvec2[as.factor(BetaEnvNew$Thaw.Stage)])
legend("topright", legend=c(levels(BetaEnvNew$Year), levels(BetaEnvNew$Thaw.Stage)), 
       pch=c(15,16,17, 15, 15, 15, 15), col=c("black", "black", "black", "#96d6e9", "#000080","#c016bc","#ff8000")
       ,bty="n", cex=1)


#with(BetaEnvNew, points(betanew.mds, display = "sites", col = colvec3[Stage.Year],
#                       pch = 21, bg = colvec3[Stage.Year]))
#with(BetaEnvNew, legend("topright", legend = levels(Thaw.Stage), bty = "n",
#                       col = colvec2, pch = 21, pt.bg = colvec2))
with(BetaEnvNew, ordiellipse(betanew.mds, Thaw.Stage, scaling=scl,
                             kind = "ehull", conf = 0.90, col=colvec2))
text(betanew.mds, display="sites")
text(betanew.mds, display="species")


###Adding in environmental vectors
envsubset<-BetaEnvNew[c(15, 19, 20, 21, 26, 30, 32, 34)]
envsubset<-rename(envsubset, "ALT"=MaxThaw)
envsubset<-rename(envsubset, "Soil Moisture"=Soil.moisture.mV)
envsubset<-rename(envsubset, "Soil Temperature"=Soil.temp.C)
envsubset<-rename(envsubset, "CH4"=JulyPeakMeth.umol.m2.min)
envsubset<-rename(envsubset, "NEE"=JulyNEE)
envsubset<-rename(envsubset, "ER"=JulyER)
envsubset<-rename(envsubset, "Species Richness"=Richness)
View(envsubset)
env<-envfit(betanew.mds, envsubset, permutations=999, na.rm=TRUE)
plot(env, col="black")

efBetaIceNMDS<-envfit(betanew.mds, BetaEnvNew$MaxThaw, permutations = 999, choices=c(1,2),
                      display = "sites", w  = weights(betanew.mds), na.rm = TRUE)
plot(efBetaIceNMDS, col="black")
efBetaMethNMDS<-envfit(betanew.mds, BetaEnvNew$JulyPeakMeth.umol.m2.min, permutations = 0, choices=c(1,2), 
                       display = "sites", w  = weights(betanew.mds), na.rm = TRUE)
plot(efBetaMethNMDS, col="black")
efBetaMoistNMDS<-envfit(betanew.mds, BetaEnvNew$Soil.moisture.mV, permutations = 0, choices=c(1,2), 
                        display = "sites", w  = weights(betanew.mds), na.rm = TRUE)
plot(efBetaMoistNMDS, col="black")
efBetaTempNMDS<-envfit(betanew.mds, BetaEnvNew$Soil.temp.C, permutations = 0, choices=c(1,2), 
                       display = "sites", w  = weights(betanew.mds), na.rm = TRUE)
plot(efBetaTempNMDS, col="black")
efBetaERNMDS<-envfit(betanew.mds, BetaEnvNew$JulyER, permutations = 0, choices=c(1,2), 
                     display = "sites", w  = weights(betanew.mds), na.rm = TRUE)
plot(efBetaERNMDS, col="black")
efBetaNEENMDS<-envfit(betanew.mds, BetaEnvNew$JulyNEE, permutations = 0, choices=c(1,2), 
                      display = "sites", w  = weights(betanew.mds), na.rm = TRUE)
plot(efBetaNEENMDS, col="black")

#New NMDS
#site.scrs<-as.data.frame(scores(betanew.mds, display="sites"))
#site.scrs<-cbind(site.scrs, Thaw.Stage=BetaEnvNew$Thaw.Stage)
#site.scrs<-cbind(site.scrs, Year=BetaEnvNew$Year)
#head(site.scrs)

#nmds.plot.beta<-ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+
# geom_point(aes(NMDS1, NMDS2, colour=factor(Thaw.Stage), shape=factor(Year)), size=2)+
#coord_fixed()+
#theme_classic()+
#theme(panel.background=element_rect(fill=NA, colour="black", size=1, linetype="solid)"))+
#labs(colour="Thaw Stage", shape="Year")+
#theme(legend.position="right", legend.text=element_text(size=12), legend.title=element_text(size=12),
#     axis.text = element_text(size=10))
#nmds.plot.beta


###cluster NMDS
plot(betanew.mds)
with(BetaEnvNew, points(betanew.mds, display = "sites", col = colvec2[WardCluster],
                        scaling = scl, pch = 21, bg = colvec2[WardCluster]))
with(BetaEnvNew, legend("topright", legend = levels(WardCluster), bty = "n",
                        col = colvec2, pch = 21, pt.bg = colvec2))
with(BetaEnvNew, ordiellipse(betanew.mds, Thaw.Stage, scaling=scl,
                             kind = "ehull", conf = 0.95))
with(BetaEnvNew, ordiellipse(betanew.mds, WardCluster, scaling=scl,
                             kind = "ehull", conf = 0.95))
text(betanew.mds, display="sites")
text(betanew.mds, display="species")

clustmod<-lm(BetaEnvNew$Stage~BetaEnvNew$WardCluster)
summary(clustmod)

clustmod2<-lm(BetaEnvNew$SepIceDepth~BetaEnvNew$WardCluster)
summary(clustmod2)
anova(clustmod2)

clustmod3<-lm(BetaEnvNew$JulyPeakMeth~BetaEnvNew$WardCluster)
summary(clustmod3)

###just lichens and bryophytes
betabryo.dist<-vegdist(BetaVegBryo)

betabryo.mds<-metaMDS(BetaVegBryo, trace=FALSE)

betabryo.mds
stressplot(betabryo.mds)

plot(betabryo.mds)
with(BetaEnvNew, points(betabryo.mds, display = "sites", col = colvec3[Stage.Year],
                        scaling = scl, pch = 21, bg = colvec3[Stage.Year]))
with(BetaEnvNew, legend("topright", legend = levels(Thaw.Stage), bty = "n",
                        col = colvec2, pch = 21, pt.bg = colvec2))
with(BetaEnvNew, ordiellipse(betabryo.mds, Thaw.Stage, scaling=scl,
                             kind = "ehull", conf = 0.95))
text(betanew.mds, display="sites")
text(betanew.mds, display="species")

###visuals
IceGraph<- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/IceGraph.csv", row.names=1)
View(IceGraph)
CGraph<- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/CGraph.csv", row.names=1)
View(CGraph)

plot(BetaEnvNew$Year, BetaEnvNew$SepIceDepth, col=colvec2[BetaEnvNew$Location], pch=21, 
     bg=colvec2[BetaEnvNew$Location])

with(BetaEnvNew, legend("topright", legend = levels(BetaEnvNew$Location), bty = "n",
                        col = colvec, pch = 21, pt.bg = colvec))

BetaEnvNew$Year<-factor(BetaEnvNew$Year)
CGraph$Year<-factor(CGraph$Year)
IceGraph$Year<-factor(IceGraph$Year)
VegGraph$Year<-factor(VegGraph$Year)
ALData.Means.Time$Year<-factor(ALData.Means.Time$Year)

##AL fig
HistAL <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/HistAL.csv")
HistData.Means.Time <-summarySE(data = HistAL,
                                groupvars =c("Treat"), 
                                measurevar="AL.depth", na.rm=TRUE)
View(HistData.Means.Time)

ALData.Means.Time <-summarySE(data = BetaEnvNew,
                              groupvars =c("Thaw.Stage"), 
                              measurevar="MaxThaw", na.rm=TRUE)
View(ALData.Means.Time)
###bargraph
ggplot(ALData.Means.Time, aes(x=Thaw.Stage, y=MaxThaw, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=MaxThaw-ci, ymax=MaxThaw+ci),
                  position=position_dodge(0.9), width=0.2)+
    geom_hline(yintercept=46)+
    geom_hline(yintercept=46+9.01, linetype="dashed")+
    geom_hline(yintercept=46-9.01, linetype="dashed")


ggplot(CGraph, aes(x=Location, y=MeanC.g.min.m2, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=MeanC.g.min.m2-SDC.g.min.m2, ymax=MeanC.g.min.m2+SDC.g.min.m2),
                  position=position_dodge(0.9), width=0.2)

ggplot(CGraph, aes(x=Location, y=MeanC.g.day.m2, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=MeanC.g.day.m2-SDC.g.day.m2, ymax=MeanC.g.day.m2+SDC.g.day.m2),
                  position=position_dodge(0.9), width=0.2)

ggplot(CGraph, aes(x=Location, y=CO2eq, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=CO2eq-SDCO2eq, ymax=CO2eq+SDCO2eq),
                  position=position_dodge(0.9), width=0.2)

###Methane Fig base
library(Rmisc)
BetaEnvNew$Year<-as.factor(BetaEnvNew$Year)
CH4Data.Means.Time <-summarySE(data = BetaEnvNew,
                               groupvars =c("Year","Thaw.Stage"), 
                               measurevar="JulyPeakMeth.umol.m2.min", na.rm=TRUE)
View(CH4Data.Means.Time)

#
#CH4Data.Means.Time.Year <-summarySE(data = BetaEnvNew,
#groupvars =c("Year"), 
#measurevar="JulyPeakMeth.umol.m2.min", na.rm=TRUE)
#View(CH4Data.Means.Time.Year)

ggplot(CH4Data.Means.Time, aes(x=Thaw.Stage, y=JulyPeakMeth.umol.m2.min, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=JulyPeakMeth.umol.m2.min-ci, ymax=JulyPeakMeth.umol.m2.min+ci),
                  position=position_dodge(0.9), width=0.2)
###controls only
BetaEnvControlsOnly <- read.csv("C:/Users/danci/Desktop/BetaEnvControlsOnly.csv")
View(BetaEnvControlsOnly)

BetaEnvControlsOnly$Year<-as.factor(BetaEnvControlsOnly$Year)
CH4Data.Means.TimeC <-summarySE(data = BetaEnvControlsOnly,
                                groupvars =c("Thaw.Stage"), 
                                measurevar="JulyPeakMeth.umol.m2.min", na.rm=TRUE)
View(CH4Data.Means.TimeC)

ggplot(CH4Data.Means.TimeC, aes(x=Thaw.Stage, y=JulyPeakMeth.umol.m2.min))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=JulyPeakMeth.umol.m2.min-ci, ymax=JulyPeakMeth.umol.m2.min+ci),
                  position=position_dodge(0.9), width=0.2)


#Methane season fig
ggplot(CGraph, aes(x=Stage, y=Cum.CH4.umol, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Cum.CH4.umol-SD.cum.CH4, ymax=Cum.CH4.umol+SD.cum.CH4),
                  position=position_dodge(0.9), width=0.2)

###CO2 fig base
CO2Data.Means.Time <-summarySE(data = BetaEnvNew,
                               groupvars =c("Year","Thaw.Stage"), 
                               measurevar="JulyER", na.rm=TRUE)
View(CO2Data.Means.Time)

ggplot(CO2Data.Means.Time, aes(x=Thaw.Stage, y=JulyER, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=JulyER-ci, ymax=JulyER+ci),
                  position=position_dodge(0.9), width=0.2)
#ER season fig
ggplot(CGraph, aes(x=Stage, y=Cum.ER.umol, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Cum.ER.umol-SD.cum.ER, ymax=Cum.ER.umol+SD.cum.ER),
                  position=position_dodge(0.9), width=0.2)

NEEData.Means.Time <-summarySE(data = BetaEnvNew,
                               groupvars =c("Year","Thaw.Stage"), 
                               measurevar="JulyNEE", na.rm=TRUE)
View(NEEData.Means.Time)
ggplot(NEEData.Means.Time, aes(x=Thaw.Stage, y=JulyNEE, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=JulyNEE-ci, ymax=JulyNEE+ci),
                  position=position_dodge(0.9), width=0.2)
###total C fig
ggplot(CGraph, aes(x=Stage, y=Total.C, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Total.C-SDC, ymax=Total.C+SDC),
                  position=position_dodge(0.9), width=0.2)

#total c season
ClossData.Means.Time <-summarySE(data = BetaEnvNew,
                                 groupvars =c("Year","Thaw.Stage"), 
                                 measurevar="Closs.gC.m2", na.rm=TRUE)
View(ClossData.Means.Time)
ggplot(CGraph, aes(x=Stage, y=Total.C.umol, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Total.C.umol-SD.total.C, ymax=Total.C.umol+SD.total.C),
                  position=position_dodge(0.9), width=0.2)
#ER vs CH4
plot(BetaEnvNew$Cum.CH4.umol, BetaEnvNew$Cum.ER.umol, col = colvec2[BetaEnvNew$Thaw.Stage],
     pch = 21, bg = colvec2[BetaEnvNew$Thaw.Stage])
ERCH4mod<-lm(BetaEnvNew$Cum.ER.umol~BetaEnvNew$Cum.CH4.umol+BetaEnvNew$Cum.CH4.umol^2)

fit2<-lm(BetaEnvNew$Cum.ER.umol~poly(BetaEnvNew$Cum.CH4.umol,2,raw=TRUE))
summary(fit2)

fit2$coefficient[1]
fit2$coefficient[2]
fit2$coefficient[3]

quadratic = fit2$coefficient[3]*BetaEnvNew$Cum.CH4.umol^2 + 
    fit2$coefficient[2]*BetaEnvNew$Cum.CH4.umol + fit2$coefficient[1]
quadratic

plot(BetaEnvNew$Cum.CH4.umol, BetaEnvNew$Cum.ER.umol, main="Scatterplot", 
     xlab="Cumulative CH4 (umol)", ylab="Cumulative ER (umol)", 
     pch=19, col=colvec2[BetaEnvNew$Thaw.Stage])
par(new = TRUE)
lines(BetaEnvNew$Cum.CH4.umol,quadratic, col="red")

##Loess fit
library(lattice)

xyplot(BetaEnvNew$Cum.ER.umol~BetaEnvNew$Cum.CH4.umol, type=c("smooth", "p"),
       pch=19, col=colvec3d[BetaEnvNew$Thaw.Stage])

CH4vsER<-lm(BetaEnvNew$JulyER~BetaEnvNew$JulyPeakMeth.umol.m2.min)
summary(CH4vsER)
plot(BetaEnvNew$JulyPeakMeth.umol.m2.min,BetaEnvNew$JulyER, pch=19, col="black")
plot(BetaEnvNew$JulyPeakMeth.umol.m2.min,BetaEnvNew$JulyER, pch=19, col=colvec3d[BetaEnvNew$Thaw.Stage])
abline(CH4vsER)

CH4vsER2<-lm(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$JulyER)
summary(CH4vsER2)
plot(BetaEnvNew$JulyER,BetaEnvNew$JulyPeakMeth.umol.m2.min, pch=19, col="black")
plot(BetaEnvNew$JulyER,BetaEnvNew$JulyPeakMeth.umol.m2.min, pch=19, col=colvec3d[BetaEnvNew$Thaw.Stage])
abline(CH4vsER2)

library(ggplot2)
ggplot(BetaEnvNew, aes(JulyPeakMeth.umol.m2.min,JulyER))+
    geom_point(col=colvec3d[BetaEnvNew$Thaw.Stage])

ggplot(BetaEnvNew, aes(JulyPeakMeth.umol.m2.min,JulyER))+
    geom_point(col=colvec3d[BetaEnvNew$Thaw.Stage])+
    stat_ellipse(aes(x=JulyPeakMeth.umol.m2.min, y=JulyER,color=colvec3d[BetaEnvNew$Thaw.Stage], 
                     group=Thaw.Stage),type = "norm")
##other figs
plot(BetaEnvNew$JulyPeakMeth.umol.m2.min, BetaEnvNew$JulyER, col=colvec2[BetaEnvNew$Thaw.Stage], 
     pch = 21, bg = colvec2[BetaEnvNew$Thaw.Stage])
plot(BetaEnvNew$Cum.CH4.umol, BetaEnvNew$GPP, col=colvec2[BetaEnvNew$Thaw.Stage], 
     pch = 21, bg = colvec2[BetaEnvNew$Thaw.Stage])

##surface plot
library(devtools)
library(rlang)
devtools::install_github("AckerDWM/gg3D")

library("gg3D")
install.packages("scatterplot3d") # Install
library("scatterplot3d")

data(BetaEnvNew)
ALdepth<-c(-1*BetaEnvNew$MaxThaw)
colvec3d<-c("#0F2080","#85C0F9","#A95AA1","#F5793A")
colvec3d <- colvec3d[as.numeric(BetaEnvNew$Thaw.Stage)]
s3d<-scatterplot3d(BetaEnvNew$JulyPeakMeth.umol.m2.min, BetaEnvNew$JulyNEE, ALdepth,
                   pch=16, color=colvec3d, type="h")
legend("right", legend = levels(BetaEnvNew$Thaw.Stage),
       col =  c("#0F2080","#85C0F9","#A95AA1","#F5793A"), pch = 16)
my.lm <- lm(ALdepth ~ BetaEnvNew$JulyPeakMeth.umol.m2.min + BetaEnvNew$JulyER)
s3d$plane3d(my.lm)

s3d2<-scatterplot3d(BetaEnvNew$Cum.CH4.gC, BetaEnvNew$Cum.NEE.gC, ALdepth,
                    pch=16, color=colvec3d, type="h")
legend("right", legend = levels(BetaEnvNew$Thaw.Stage),
       col =  c("#0F2080","#85C0F9","#A95AA1","#F5793A"), pch = 16)
my.lm2<- lm(ALdepth ~ BetaEnvNew$Cum.CH4.gC + BetaEnvNew$Cum.NEE.gC)
s3d2$plane3d(my.lm2)

#new 3d plot
library("plotly")
library(reshape2)
plot_matrix <- t(acast(BetaEnvNew, JulyPeakMeth.umol.m2.min~JulyER, value.var="MaxThaw"))
plot_matrix
write.csv(plot_matrix, "plot_matrix.csv")

library(interp)
####this code actually works
x<-c(BetaEnvNew$JulyPeakMeth.umol.m2.min)
y<-c(BetaEnvNew$JulyNEE)
z<-c(BetaEnvNew$MaxThaw)
interp(x,y,z)
interpmatrix<-interp(x,y,z, linear=FALSE)
View(interpmatrix)
fig <- plot_ly(x = interpmatrix$x,
               y = interpmatrix$y,
               z = interpmatrix$z*-1) %>% add_surface()
fig

x<-c(BetaEnvNew$Cum.CH4.gC)
y<-c(BetaEnvNew$Cum.NEE.gC)
z<-c(BetaEnvNew$MaxThaw)
interp(x,y,z)
interpmatrix<-interp(x,y,z, linear=FALSE)
View(interpmatrix)
fig <- plot_ly(x = interpmatrix$x,
               y = interpmatrix$y,
               z = interpmatrix$z) %>% add_surface()
fig

balancemod<-lme(Closs.gC.m2~Thaw.Stage*Year, random = ~1|Treatment/PlotName, data=BetaEnvNew)
summary(balancemod)
anova.lme(balancemod)

###new 3d surface plot
xmod<-lm(BetaEnvNew$JulyPeakMeth.umol.m2.min~BetaEnvNew$MaxThaw*-1)
ymod<-lm(BetaEnvNew$JulyER~BetaEnvNew$MaxThaw*-1)
dummymod<-lm(BetaEnvNew$MaxThaw*-1~BetaEnvNew$JulyPeakMeth.umol.m2.min+BetaEnvNew$JulyER)
summary(dummymod)
plot(BetaEnvNew$JulyPeakMeth.umol.m2.min, BetaEnvNew$MaxThaw*-1)
plot(BetaEnvNew$JulyER, BetaEnvNew$MaxThaw*-1)
CH4vect<-predict()

library(rgl)
library(ggplot2)
plot3d(BetaEnvNew$JulyPeakMeth.umol.m2.min, BetaEnvNew$JulyER, BetaEnvNew$MaxThaw*-1)
ggplot(BetaEnvNew,aes(y=JulyER,x=JulyPeakMeth.umol.m2.min,color=MaxThaw*-1))
+geom_point()+stat_smooth(method="lm",se=FALSE)

#veg figs
GramsData.Means.Time <-summarySE(data = BetaEnvNew,
                                 groupvars =c("Year","Thaw.Stage"), 
                                 measurevar="Graminoids", na.rm=TRUE)
View(GramsData.Means.Time)
ggplot(GramsData.Means.Time, aes(x=Thaw.Stage, y=Graminoids, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Graminoids-ci, ymax=Graminoids+ci),
                  position=position_dodge(0.9), width=0.2)

ErisData.Means.Time <-summarySE(data = BetaEnvNew,
                                groupvars =c("Year","Thaw.Stage"), 
                                measurevar="Ericoids", na.rm=TRUE)
View(ErisData.Means.Time)
ggplot(ErisData.Means.Time, aes(x=Thaw.Stage, y=Ericoids, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Ericoids-ci, ymax=Ericoids+ci),
                  position=position_dodge(0.9), width=0.2)

SphagData.Means.Time <-summarySE(data = BetaEnvNew,
                                 groupvars =c("Year","Thaw.Stage"), 
                                 measurevar="Sphagna", na.rm=TRUE)
View(SphagData.Means.Time)
ggplot(SphagData.Means.Time, aes(x=Thaw.Stage, y=Sphagna, fill=Year))+
    geom_bar(position="dodge", stat="identity")+
    geom_errorbar(aes(ymin=Sphagna-ci, ymax=Sphagna+ci),
                  position=position_dodge(0.9), width=0.2)


####Cluster analysis
library(cluster)
View(BetaVegCombo)

wss <- (nrow(BetaVegCombo)-1)*sum(apply(BetaVegCombo,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(BetaVegCombo,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

fit<-kmeans(BetaVegCombo, 4)
aggregate(BetaVegCombo, by=list(fit$cluster), FUN=mean)
clusplot(BetaVegCombo, fit$cluster, color=TRUE, shade=TRUE,
         labels=4, lines=0)

betaclust<-hclust(betanew.dist,"ward")
plot(betaclust)
rect.hclust(fit, k=4, border="red")
cut<-cutree(betaclust,k=4)
plot(betaclust, labels = as.character(cut))

betaVP.dist<-vegdist(BetaVegNoBryo)
betaclustVP<-hclust(betaVP.dist,"ward")
plot(betaclustVP)

library(fpc)
plotcluster(BetaVegCombo, fit$cluster)

### ward cluster analysis
library(pvclust)
fit <- pvclust(BetaVegCombo, method.hclust="ward",
               method.dist="euclidean")
plot(fit)
pvrect(fit, alpha=.95)

# Ward Hierarchical Clustering
d <- dist(BetaVegCombo, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters
rect.hclust(fit, k=4, border="red") 

#comparison of Ward and predefined
WardCor<-cor.test(BetaEnvNew$StageNum, BetaEnvNew$ClustNum, method="pearson")
print(WardCor)
plot(BetaEnvNew$StageNum, BetaEnvNew$ClustNum)

wardexpl<-lm((StageNum-ClustNum)~IceChange17to19, data=BetaEnvNew)
summary(wardexpl)
plot(wardexpl)
plot((StageNum-ClustNum)~IceChange17to19, data=BetaEnvNew)

#### Mixed effects models
library(nlme)
library(multcomp)
library(multcompView)
library(emmeans)
library(rcompanion)
mixedmeth<-lme(JulyPeakMeth.umol.m2.min~Thaw.Stage*Year, random = ~1|Treatment/PlotName, data=BetaEnvNew)
summary(mixedmeth)
anova.lme(mixedmeth)

##new corrections
#methane
NALT<-interaction(BetaEnvNew$Treatment, BetaEnvNew$MaxThaw)
Nthaw<-interaction(BetaEnvNew$Treatment, BetaEnvNew$Thaw.Stage)
BetaEnvNew$Year<-factor(BetaEnvNew$Year)

correctedmixmeth<-lme(JulyPeakMeth.umol.m2.min~Thaw.Stage*Year, 
                      random=list(Nthaw=~1, PlotName=~1), 
                      data=BetaEnvNew)
summary(correctedmixmeth)
anova.lme(correctedmixmeth)

CH4.emm<-emmeans(correctedmixmeth,
                 pairwise ~ Thaw.Stage|Year,
                 adjust="tukey")
class(CH4.emm)
cld(CH4.emm$emmeans, alpha=0.05, Letters=letters, adjust="tukey")

summary(CH4.emm)


#ER
correctedmixER<-lme(JulyER~Thaw.Stage*Year, 
                    random=list(Nthaw=~1, PlotName=~1), 
                    data=BetaEnvNew)
summary(correctedmixER)
anova.lme(correctedmixER)

ER.emm<-emmeans(correctedmixER,
                pairwise ~ Thaw.Stage|Year,
                adjust="tukey")
class(ER.emm)
cld(ER.emm$emmeans, alpha=0.05, Letters=letters, adjust="tukey")

summary(ER.emm)

#GPP
correctedmixGPP<-lme(GPP~Thaw.Stage*Year, 
                     random=list(Nthaw=~1, PlotName=~1), 
                     data=BetaEnvNew)
summary(correctedmixGPP)
anova.lme(correctedmixGPP)

GPP.emm<-emmeans(correctedmixGPP,
                 pairwise ~ Thaw.Stage|Year,
                 adjust="tukey")
class(GPP.emm)
cld(GPP.emm$emmeans, alpha=0.05, Letters=letters, adjust="tukey")

summary(GPP.emm)

#NEE
correctedmixNEE<-lme(JulyNEE~Thaw.Stage*Year, 
                     random=list(Nthaw=~1, PlotName=~1), 
                     data=BetaEnvNew)
summary(correctedmixNEE)
anova.lme(correctedmixNEE)

NEE.emm<-emmeans(correctedmixNEE,
                 pairwise ~ Thaw.Stage|Year,
                 adjust="tukey")
class(NEE.emm)
cld(NEE.emm$emmeans, alpha=0.05, Letters=letters, adjust="tukey")

summary(NEE.emm)




##N scatter
#meth
ggplot(BetaEnvNew, aes(x=MaxThaw, y=JulyPeakMeth.umol.m2.min, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")

#ER
ggplot(BetaEnvNew, aes(x=MaxThaw, y=JulyER, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")

#GPP
ggplot(BetaEnvNew, aes(x=MaxThaw, y=GPP, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")

#Ericoids
ggplot(BetaEnvNew, aes(x=MaxThaw, y=Ericoids, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")
#Graminoids
ggplot(BetaEnvNew, aes(x=MaxThaw, y=Graminoids, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")
#Sphagna
ggplot(BetaEnvNew, aes(x=MaxThaw, y=Sphagna, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")
#Richness
ggplot(BetaEnvNew, aes(x=MaxThaw, y=Richness, colour=factor(Ntreat)))+
    geom_point(aes(color=factor(Ntreat)))+
    stat_smooth(method="lm")

#cumulative
seasonmeth<-lme(Cum.CH4.umol~Thaw.Stage*Year, random = ~1|Treatment/PlotName, data=BetaEnvNew)
summary(seasonmeth)
anova.lme(seasonmeth)

mixedER<-lme(PeakER~Thaw.Stage*Year, random = ~1|Treatment/PlotName, data=BetaEnvNew)
summary(mixedER)
anova.lme(mixedER)

seasonER<-lme(Cum.ER.umol~Thaw.Stage*Year, random = ~1|Treatment/PlotName, data=BetaEnvNew)
summary(seasonER)
anova.lme(seasonER)

grandmixedmeth<-lme(JulyPeakMeth.umol.m2.min~Thaw.Stage*Year+Soil.moisture.mV+
                        Soil.temp.C+MaxThaw, random = ~1|Treatment/PlotName, data=BetaEnvNew)
summary(grandmixedmeth)
anova.lme(grandmixedmeth)

newmixmeth<-lme(JulyPeakMeth.umol.m2.min~Thaw.Stage+Year
                , random = ~1|Treatment, data=BetaEnvNew)
summary(newmixmeth)
anova.lme(newmixmeth)

tempmod<-lme(JulyPeakMeth.umol.m2.min~
                 Soil.temp.C, random = ~1|Treatment, data=BetaEnvNew)
summary(tempmod)
plot(BetaEnvNew$Soil.temp.C,BetaEnvNew$JulyPeakMeth.umol.m2.min)
tempmodsimple<-lm(JulyPeakMeth.umol.m2.min~
                      Soil.temp.C, data=BetaEnvNew)
abline(tempmodsimple)

BetaEnvNew$JulyPeakMeth.umol.m2.min.log<- transformTukey(BetaEnvNew$JulyPeakMeth.umol.m2.min, plotit = FALSE)

#mixed veg groups
gramsmod<-lme(Graminoids~Thaw.Stage*Year
              , random = ~1|Treatment, data=BetaEnvNew)
summary(gramsmod)
anova.lme(gramsmod)

ericoidsmod<-lme(Ericoids~Thaw.Stage*Year
                 , random = ~1|Treatment, data=BetaEnvNew)
summary(ericoidsmod)
anova.lme(ericoidsmod)

sphagnamod<-lme(Sphagna~Thaw.Stage*Year
                , random = ~1|Treatment, data=BetaEnvNew)
summary(sphagnamod)
anova.lme(sphagnamod)

richnessmod<-lme(Richness~Thaw.Stage*Year
                 , random = ~1|Treatment, data=BetaEnvNew)
summary(richnessmod)
anova.lme(richnessmod)

NEEmod<-lme(JulyNEE~Thaw.Stage*Year
            , random = ~1|Treatment, data=BetaEnvNew)
summary(NEEmod)
anova.lme(NEEmod)

library(nlme)
#assessing models
qqnorm(resid(newmixmeth))
qqline(resid(newmixmeth))

qqnorm(resid(mixedER))
qqline(resid(mixedER))

qqnorm(resid(icemodnew))
qqline(resid(icemodnew))

qqnorm(resid(seasonmeth))
qqline(resid(seasonmeth))

qqnorm(resid(seasonER))
qqline(resid(seasonER))

#Checking what we would expect this qqplot to look like if normal with random data
qqnorm(rnorm(400, mean=0, sd=1))
qqline(rnorm(400, mean=0, sd=1))
plot(newmixmeth)
plot(newmixmeth, Treatment~resid(.))
plot(newmixmeth, PlotName~resid(.))
qqnorm (newmixmeth, ~ranef(., level=1))
qqnorm (newmixmeth, ~ranef(., level=2))
plotNormalHistogram(residuals(newmixmeth))
intervals(newmixmeth)
shapiro.test(residuals(newmixmeth))

##Tukey ladder transformation
library(desc)
library(DescToolsAddIns)
library(rcompanion)
BetaEnvNew$JulyMeth.Tuk<- transformTukey(BetaEnvNew$JulyPeakMeth.umol.m2.min, plotit = FALSE)
transformTukey(BetaEnvNew$JulyPeakMeth.umol.m2.min, plotit = FALSE)

newmixmethtukey<-lme(JulyMeth.Tuk~Thaw.Stage*Year
                     , random = ~1|Treatment, data=BetaEnvNew)
summary(newmixmethtukey)
anova.lme(newmixmethtukey)
plot(newmixmethtukey)
qqnorm (newmixmethtukey, ~ranef(., level=1))
qqnorm(resid(newmixmethtukey))
qqline(resid(newmixmethtukey))

BetaEnvNew$Cum.CH4.umol.Tuk<-transformTukey(BetaEnvNew$Cum.CH4.umol, plotit = FALSE)
seasonCH4.Tuk<-lme(Cum.CH4.umol.Tuk~Thaw.Stage*Year
                   , random = ~1|Treatment, data=BetaEnvNew)
summary(seasonCH4.Tuk)
anova.lme(seasonCH4.Tuk)
qqnorm(resid(seasonCH4.Tuk))
qqline(resid(seasonCH4.Tuk))

#transforming ER model
BetaEnvNew$JulyER.Tuk<- transformTukey(BetaEnvNew$PeakER, plotit = FALSE)
transformTukey(BetaEnvNew$PeakER, plotit = FALSE)
untransER<-lme(PeakER~Thaw.Stage*Year
               , random = ~1|Treatment, data=BetaEnvNew)
newmixERtukey<-lme(JulyER.Tuk~Thaw.Stage*Year
                   , random = ~1|Treatment, data=BetaEnvNew)
summary(newmixERtukey)
anova.lme(newmixERtukey)
plot(newmixERtukey)
qqnorm (newmixERtukey, ~ranef(., level=1))
qqnorm(resid(newmixERtukey))
qqline(resid(newmixERtukey))

BetaEnvNew$Cum.ER.umol.Tuk<-transformTukey(BetaEnvNew$Cum.ER.umol, plotit = FALSE)
seasonER.Tuk<-lme(Cum.ER.umol.Tuk~Thaw.Stage*Year
                  , random = ~1|Treatment, data=BetaEnvNew)
summary(seasonER.Tuk)
anova.lme(seasonER.Tuk)
qqnorm(resid(seasonER.Tuk))
qqline(resid(seasonER.Tuk))

BetaEnvNew$AL.Tuk<-transformTukey(BetaEnvNew$MaxThaw, plotit=FALSE)
icemodrand<-lme(MaxThaw~Thaw.Stage*Year
                , random = ~1|Treatment, data=BetaEnvNew)
summary(icemodrand)
anova.lme(icemodrand)
qqnorm(resid(icemodrand))
qqline(resid(icemodrand))

icemodTuk<-lme(AL.Tuk~Thaw.Stage*Year
               , random = ~1|Treatment, data=BetaEnvNew)
qqnorm(resid(icemodTuk))
qqline(resid(icemodTuk))
anova.lme(icemodTuk)

BetaEnvNew$NEE.Tuk<-transformTukey(BetaEnvNew$JulyNEE, plotit=FALSE)
NEEmodrand<-lme(Cum.NEE~Thaw.Stage*Year
                , random = ~1|Treatment, data=BetaEnvNew)
summary(NEEmodrand)
anova.lme(NEEmodrand)

BetaEnvNew$GPP.Tuk<-transformTukey(BetaEnvNew$GPP, plotit=FALSE)
GPPmodTuk<-lme(GPP~Thaw.Stage*Year
               , random = ~1|Treatment, data=BetaEnvNew)
qqnorm(resid(GPPmodTuk))
qqline(resid(GPPmodTuk))
summary(GPPmodTuk)
anova.lme(GPPmodTuk)

##Predictor species new analysis
#Sitewide
library(MuMIn)
MixedBetaVeg<-read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/MixedBetaVeg.csv", row.names=1)
View(MixedBetaVeg)

grandspmodel<-lme(MaxThaw~BEGL+CACA+CAAQ+
                      CAGY+CHCA+DRRO+EQSC+
                      ERCH+ERVA+LALA+LEDE+
                      LEGR+OXMI+PIMA+RUCH+
                      VAUL+VAVI+AUPA+CAGI+
                      CAST1+CAST2+CLAM+CLCE+
                      CLCR+CLGR+CLMI+DIUN+
                      HYSP+MYAN+PEAP+PENE+
                      PLSC+POST+PTCR+RHPS+
                      SPAN+SPCA+SPFU+SPGI+
                      SPMA+SPRI+SPSQ+SPWA+
                      TONI,data=MixedBetaVeg)
stepAIC(grandspmodel, direction = "both", trace=TRUE)
summary(grandspmodel)

grandmixmodel<-lme(MaxThaw~BEGL+CACA+CAAQ+
                       CAGY+CHCA+DRRO+EQSC+
                       ERCH+ERVA+LALA+LEDE+
                       LEGR+OXMI+PIMA+RUCH+
                       VAUL+VAVI+AUPA+CAGI+
                       CAST1+CAST2+CLAM+CLCE+
                       CLCR+CLGR+CLMI+DIUN+
                       HYSP+MYAN+PEAP+PENE+
                       PLSC+POST+PTCR+RHPS+
                       SPAN+SPCA+SPFU+SPGI+
                       SPMA+SPRI+SPSQ+SPWA+
                       TONI, random = list (~1|Year, ~1|Treatment, ~1|PlotName),data=MixedBetaVeg)
summary(grandmixmodel)
anova.lme(grandmixmodel)

grandpspmodel<-lme(MaxThaw ~ CAGY + CHCA + 
                       DRRO + ERVA + LALA + 
                       LEDE + VAUL + VAVI + 
                       CAST2 + CLCE + CLCR + 
                       HYSP + PEAP + PLSC + 
                       POST + SPGI + SPMA + 
                       SPRI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data=MixedBetaVeg)

summary(grandpspmodel)
anova.lme(grandpspmodel)


fullvspmodel<-lm(MaxThaw~BEGL+CACA+CAAQ+
                     CAGY+CHCA+DRRO+EQSC+ERCH
                 +ERVA+LALA+LEDE+LEGR+OXMI+
                     PIMA+RUCH+VAUL+VAVI, data=MixedBetaVeg)
stepAIC(fullvspmodel, direction = "both", trace=FALSE )
summary(fullvspmodel)

pvspmodel<-lme(MaxThaw~CAAQ+CHCA+DRRO+
                   LEDE+VAUL+VAVI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data=MixedBetaVeg)
summary(pvspmodel)
anova.lme(pvspmodel)


fullcspmodel<-lm(MaxThaw~AUPA+CAGI+CAST1+
                     CAST2+CLAM+CLCE+CLCR+CLGR
                 +CLMI+DIUN+HYSP+MYAN+PEAP+
                     PENE+PLSC+POST+PTCR+RHPS
                 +SPAN+SPCA+SPFU+SPGI+SPMA+
                     SPRI+SPSQ+SPWA+TONI, data=MixedBetaVeg)
stepAIC(fullcspmodel, direction = "both", trace=FALSE )
summary(fullcspmodel)

pcspmodel<-lme(MaxThaw~CLAM+DIUN+PENE+
                   SPAN+SPCA+SPGI+SPRI+
                   SPWA, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data=MixedBetaVeg)
summary(pcspmodel)
anova.lme(pcspmodel)

fullpspmodel<-lm(MaxThaw~CAAQ+CHCA+DRRO+
                     LEDE+VAUL+VAVI+CLAM+
                     DIUN+PENE+SPAN+SPCA+
                     SPGI+SPRI+SPWA, data=MixedBetaVeg)
stepAIC(fullpspmodel, direction = "both", trace=FALSE)
summary(fullpspmodel)

pspmodel<-lme(MaxThaw~CHCA+DRRO+
                  LEDE+VAUL+VAVI+CLAM+
                  SPAN+SPCA+SPRI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data=MixedBetaVeg)
summary(pspmodel)
anova.lme(pspmodel)

hybridpspmodel<-lm(MaxThaw~CAAQ+CHCA+DRRO+LEDE+
                       VAUL+VAVI+CLAM+DIUN+PENE+
                       SPAN+SPCA+SPGI+SPRI+SPWA+
                       CAGY+ERVA+LALA+CAST2+CLCE+
                       CLCR+HYSP+PEAP+PLSC+POST+
                       SPMA, data=MixedBetaVeg)
stepAIC(hybridpspmodel, direction="both", trace=FALSE)
summary(hybridpspmodel)

MixHybridpsp<-lme(MaxThaw ~ CHCA + DRRO + LEDE + VAUL + VAVI + SPAN + 
                      SPCA + SPRI + CLCR , random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedBetaVeg)
summary(MixHybridpsp)
anova.lme(MixHybridpsp)

##Just the PP predictor species
MixedPPVeg<-read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/MixedPPVeg.csv", row.names=1)
View(MixedPPVeg)

PPgrandspmodel<-lm(MaxThaw~BEGL+CACA+CAAQ+
                       CAGY+CHCA+DRRO+EQSC+
                       ERCH+ERVA+LALA+LEDE+
                       LEGR+OXMI+PIMA+RUCH+
                       VAUL+VAVI+AUPA+CAGI+
                       CAST1+CAST2+CLAM+CLCE+
                       CLCR+CLGR+CLMI+DIUN+
                       HYSP+MYAN+PEAP+PENE+
                       PLSC+POST+PTCR+RHPS+
                       SPAN+SPCA+SPFU+SPGI+
                       SPMA+SPRI+SPSQ+SPWA+
                       TONI,data=MixedPPVeg)
stepAIC(PPgrandspmodel, direction = "both", trace=TRUE)
summary(PPgrandspmodel)

PPgrandpspmodel<-lme(MaxThaw ~ BEGL + CACA + CHCA + DRRO + EQSC + LALA + 
                         LEDE + LEGR + RUCH + VAUL + VAVI + CLAM + CLCE + CLGR + HYSP + 
                         PLSC + PTCR + SPCA + SPGI + SPRI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(PPgrandpspmodel)
anova.lme(PPgrandpspmodel)

PPgrandpspmodel2<-lme(MaxThaw ~ BEGL + EQSC + LEDE + RUCH + VAUL +
                          PTCR + SPCA + SPRI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(PPgrandpspmodel2)
anova.lme(PPgrandpspmodel2)

#Just vascular
PPfullvspmodel<-lm(MaxThaw~BEGL+CACA+CAAQ+
                       CAGY+CHCA+DRRO+EQSC+ERCH
                   +ERVA+LALA+LEDE+LEGR+OXMI+
                       PIMA+RUCH+VAUL+VAVI, data=MixedPPVeg)
stepAIC(PPfullvspmodel, direction = "both", trace=FALSE )
summary(PPfullvspmodel)

PPpvspmodel<-lme(MaxThaw ~ BEGL + DRRO + EQSC + ERCH + LEDE + RUCH + 
                     VAVI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(PPpvspmodel)
anova.lme(PPpvspmodel)

PPpvspmodel2<-lme(MaxThaw ~ BEGL + EQSC + ERCH + LEDE 
                  , random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(PPpvspmodel2)
anova.lme(PPpvspmodel2)

#Just cryptogams
PPfullcspmodel<-lm(MaxThaw~AUPA+CAGI+CAST1+
                       CAST2+CLAM+CLCE+CLCR+CLGR
                   +CLMI+DIUN+HYSP+MYAN+PEAP+
                       PENE+PLSC+POST+PTCR+RHPS
                   +SPAN+SPCA+SPFU+SPGI+SPMA+
                       SPRI+SPSQ+SPWA+TONI, data=MixedPPVeg)
stepAIC(PPfullcspmodel, direction = "both", trace=FALSE )
summary(PPfullcspmodel)

PPcspmodel<-lme(MaxThaw ~ POST + SPAN + SPFU + SPMA + SPRI + PEAP, 
                random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(PPcspmodel)
anova.lme(PPcspmodel)

PPcspmodel2<-lme(MaxThaw ~ POST + SPAN + SPFU + SPMA, 
                 random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(PPcspmodel2)
anova.lme(PPcspmodel2)

#hybrid model
pphybridmodel<-lm(MaxThaw~POST + SPAN + SPFU + SPMA + SPRI + PEAP+
                      BEGL + DRRO + EQSC + ERCH + LEDE + RUCH + 
                      VAVI+ CACA + CHCA + LEGR + VAUL+
                      CLAM+CLGR+CLCE+HYSP+PLSC+PTCR+SPCA+SPGI, data=MixedPPVeg)
stepAIC(pphybridmodel, direction = "both", trace=FALSE )
summary(pphybridmodel)

pphybridmix<-lme(MaxThaw ~ SPRI + BEGL + DRRO + EQSC + ERCH + LEDE + 
                     RUCH + VAVI + CACA + CLAM + CLGR + CLCE + HYSP + PLSC + PTCR + 
                     SPCA + SPGI, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(pphybridmix)
anova.lme(pphybridmix)

pphybrid2<-lme(MaxThaw ~ SPRI + BEGL + EQSC + ERCH + LEDE + 
                   RUCH + VAVI + CLGR + PLSC + PTCR + 
                   SPCA, random = list (~1|Year, ~1|Treatment, ~1|PlotName), data = MixedPPVeg)
summary(pphybrid2)
anova.lme(pphybrid2)



### SIMPER
library(vegan)
betacombo.dist<-vegdist(BetaVegCombo)
betaSIMPER<-simper(BetaVegCombo, BetaEnvNew$Thaw.Stage)
summary(betaSIMPER)

simpcombo<-simper(BetaVegCombo, group=BetaEnvNew$Thaw.Stage)
summary(simpcombo)

lapply(simpcombo, FUN=function(x){x$overall})

betawardSIMPER<-simper(BetaVegCombo, BetaEnvNew$WardCluster)
summary(betawardSIMPER)



### Indicators
library(indicspecies)
#individual species
multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage)
indval = multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage, control = how(nperm=999))
summary.multipatt(indval, indvalcomp = TRUE)
indval$sign

indvalori = multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage, duleg = TRUE,             
                      control = how(nperm=999))
summary.multipatt(indvalori, indvalcomp = TRUE)

indvaloriward = multipatt(BetaVegCombo, BetaEnvNew$WardCluster, duleg = TRUE,             
                          control = how(nperm=999))
summary(indvaloriward)

#including species pairs
BetaComb = combinespecies(BetaVegCombo, max.order = 2)$XC
dim(BetaComb)
pairs<-multipatt(BetaComb, BetaEnvNew$Thaw.Stage, duleg = TRUE, 
                 control = how(nperm=999))
summary.multipatt(pairs, indvalcomp = TRUE)

indvalspcomb = multipatt(BetaVegCombo, BetaEnvNew$Thaw.Stage, duleg = TRUE, 
                         control = how(nperm=999))
summary.multipatt(indvalspcomb, indvalcomp = TRUE)

#indicators function for pairs and trios
#stable
scstable= indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=0,           
                     max.order = 3, verbose=TRUE,             
                     At=0.5, Bt=0.2, enableFixed = TRUE)
print(scstable, sqrtIVt = 0.6)
sc2stable=pruneindicators(scstable, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2stable)

#early
scearly= indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=1,           
                    max.order = 3, verbose=TRUE,             
                    At=0.5, Bt=0.2, enableFixed = TRUE)
print(scearly, sqrtIVt = 0.6)
sc2early=pruneindicators(scearly, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2early)

#inter
scint= indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=2,           
                  max.order = 3, verbose=TRUE,             
                  At=0.5, Bt=0.2, enableFixed = TRUE)
print(scint, sqrtIVt = 0.6)
sc2int=pruneindicators(scint, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2int)

#adv
scadv= indicators(X=BetaVegCombo, cluster=BetaEnvNew$Thaw.Stage, group=3,           
                  max.order = 3, verbose=TRUE,             
                  At=0.5, Bt=0.2, enableFixed = TRUE)
print(scadv, sqrtIVt = 0.6)
sc2adv=pruneindicators(scadv, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2adv)


###Interpolation for C calcs


#######################Linear Interpolation of Flux Data#############################################################
#Upload datasets
Interpolation.Matrix.CH4.2017 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix CH4 2017.csv")
Interpolation.Matrix.CH4.2018 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix CH4 2018.csv")
Interpolation.Matrix.CH4.2019 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix CH4 2019.csv")

Interpolation.Matrix.ER.2017 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix ER 2017.csv")
Interpolation.Matrix.ER.2018 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix ER 2018.csv")
Interpolation.Matrix.ER.2019 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix ER 2019.csv")

Interpolation.Matrix.NEE.2017 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix NEE 2017.csv")
Interpolation.Matrix.NEE.2018 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix NEE 2018.csv")
Interpolation.Matrix.NEE.2019 <- read.csv("C:/Users/danci/Desktop/PhD/R directory/BetaVeg/Interpolation Matrix NEE 2019.csv")


#Convert Dataframe to time series before interpolation
Interpolation_MatrixCH42017_TS<-as.ts(Interpolation.Matrix.CH4.2017)
Interpolation_MatrixCH42018_TS<-as.ts(Interpolation.Matrix.CH4.2018)
Interpolation_MatrixCH42019_TS<-as.ts(Interpolation.Matrix.CH4.2019)

Interpolation_MatrixER2017_TS<-as.ts(Interpolation.Matrix.ER.2017)
Interpolation_MatrixER2018_TS<-as.ts(Interpolation.Matrix.ER.2018)
Interpolation_MatrixER2019_TS<-as.ts(Interpolation.Matrix.ER.2019)

Interpolation_MatrixNEE2017_TS<-as.ts(Interpolation.Matrix.NEE.2017)
Interpolation_MatrixNEE2018_TS<-as.ts(Interpolation.Matrix.NEE.2018)
Interpolation_MatrixNEE2019_TS<-as.ts(Interpolation.Matrix.NEE.2019)
#Now interpolate matrix using linear
install.packages("imputeTS")
library("imputeTS")

Interpolation.Matrix.CH4.2017<-na_interpolation(Interpolation_MatrixCH42017_TS, option = "linear")
Interpolation.Matrix.CH4.2018<-na_interpolation(Interpolation_MatrixCH42018_TS, option = "linear")
Interpolation.Matrix.CH4.2019<-na_interpolation(Interpolation_MatrixCH42019_TS, option = "linear")

Interpolation.Matrix.ER.2017<-na_interpolation(Interpolation_MatrixER2017_TS, option = "linear")
Interpolation.Matrix.ER.2018<-na_interpolation(Interpolation_MatrixER2018_TS, option = "linear")
Interpolation.Matrix.ER.2019<-na_interpolation(Interpolation_MatrixER2019_TS, option = "linear")

Interpolation.Matrix.NEE.2017<-na_interpolation(Interpolation_MatrixNEE2017_TS, option = "linear")
Interpolation.Matrix.NEE.2018<-na_interpolation(Interpolation_MatrixNEE2018_TS, option = "linear")
Interpolation.Matrix.NEE.2019<-na_interpolation(Interpolation_MatrixNEE2019_TS, option = "linear")
#Now sum each column to calculate the cumulative CO2 released per unit
Cumulative_CH42017<-colSums(Interpolation.Matrix.CH4.2017)
Cumulative_CH42017<-as.data.frame(Cumulative_CH42017)

Cumulative_CH42018<-colSums(Interpolation.Matrix.CH4.2018)
Cumulative_CH42018<-as.data.frame(Cumulative_CH42018)

Cumulative_CH42019<-colSums(Interpolation.Matrix.CH4.2019)
Cumulative_CH42019<-as.data.frame(Cumulative_CH42019)

Cumulative_ER2017<-colSums(Interpolation.Matrix.ER.2017)
Cumulative_ER2017<-as.data.frame(Cumulative_ER2017)

Cumulative_ER2018<-colSums(Interpolation.Matrix.ER.2018)
Cumulative_ER2018<-as.data.frame(Cumulative_ER2018)

Cumulative_ER2019<-colSums(Interpolation.Matrix.ER.2019)
Cumulative_ER2019<-as.data.frame(Cumulative_ER2019)

Cumulative_NEE2017<-colSums(Interpolation.Matrix.NEE.2017)
Cumulative_NEE2017<-as.data.frame(Cumulative_NEE2017)

Cumulative_NEE2018<-colSums(Interpolation.Matrix.NEE.2018)
Cumulative_NEE2018<-as.data.frame(Cumulative_NEE2018)

Cumulative_NEE2019<-colSums(Interpolation.Matrix.NEE.2019)
Cumulative_NEE2019<-as.data.frame(Cumulative_NEE2019)

write.csv(Cumulative_CH42017, "Cumulative_CH42017.csv")
write.csv(Cumulative_CH42018, "Cumulative_CH42018.csv")
write.csv(Cumulative_CH42019, "Cumulative_CH42019.csv")

write.csv(Cumulative_ER2017, "Cumulative_ER2017.csv")
write.csv(Cumulative_ER2018, "Cumulative_ER2018.csv")
write.csv(Cumulative_ER2019, "Cumulative_ER2019.csv")

write.csv(Cumulative_NEE2017, "Cumulative_NEE2017.csv")
write.csv(Cumulative_NEE2018, "Cumulative_NEE2018.csv")
write.csv(Cumulative_NEE2019, "Cumulative_NEE2019.csv")
#Move data to excel file

library(xlsx)

write.xlsx(Interpolation.Matrix.CH4.2017, "Interpolation.Matrix.CH4.2017.Full.xlsx")
write.xlsx(Interpolation.Matrix.CH4.2018, "Interpolation.Matrix.CH4.2018.Full.xlsx")
write.xlsx(Interpolation.Matrix.CH4.2019, "Interpolation.Matrix.CH4.2019.Full.xlsx")

write.xlsx(Interpolation.Matrix.ER.2017, "Interpolation.Matrix.ER.2017.Full.xlsx")
write.xlsx(Interpolation.Matrix.ER.2018, "Interpolation.Matrix.ER.2018.Full.xlsx")
write.xlsx(Interpolation.Matrix.ER.2019, "Interpolation.Matrix.ER.2019.Full.xlsx")

write.xlsx(Interpolation.Matrix.NEE.2017, "Interpolation.Matrix.NEE.2017.Full.xlsx")
write.xlsx(Interpolation.Matrix.NEE.2018, "Interpolation.Matrix.NEE.2018.Full.xlsx")
write.xlsx(Interpolation.Matrix.NEE.2019, "Interpolation.Matrix.NEE.2019.Full.xlsx")

###Boxplots
#AL depth
library(ggplot2)
library(cowplot)
library(bigleaf)
library(multcompView)
library(plyr)
ALBox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=MaxThaw*-1, fill=Year)) + geom_boxplot()+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept=-46)+
    geom_hline(yintercept=-46+9.01, linetype="dashed")+
    geom_hline(yintercept=-46-9.01, linetype="dashed")+
    geom_hline(yintercept = -120)
ALBox+theme_bw()
#methane
methbox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=JulyPeakMeth.umol.m2.min, fill=Year))+ geom_boxplot()
methbox+theme_bw()
#ER
ERbox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=JulyER, fill=Year)) + geom_boxplot()+ylim(-15,15)+
    geom_hline(yintercept=0, linetype="dashed")
ERbox+theme_bw()
#NEE
NEEbox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=JulyNEE, fill=Year)) + geom_boxplot()+ylim(-15,15)+
    geom_hline(yintercept=0, linetype="dashed")
NEEbox+theme_bw()
##NEE plus meth mass
NEEbyMass.g<-umolCO2.to.gC(BetaEnvNew$JulyNEE)
CH4byMass.CO2Eq<-((BetaEnvNew$JulyPeakMeth.umol.m2.s*(10^-6))*16.04*84)
TotalCO2Equiv<-(NEEbyMass.g+(CH4byMass.CO2Eq))
totalCbox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=TotalCO2Equiv, fill=Year)) + geom_boxplot()+
    geom_hline(yintercept=0, linetype="dashed")
totalCbox+theme_bw()
#ericoids
EriBox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=Ericoids, fill=Year)) + geom_boxplot()+ylim(0,100)
EriBox+theme_bw()
#graminoids
GramBox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=Graminoids, fill=Year)) + geom_boxplot()+ylim(0,100)
GramBox+theme_bw()
#sphagna
SphagBox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=Sphagna, fill=Year)) + geom_boxplot()+ylim(0,100)
SphagBox+theme_bw()
#richness
RichBox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=Richness, fill=Year)) + geom_boxplot()
RichBox+theme_bw()
#multipanel
library(ggpubr)
ggarrange(RichBox+theme_bw(), GramBox+theme_bw(), EriBox+theme_bw(), SphagBox+theme_bw(), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

ggarrange(ERbox+theme_bw(), GPPbox+theme_bw(), NEEbox+theme_bw(),
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
#ClossSeason
ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=Closs.gC.m2, fill=Year)) + geom_boxplot()+
    geom_hline(yintercept=0, linetype="dashed")
#geom_hline(yintercept=106.6, linetype="dashed")+
#geom_hline(yintercept=-27.3, linetype="dashed")

#gpp
GPPbox<-ggplot(BetaEnvNew, aes(x=Thaw.Stage, y=GPP*-1, fill=Year)) + geom_boxplot()+ylim(-15,15)+
    geom_hline(yintercept=0, linetype="dashed")
GPPbox+theme_bw()
###gpp model
gppmod<-lme(GPP~Thaw.Stage*Year, random = ~1|Treatment, data=BetaEnvNew)
summary(gppmod)
anova.lme(gppmod)


#####Interpolation workspace
library(chebpol)

ipol(plot_matrixfiddling)

####Citations
citation(package="vegan")
citation(package="nlme")
citation(package="indicspecies")
citation(package="rcompanion")
citation(package="pvclust")

