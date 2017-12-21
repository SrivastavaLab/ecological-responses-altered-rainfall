source("01_datasets_for_paper1.R")
source("02_functions.R")

#=====combined effects of taxonomy and hydrology?-----------

#piercers
#stats notes: I was careful to scale before rounding such that theta was 0.2-0.5 and all models ha similar theta

bestmodpi<-glm.nb(round(piercer_bio*1)~log(maxvol)+site*prop.overflow.days, data = nococrprdata); Anova(bestmodpi, type=2)
visreg(bestmodpi, "prop.overflow.days", by="site",ylab="Piercer biomass", overlay=FALSE, partial=TRUE, band=FALSE)
m1<-glm.nb(round(piercer_bio*1)~log(maxvol)+prop.overflow.days, data = ardata); Anova(m1, type=2)#ns argentina, but looks like a trend!
m2<-glm.nb(round(piercer_bio*1)~log(maxvol)+prop.overflow.days, data = madata); Anova(m2, type=2)#sig neg macae
m3<-glm.nb(round(piercer_bio*1)~log(maxvol)+prop.overflow.days, data = fgdata); Anova(m3, type=2)#ns FG
visreg(m2, "prop.overflow.days", ylab="Piercer biomass",  partial=TRUE, band=FALSE)

m10<-glm.nb(round(Tabanidae_bio*0.1)~log(maxvol)+site*prop.overflow.days, data = nococrprdata); Anova(m10, type=2)#not sig: note that needed to dial down scalar or else theta way too big
visreg(m10, "prop.overflow.days", by="site",ylab="Tabanidae_bio", overlay=FALSE, partial=TRUE, band=FALSE)
summary(m10)

m1<-glm.nb(round(Tabanidae_bio*1)~log(maxvol)+prop.overflow.days, data = ardata); Anova(m1, type=2)#ns
plot(ardata$Tabanidae_bio~ardata$prop.overflow.days)
m4<-glm.nb(round(Tabanidae_bio*0.1)~log(maxvol)+prop.overflow.days, data = madata); Anova(m4, type=2)#only in 5 brom, ns once scaled appropriately
#Ah, the zeroes, when corrected for log maxvol, drive the model
plot(madata$Tabanidae_bio~madata$prop.overflow.days)
#not enough Tabanid in fg to run model
visreg(m4, "prop.overflow.days", ylab="Tabanid biomass",  partial=TRUE, band=TRUE)

#There are a ton of other piercers in ma and fg...what are they? fg must be insensitive, and ma sensitive
tapply(fulldata$Tabanidae_bio, fulldata$site, datacheck)
#perhaps Tanypodinae in ma? "Empididae_bio" "Tanypodinae_bio" "predCerato_bio" "Ephydridae_bio"  Lampyridae_bio
tapply(fulldata$Tanypodinae_bio, fulldata$site, datacheck)
#Tanypodinae are in Macae (23 brom, mean 0.355) cardoso (30, mean 1.10) and puerto rico (8 brom, mean 0.067)...but CA no hydrology
tapply(fulldata$Empididae_bio, fulldata$site, datacheck)
#Empididae are in Macae (16, mean 0.025) cardoso (23, 0.053)...but CA no hydrology
tapply(fulldata$predCerato_bio, fulldata$site, datacheck)
#predCerato are in Macae (6, mean 0.04) cardoso (22, 0.058),colombia (6), french guiana (28, mean 0.348)...but CA no hydrology
tapply(fulldata$Ephydridae_bio, fulldata$site, datacheck)
#Ephrydidae are entirely in cardoso (29, mean 0.279) and 1 macae brom
tapply(fulldata$Lampyridae_bio, fulldata$site, datacheck)
#a few lmpyridae in macae and cardoso

#Looks like Tanypodinae in Macae is driving it 
m3<-glm.nb(round(Tanypodinae_bio*50)~log(maxvol)+prop.overflow.days, data = madata); Anova(m3, type=2) #ns (p=0.07) negative slope 
summary(m3)
visreg(m3, "prop.overflow.days", ylab="Tanypodinae biomass",  partial=TRUE, band=FALSE)
plot(madata$Tanypodinae_bio~madata$prop.overflow.days)#constrained relationship
plot(fulldata$Tanypodinae_bio~fulldata$prop.overflow.days)
m5<-glm.nb(round(Tanypodinae_bio*50)~log(maxvol)+site*prop.overflow.days, data = maprdata); Anova(m5, type=3) #p=0.05
summary(m5)



#not Empididae in Macae
m3<-glm.nb(round(Empididae_bio*200)~log(maxvol)+prop.overflow.days, data = madata); Anova(m3, type=2) #ns (p=0.65) 
summary(m3)

#and predCerato in FG are not responsive
m3<-glm.nb(round(predCerato_bio*10)~log(maxvol)+prop.overflow.days, data = madata); Anova(m3, type=2) #ns (p=0.20) 
summary(m3)



#though Dsqured suggests that last_wet might be even better than prop.overflow
m5<-glm.nb(round(Tanypodinae_bio*50)~log(maxvol)+site*last_wet, data = maprdata); Anova(m5, type=2) #only sig is last wet p=0.00001
plot(maprdata$Tanypodinae_bio~maprdata$last_wet)#but no pattern...what is going on?!
visreg(m5, "last_wet", ylab="Tanypodinae_bio",  by="site",overlay=FALSE, partial=TRUE, band=FALSE)#nothing there!


bestpi<-glm.nb(round(piercer_bio*10)~log(maxvol)+site+log(mu.scalar), data = nococrdata); Anova(bestpi, type=2)
summary(bestpi) #negative effects of mu (and k if put it in the model too...)
tapply(fulldata$piercer_bio, fulldata$site, datacheck)

#so general and NEGATIVE effects of mu, perhaps due to overflow
bestpi<-glm.nb(round(piercer_bio*10)~log(maxvol)+site+log(mu.scalar), data = nocacocrdata); Anova(bestpi, type=2)
#still there without cardoso
bestpi<-glm.nb(round(piercer_bio*1)~log(maxvol)+site*prop.overflow.days+log(mu.scalar), data = nocacocrdata); Anova(bestpi, type=2)
summary(bestpi)
#so,now site x prop.overflow is sig, and mu has become nonsig...suggesting proximite driver is site-specific effects of mu on prop.overflow

#this was the original model...
bestpi<-glm.nb(round(piercer_bio*10)~log(maxvol)+site+log(mu.scalar)+log(k.scalar), data = nococrprdata); Anova(bestpi, type=2)
bestpi<-glm.nb(round(piercer_bio*1)~log(maxvol)+site*prop.overflow.days+log(mu.scalar)+log(k.scalar), data = nococrprdata); Anova(bestpi, type=2)
#also rainfall become ns after incl site x prop.overflow

aic.lmx(sqrt(nocacocrdata$prop.overflow.days),  gaussian, nocacocrdata)#hmm, but site x k+k2 drives prop.oveflow
m11<-glm((prop.overflow.days)^0.33~log(maxvol)+site*log(mu.scalar), data = nocacocrdata); Anova(m11, type=2)#not mu
par(mfrow=c(2,2)); plot(m11)
