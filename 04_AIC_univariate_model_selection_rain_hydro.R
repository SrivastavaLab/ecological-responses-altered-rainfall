source("01_datasets_for_paper1.R")
source("02_functions.R")

###==lets make this flow

#an AICc approach
#rules: look at all models within delta= 2 of best; if null in the top choose the null
#examine significance of terms within top models, remove models with insignificant terms 
#for very close models, also considered leverage effects of x variable
#otherwise choose top model
#note that final model selection was based on a nd model with fixed theta

#for aic.percent, considered hydro models with both hydrology and temperature
#since temperature never import, final model run for all brom with hydro data (even if they were missing temperature)
#removed effects of maxvol and site to see how much of residual variance could be explained by any model in top lot

#to calculate deviance explained, used best non-null rainfall 9or hydrology) model, even if ns.

#model selection involving hydrology variables needed to omit cardoso (no measurements), and involving temperature also those bromeliads with no ibutton data
#outlier tests suggested removal of a couple of bromeliads for a few variables


response<-as.vector(c("Decomposition", "Nitrogen uptake", "CO2 flux", "Shredder", "Filter feeder", "Scraper", "Gatherer", "Engulfer", "Piercer", "Bacterial density", "Total Invertebrates"))
aic.percent<-as.data.frame(response)

#decomp - rain
aic.lmx(sqrt(fulldata$decomp), gaussian, fulldata) #m2, m0, m9 #choose m0 as simplest
aic.lmx.x(sqrt(fulldata$decomp), gaussian, fulldata)#m17: just site, no vol, no rain
#m2<-glm(sqrt(decomp)~log(maxvol)+site+log(k.scalar), family = gaussian, data = fulldata) confirms k not sig
bestdecomp<-glm(sqrt(decomp)~log(maxvol)+site, family=gaussian, data=fulldata)
fulldata$resid.decomp<-resid(glm(sqrt(decomp)~log(maxvol)+site,data=fulldata, na.action=na.exclude))
bestdecomp.r<-glm(resid.decomp~log(k.scalar), family=gaussian, data = fulldata)
par(mfrow=c(2,2)); plot(bestdecomp.r); Anova(bestdecomp.r, type=2) #k not sig
aic.percent$draintrue[aic.percent$response=="Decomposition"]<-Dsquared(bestdecomp.r, adjust=TRUE)#0
aic.percent$drainfalse[aic.percent$response=="Decomposition"]<-Dsquared(bestdecomp.r, adjust=FALSE)#0
aic.percent$raintype[aic.percent$response=="Decomposition"]<-"ns"

#n15 - eighth root - rain

aic.lmx((no126data$n15.bromeliad.final+4)^0.125, gaussian, no126data)#m8 = site x k +k2
aic.lmx.x((no126data$n15.bromeliad.final+4)^0.125, gaussian, no126data)#same
bestnit<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), family=gaussian, data = no126data)
#k and k2 sig, but not in macae and french guiana
absnit<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site+(log(intended.k)+I(log(intended.k)^2)), family=gaussian, data = no126data)
anova(bestnit,absnit, test="LRT") #Deviance change -0.240, p = 0.01321 *
abs.check(bestnit,absnit)#absolute better !
visreg(absnit, "intended.k", by="site", ylab="Nitrogen uptake", xtrans=log, overlay=TRUE, partial=TRUE, rug=TRUE, band=FALSE)
visreg(absnit, "intended.k", by="site", ylab="Nitrogen uptake",  overlay=TRUE, partial=TRUE, rug=TRUE, band=FALSE)

no126data$tribe<-"Tillansoidae"
no126data$tribe[no126data$site=="argentina"]<-"Bromelioideae"
no126data$tribe[no126data$site=="macae"]<-"Bromelioideae"
no126data$tribe[no126data$site=="cardoso"]<-"Bromelioideae"
no126data$tribe<-as.factor(no126data$tribe)
absnit2<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site+tribe*(log(intended.k)+I(log(intended.k)^2)), family=gaussian, data = no126data)
absnit3<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site+tribe*(log(intended.k)), family=gaussian, data = no126data)
absnit4<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site+tribe+(log(intended.k)+I(log(intended.k)^2)), family=gaussian, data = no126data)
print(aicset<-model.sel(absnit,absnit2, absnit3, absnit4))
visreg(absnit2, "intended.k", by="tribe", ylab="Nitrogen uptake", xtrans=log, overlay=FALSE, partial=TRUE, rug=TRUE, band=FALSE)
no126data$resid.nit<-resid(glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site,data=no126data, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Nitrogen uptake"]<-Dsquared(glm(resid.nit~site*(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, no126data), adjust=TRUE) #0.05186261
aic.percent$drainfalse[aic.percent$response=="Nitrogen uptake"]<-Dsquared(glm(resid.nit~site*(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, no126data), adjust=FALSE) #0.05186261
aic.percent$raintype[aic.percent$response=="Nitrogen uptake"]<-"contingent"



#co2 - rain
aic.lmx(log(nocaprdata$co2.final), gaussian, nocaprdata)#m0, m4, m1
aic.lmx.x(log(nocaprdata$co2.final), gaussian, nocaprdata)#m17..no rain
bestco2<-glm(log(co2.final)~log(maxvol)+site, gaussian, nocaprdata)
nocaprdata$resid.co2<-resid(glm(log(nocaprdata$co2.final)~log(maxvol)+site,data=nocaprdata, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="CO2 flux"]<-Dsquared(glm(resid.co2~(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, nocaprdata), adjust=TRUE) #0
aic.percent$drainfalse[aic.percent$response=="CO2 flux"]<-Dsquared(glm(resid.co2~(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, nocaprdata), adjust=FALSE) #0
aic.percent$raintype[aic.percent$response=="CO2 flux"]<-"ns"

#shredder - rain
aic.lmxnb(round(fulldata$shredder_bio*10), fulldata)#m4 (k+k2) m2
aic.lmxnb.x(round(fulldata$shredder_bio*10), fulldata)#m21 (k+k2) no vol
bestsh<-glm.nb(round(shredder_bio*10)~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), data = fulldata)#yes both k and k2 very sig
aic.lmx(round(fulldata$shredder_bio*10), family=negative.binomial(theta = 1.6521), fulldata)#m4 (k+k2)
fulldata$resid.shr<-resid(glm.nb(round(shredder_bio*10)~log(maxvol)+site, data=fulldata, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Shredder"]<-Dsquared(glm(resid.shr~log(k.scalar)+I(log(k.scalar)^2),family=gaussian, fulldata), adjust=TRUE) #0.0467
aic.percent$drainfalse[aic.percent$response=="Shredder"]<-Dsquared(glm(resid.shr~log(k.scalar)+I(log(k.scalar)^2),family=gaussian, fulldata), adjust=FALSE) #0.0467
aic.percent$raintype[aic.percent$response=="Shredder"]<-"general"

#filter feeder- rain
aic.lmxnb(round(nococr140data$filter.feeder_bio*100), nococr140data)#m7 = site x mu+mu2
aic.lmxnb.x(round(nococr140data$filter.feeder_bio*100), nococr140data)#m7
aic.lmx(round(nococr140data$filter.feeder_bio*100), family=negative.binomial(theta = 1.3407), nococr140data)#m7
bestff<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), data = nococr140data)
betterff<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site*log(mu.scalar)+I(log(mu.scalar)^2), data = nococr140data)
#looks like the effects of mu are entirely driven by cardoso...
bestff2<-glm(round(filter.feeder_bio*100)~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)), family=negative.binomial(theta=bestff$theta),data = nococr140data)
absff<-glm(round(filter.feeder_bio*100)~log(maxvol)+site+(log(intended.mu)+I(log(intended.mu)^2)), family=negative.binomial(theta=bestff$theta),data = nococr140data)
anova(bestff2,absff, test="LRT") #Deviance change -34.241 4.093e-06 ***
abs.check(bestff2,absff)#relative better

nococr140data$resid.ff<-resid(glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site, data=nococr140data, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Filter feeder"]<-Dsquared(glm(resid.ff~site*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, nococr140data), adjust=TRUE) #0.1385931
aic.percent$drainfalse[aic.percent$response=="Filter feeder"]<-Dsquared(glm(resid.ff~site*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, nococr140data), adjust=FALSE) #0.1385931
aic.percent$raintype[aic.percent$response=="Filter feeder"]<-"contingent"

#strength of contingency
Dsquared(glm(resid.ff~site*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, nococr140data), adjust=FALSE) #total rainfall dsq is 22.6%
Dsquared(glm(resid.ff~site+(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, nococr140data),adjust=FALSE) #ignoring interactio gives 10.1%
0.1013384/0.2258979 #so 45% of total rainfall effect is main effect, leaving 55% is contingent

#scraper - rain --------------
aic.lmxnb(round(fulldata$scraper_bio*10), fulldata)#m8 = site x k+k2
aic.lmxnb.x(round(fulldata$scraper_bio*10), fulldata)#m8
aic.lmx(round(fulldata$scraper_bio*10), family=negative.binomial(theta = 1.834), fulldata)#m8
bestsc<-glm.nb(round(scraper_bio*10)~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), data = fulldata)
#scrapers all go down as k does, curve mainly due to PR
bestsc2<-glm(round(scraper_bio*10)~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), family=negative.binomial(theta=bestsc$theta), data = fulldata)
abssc2<-glm(round(scraper_bio*10)~log(maxvol)+site+(log(intended.k)+I(log(intended.k)^2)), family=negative.binomial(theta=bestsc$theta), data = fulldata)
anova(bestsc2,abssc2, test="LRT") #Deviance change -42.037 0.0002562 ***
abs.check(bestsc2,abssc2)#relative better
Dsquared(bestsc, adjust=TRUE)#0.60
fulldata$resid.scraper<-resid(glm.nb(round(scraper_bio*10)~log(maxvol)+site, data=fulldata, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Scraper"]<-Dsquared(glm(resid.scraper~site*(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, fulldata), adjust=TRUE) #0.09918389
aic.percent$drainfalse[aic.percent$response=="Scraper"]<-Dsquared(glm(resid.scraper~site*(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, fulldata), adjust=FALSE) #
aic.percent$raintype[aic.percent$response=="Scraper"]<-"contingent"

#gatherer - rain
aic.lmxnb(round(no67185data$gatherer_bio*10), no67185data)#no convergence with x 100, but with x 10 get m8 = site x k+k2
aic.lmxnb.x(round(no67185data$gatherer_bio*10), no67185data)#m8
aic.lmx(round(no67185data$gatherer_bio*10),family=negative.binomial(theta = 1.2911), no67185data)# m8
bestga<-glm.nb(round(gatherer_bio*10)~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), data = no67185data)
bestga2<-glm(round(gatherer_bio*10)~log(maxvol)+site*(log(k.scalar)+I(log(k.scalar)^2)), family=negative.binomial(theta=bestga$theta), data = no67185data)
absga<-glm(round(gatherer_bio*10)~log(maxvol)+site+(log(intended.k)+I(log(intended.k)^2)), family=negative.binomial(theta=bestga$theta), data = no67185data)
anova(bestga2,absga, test="LRT") #-12  -41.918 0.001899 **
abs.check(bestga2,absga)#relative better
no67185data$resid.ga<-resid(glm.nb(round(gatherer_bio*10)~log(maxvol)+site, data=no67185data, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Gatherer"]<-Dsquared(glm(resid.ga~site*(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, no67185data), adjust=TRUE) #0.06993203
aic.percent$drainfalse[aic.percent$response=="Gatherer"]<-Dsquared(glm(resid.ga~site*(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, no67185data), adjust=FALSE) #0.06993203
aic.percent$raintype[aic.percent$response=="Gatherer"]<-"contingent"

#engulfer - rain
aic.lmxnb(round(noargco13data$engulfer_bio*100), noargco13data)# m0 m3 m1...NO RAINFALL EFFECT
aic.lmxnb.x(round(noargco13data$engulfer_bio*100), noargco13data) #m17..no rain effect, no vol effect
aic.lmx(round(noargco13data$engulfer_bio*100), family=negative.binomial(theta = 0.3697), noargco13data)#m0 m3 m1
noargco13data$resid.en<-resid(glm.nb(round(noargco13data$engulfer_bio*100)~log(maxvol)+site, data=noargco13data, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Engulfer"]<-Dsquared(glm(resid.en~(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, noargco13data), adjust=TRUE) #0
aic.percent$drainfalse[aic.percent$response=="Engulfer"]<-Dsquared(glm(resid.en~(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, noargco13data), adjust=FALSE) 
aic.percent$raintype[aic.percent$response=="Engulfer"]<-"ns"

#piercer - rain -----------
aic.lmxnb(round(nococrdata$piercer_bio*10), nococrdata) #m1,m9, m2
aic.lmxnb.x(round(nococrdata$piercer_bio*10), nococrdata)#m18: mu, no vol
aic.lmx(round(nococrdata$piercer_bio*10), family=negative.binomial(theta = 0.364), nococrdata)#m1,m9,m2
#I decided to include puerto rico with 8 brom, make sure that this is in concordance too (note: results changed after corrected piercer definition)
bestpi<-glm.nb(round(piercer_bio*10)~log(maxvol)+site+log(mu.scalar), data = nococrdata); Anova(bestpi, type=2)
nextpi<-glm.nb(round(piercer_bio*10)~log(maxvol)+site+log(mu.scalar)*log(k.scalar), data = nococrdata); Anova(nextpi, type=2)
nococrdata$resid.pi<-resid(glm.nb(round(nococrdata$piercer_bio*10)~log(maxvol)+site, data=nococrdata, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Piercer"]<-Dsquared(glm(resid.pi~log(mu.scalar),family=gaussian, nococrdata), adjust=TRUE) 
aic.percent$drainfalse[aic.percent$response=="Piercer"]<-Dsquared(glm(resid.pi~log(mu.scalar),family=gaussian, nococrdata), adjust=FALSE) 
aic.percent$raintype[aic.percent$response=="Piercer"]<-"general"

aic.lmxnb(round(nococrprdata$piercer_bio*10), nococrprdata)

#bacteria - rain ---------------
aic.lmxnb(round(noargco123data$bacteria.per.nl.final*100), noargco123data) #m1 m5 ...mu or site x mu
aic.lmxnb.x(round(noargco123data$bacteria.per.nl.final*100), noargco123data)#m1
aic.lmx(round(noargco123data$bacteria.per.nl.final*100), family=negative.binomial(theta = 2.6407), noargco123data) #m1 m5 ...mu or site x mu
m5<-glm.nb(round(bacteria.per.nl.final*100)~log(maxvol)+site*log(mu.scalar),  data = noargco123data)#intercation ns at p=0.097
m1<-glm.nb(round(bacteria.per.nl.final*100)~log(maxvol)+site+log(mu.scalar),  data = noargco123data)
anova(m1,m5, test="Chisq")#not sig different fit, and site:log mu not sig in m5, so m1 best
bestba<-glm.nb(round(bacteria.per.nl.final*100)~log(maxvol)+site+log(mu.scalar),  data = noargco123data)
visreg(bestba, "mu.scalar", by = "site")
noargco123data$resid.ba<-resid(glm.nb(round(bacteria.per.nl.final*100)~log(maxvol)+site, data=noargco123data, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Bacterial density"]<-Dsquared(glm(resid.ba~log(mu.scalar),family=gaussian, noargco123data), adjust=TRUE) #0.0346
aic.percent$drainfalse[aic.percent$response=="Bacterial density"]<-Dsquared(glm(resid.ba~log(mu.scalar),family=gaussian, noargco123data), adjust=FALSE)
aic.percent$raintype[aic.percent$response=="Bacterial density"]<-"general"

#total bio - rain ------------------
aic.lmxnb(round(fulldata$totalbio*10), fulldata)#m0 m1...NO RAINFALL EFFECT 
aic.lmxnb.x(round(fulldata$totalbio*10), fulldata)#m0
aic.lmx(round(fulldata$totalbio*10), family=negative.binomial(theta = 2.5231), fulldata) #mo m1
besttb<-glm.nb(round(totalbio*10)~log(maxvol)+site,  data = fulldata)
fulldata$resid.tb<-resid(glm.nb(round(totalbio*10)~log(maxvol)+site, data=fulldata, na.action=na.exclude))
aic.percent$draintrue[aic.percent$response=="Total Invertebrates"]<-Dsquared(glm(resid.tb~log(mu.scalar),family=gaussian, fulldata), adjust=TRUE) #0
aic.percent$drainfalse[aic.percent$response=="Total Invertebrates"]<-Dsquared(glm(resid.tb~log(mu.scalar),family=gaussian, fulldata), adjust=FALSE) #0
aic.percent$raintype[aic.percent$response=="Total Invertebrates"]<-"ns"

#prey_biomass
aic.lmxnb(round(fulldata$prey_biomass*10), fulldata)#site *k+k2
aic.lmxnb.x(round(fulldata$prey_biomass*10), fulldata)#same: m8
aic.lmx(round(fulldata$prey_biomass*10), family=negative.binomial(theta = 2.77), fulldata) #m8 #m4

#predator_biomass
aic.lmxnb(round(fulldata$predator_biomass*10), fulldata)#m0 m1...NO RAINFALL EFFECT 
aic.lmxnb.x(round(fulldata$predator_biomass*10), fulldata)#now m19 (site = k, no vol)
aic.lmx(round(fulldata$predtor_biomass*10), family=negative.binomial(theta = 0.40), fulldata) #mo m1


#plots ---------------
visreg(bestnit, "k.scalar", by="site", ylab="Nitrogen uptake", overlay=TRUE, partial=FALSE, band=FALSE)
visreg(bestsh, "k.scalar", by="site", ylab="Shredder biomass", overlay=TRUE, partial=FALSE, band=FALSE)
visreg(bestff, "mu.scalar", by="site", ylab="Filter feeder biomass",overlay=TRUE, partial=TRUE, band=FALSE)
visreg(bestsc, "k.scalar", by="site", ylab="Scraper biomass",overlay=TRUE, partial=FALSE, band=FALSE)
visreg(bestga, "k.scalar", by="site", ylab="Gatherer biomass",overlay=TRUE, partial=FALSE, band=FALSE)
visreg(bestpi, "mu.scalar", by="site", ylab="Piercer biomass",overlay=TRUE, partial=FALSE, band=FALSE)
visreg(bestba, "mu.scalar", by="site", ylab="Bacterial biomass", overlay=TRUE, partial=TRUE, band=FALSE)

#+ hydrology
#===best Hydrology models-------------------------

#here I made the decision to filter out argentina_15 because it's cv.depth value is so much greater than that of other bromeliads
#this bromeliad also has 100% proportion days dry, also the dataset max, and it is almost dataset max for last_wet
#but it is average for long_dry, a bit oddly...because this is based on #records, not #days.

#decomp- hydro
nocadatatemp<-filter(nocadata, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)%>%filter(site_brom.id%nin%"argentina_15")
nocadata.noout<-filter(nocadata, site_brom.id%nin%"argentina_15")
aic.hydro(sqrt(nocadatatemp$decomp), sqrt(nocadatatemp$decomp)~log(nocadatatemp$maxvol)+nocadatatemp$site, gaussian, nocadatatemp)
#m1 m3
aic.hydro.pure(sqrt(nocadata.noout$decomp), gaussian, nocadata.noout)#m1...cv.depth now best!
bestmoddecomp<-glm(sqrt(decomp)~log(maxvol)+site+cv.depth, family=gaussian, data = nocadata.noout)#m3
Anova(bestmoddecomp, type=2)#very sig
nextmoddecomp<-glm(sqrt(decomp)~log(maxvol)+site+prop.driedout.days, family=gaussian, data = nocadata.noout); Anova(nextmoddecomp, type=2)
#both m1 and m3 are very good fits
par(mfrow=c(1,1))
visreg(bestmoddecomp, "cv.depth", by="site", ylab="Decomposition", overlay=TRUE, partial=TRUE, band=FALSE)
visreg(nextmoddecomp, "prop.driedout.days", by="site", ylab="Decomposition", overlay=TRUE, partial=TRUE, band=FALSE)
#I think m3 has better coverage of data over full range of x variable, better supported
nocadata.noout$resid.decomp.hydro<-resid(glm(sqrt(decomp)~log(maxvol)+site,data=nocadata.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Decomposition"]<-Dsquared(glm(resid.decomp.hydro~cv.depth,family=gaussian, nocadata.noout), adjust=TRUE) #0.02556
aic.percent$dhydrofalse[aic.percent$response=="Decomposition"]<-Dsquared(glm(resid.decomp.hydro~cv.depth,family=gaussian, nocadata.noout), adjust=FALSE)
aic.percent$hydrotype[aic.percent$response=="Decomposition"]<-"general"
aic.percent$bestmodel[aic.percent$response=="Decomposition"]<-"hydrology"

#nitrogen - hydro---------------------------
no126datatemp<-filter(no126data, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)%>%filter(site_brom.id%nin%"argentina_15")
no126data.noout<-filter(no126data, site_brom.id%nin%"argentina_15")%>%filter(site%nin%"cardoso")
aic.hydro((no126datatemp$n15.bromeliad.final+4)^0.125, (no126datatemp$n15.bromeliad.final+4)^0.125~log(no126datatemp$maxvol)+no126datatemp$site*(log(no126datatemp$k.scalar)+I(log(no126datatemp$k.scalar)^2)), gaussian, no126datatemp)
#best is rain, then chng_mean_temp but this is not sig factor as seen here:
onemodN2<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site+change_mean_temp, family=gaussian, data = no126datatemp); Anova(onemodN2, type=2)#not sig
#therefore, because of misleading effect of temp, used this model for rain vs hydro:
aic.hydro.pure((no126data.noout$n15.bromeliad.final+4)^0.125, gaussian, no126data.noout) #best is m21 long_dry m19 site*prop.dry.dsys as temp not in candidate set
bestmodN2<-glm((n15.bromeliad.final+4)^0.125~log(maxvol)+site*prop.driedout.days, family=gaussian, data = no126data.noout); Anova(bestmodN2, type=2)
#sig site x long dry interaction
visreg(bestmodN2, "long_dry", by="site")#driven mainly by neg effect in argentina and macae vs tiny positive effect pr; other sites have hardly any variance in longest dry
visreg(bestmodN2, "long_dry", by="site",ylab="Nitrogen uptake", overlay=FALSE, partial=TRUE, band=FALSE)
Anova(bestmodN2, type=2)#most sites increase N uptake with long dry, only FG and one other goes down
#the fit of these models is very similar, but the mean.depth model has much better coverage in x axis so is more convincing
#note that CR, CO and MA decrease with increasing mean depth, AR and PR opposite pattern...no link to bromeliad tribe!
no126data.noout$resid.nit.hydro<-resid(glm((no126data.noout$n15.bromeliad.final+4)^0.125~log(no126data.noout$maxvol)+no126data.noout$site, family=gaussian, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Nitrogen uptake"]<-Dsquared(glm(resid.nit.hydro~site*long_dry,family=gaussian, no126data.noout), adjust=TRUE) #0.0269
aic.percent$dhydrofalse[aic.percent$response=="Nitrogen uptake"]<-Dsquared(glm(resid.nit.hydro~site*long_dry,family=gaussian, no126data.noout), adjust=FALSE)
aic.percent$hydrotype[aic.percent$response=="Nitrogen uptake"]<-"contingent"
aic.percent$bestmodel[aic.percent$response=="Nitrogen uptake"]<-"rain"

#co2 - hydro-----------------------
nocaprdatatemp<-filter(nocaprdata, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)
nocaprdata.noout<-filter(nocaprdata, site_brom.id%nin%"argentina_15")
aic.hydro(log(nocaprdatatemp$co2.final), log(nocaprdatatemp$co2.final)~log(nocaprdatatemp$maxvol)+nocaprdatatemp$site, gaussian, nocaprdatatemp)
aic.hydro.pure(nocaprdata.noout$co2.final, gaussian, nocaprdata.noout)#null, followed by mean.depth
#m0 m6 m7 can confirm that neither last_wet or change cv temp are sig variables in these models..so rainfall model is best (which was NULL)
bestmodco2<-glm(log(co2.final)~log(maxvol)+site, family=gaussian, data = nocaprdata.noout); Anova(bestmodco2, type=2)
nocaprdata.noout$resid.co2.hydro<-resid(glm(log(co2.final)~log(maxvol)+site,data=nocaprdata.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="CO2 flux"]<-Dsquared(glm(resid.co2.hydro~mean.depth,family=gaussian, nocaprdata.noout), adjust=TRUE) #0.0127
aic.percent$dhydrofalse[aic.percent$response=="CO2 flux"]<-Dsquared(glm(resid.co2.hydro~mean.depth,family=gaussian, nocaprdata.noout), adjust=FALSE) #0.0127
aic.percent$hydrotype[aic.percent$response=="CO2 flux"]<-"ns"
aic.percent$bestmodel[aic.percent$response=="CO2 flux"]<-"neither"

#filter feeder - hydro ---------------------------
nococr140datatemp<-filter(nococr140data, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)%>%filter(site_brom.id%nin%"argentina_15")
nococr140data.noout<-filter(nococr140data, site_brom.id%nin%"argentina_15")

aic.hydro.nb(round(nococr140datatemp$filter.feeder_bio*100),round(nococr140datatemp$filter.feeder_bio*100)~log(nococr140datatemp$maxvol)+nococr140datatemp$site*(log(nococr140datatemp$mu.scalar)+I(log(nococr140datatemp$mu.scalar)^2)),nococr140datatemp)
#best model is m17: cv.depth*site
aic.hydro(round(nococr140datatemp$filter.feeder_bio*100),round(nococr140datatemp$filter.feeder_bio*100)~log(nococr140datatemp$maxvol)+nococr140datatemp$site*(log(nococr140datatemp$mu.scalar)+I(log(nococr140datatemp$mu.scalar)^2)),family=negative.binomial(theta = 0.856),nococr140datatemp)
#again m17 is best
aic.hydro.nb.best("filter.feeder_bio", nococr140data.noout,100) #pure hydrology: cv.depth x site
aic.hydro.nb.purex(round(nococr140data.noout$filter.feeder_bio*100), nococr140data.noout) #same, without max vol 
bestmodff<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site*cv.depth,  data = nococr140data.noout)
par(mfrow=c(1,1)); visreg(bestmodff, "cv.depth", by="site", ylab="Filter feeder biomass", overlay=TRUE, partial=TRUE, band=FALSE)

#interestingly, all three sites go down roughly the same rate, but macae is more resistent, there are still 2 bromeliads driving this after leaky issues fixed!
#macae_B24, macae_B11 driving this pattern; they both seemed legit looking at depth data...
Anova(bestmodff, type=3) #interaction sig
tapply(fulldata$cv.depth, fulldata$site, max)
nococr140data.noout$resid.ff.hydro<-resid(glm.nb(round(nococr140data.noout$filter.feeder_bio*100)~log(maxvol)+site,data=nococr140data.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Filter feeder"]<-Dsquared(glm(resid.ff.hydro~site*cv.depth,family=gaussian, nococr140data.noout), adjust=TRUE) #0.1954
aic.percent$dhydrofalse[aic.percent$response=="Filter feeder"]<-Dsquared(glm(resid.ff.hydro~site*cv.depth,family=gaussian, nococr140data.noout), adjust=FALSE) #0.1954
aic.percent$hydrotype[aic.percent$response=="Filter feeder"]<-"contingent"
aic.percent$bestmodel[aic.percent$response=="Filter feeder"]<-"hydrology"


m1<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site*cv.depth+site*(log(mu.scalar)+I(log(mu.scalar)^2)),  data = nococr140data.noout)
m1a<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site*cv.depth+(log(mu.scalar)+I(log(mu.scalar)^2)),  data = nococr140data.noout)
anova(m1, m1a)# effect of rain x site interaction: df = 6 , likelihood ratio = 10.41363  = 0.1082794

m2<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)),  data = nococr140data.noout)
m2a<-glm.nb(round(filter.feeder_bio*100)~log(maxvol)+site+(log(mu.scalar)+I(log(mu.scalar)^2)),  data = nococr140data.noout)
anova(m2, m2a) #effect of rain x site interaction: df = 8, likehood ratio = 34.9, p =0.00003 

#strength of contingency
Dsquared(glm(resid.ff.hydro~site*cv.depth,family=gaussian, nococr140data.noout), adjust=FALSE) #total rainfall dsq is 31.6%
Dsquared(glm(resid.ff.hydro~site+cv.depth,family=gaussian, nococr140data.noout), adjust=FALSE) #ignoring interactio gives 26.7%
0.2671396/0.316197 #so 85% of total rainfall effect is main effect, only 15% is contingent

#explored extensively other version of filter feeder hydro model, site x cv.depth still best (non-convergent models removed)


#shredder - hydro -------------------
aic.hydro.nb(round(nocadatatemp$shredder_bio*10), round(nocadatatemp$shredder_bio*10)~log(nocadatatemp$maxvol)+nocadatatemp$site+log(nocadatatemp$k.scalar)+I(log(nocadatatemp$k.scalar)^2), nocadatatemp)
aic.hydro(round(nocadatatemp$shredder_bio*10), round(nocadatatemp$shredder_bio*10)~log(nocadatatemp$maxvol)+nocadatatemp$site+log(nocadatatemp$k.scalar)+I(log(nocadatatemp$k.scalar)^2), family=negative.binomial(theta = 1.2108), nocadatatemp)
aic.hydro(round(nocadatatemp$shredder_bio*10), round(nocadatatemp$shredder_bio*10)~log(nocadatatemp$maxvol)+nocadatatemp$site+log(nocadatatemp$k.scalar)+I(log(nocadatatemp$k.scalar)^2), family=negative.binomial(theta = 1.27), nocadatatemp)
#because temperature was subsequently found to be nonsig, we switched to a no temp series of models for selection:
aic.hydro.nb.notemp(round(nocadata$shredder_bio*10), round(nocadata$shredder_bio*10)~log(nocadata$maxvol)+nocadata$site+log(nocadata$k.scalar)+I(log(nocadata$k.scalar)^2), nocadata)
aic.hydro.notemp(round(nocadata$shredder_bio*10), round(nocadata$shredder_bio*10)~log(nocadata$maxvol)+nocadata$site+log(nocadata$k.scalar)+I(log(nocadata$k.scalar)^2), family=negative.binomial(theta = 1.22), nocadata)

aic.hydro.nb.best("shredder_bio", nocadata.noout,10) # site + proportion overflow best now
aic.hydro.nb.purex(round(nocadata.noout$shredder_bio*10), nocadata.noout)#same, without mxvol
#original rainfall model (m0) best
bestmodsh<-glm.nb(round(shredder_bio*10)~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), data = nocadata.noout)
besthydrosh<-glm.nb(round(shredder_bio*10)~log(maxvol)+site+prop.overflow.days, data = nocadata.noout); Anova(besthydrosh, type=2)

nocadata.noout$resid.ff.hydro<-resid(glm.nb(round(shredder_bio*10)~log(maxvol)+site,data=nocadata.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Shredder"]<-Dsquared(glm(resid.ff.hydro~prop.overflow.days,family=gaussian, nocadata.noout), adjust=TRUE) #0
aic.percent$dhydrofalse[aic.percent$response=="Shredder"]<-Dsquared(glm(resid.ff.hydro~prop.overflow.days,family=gaussian, nocadata.noout), adjust=FALSE) #0
aic.percent$hydrotype[aic.percent$response=="Shredder"]<-"general"
aic.percent$bestmodel[aic.percent$response=="Shredder"]<-"rain"


#scraper - hydro -------------------

aic.hydro.nb(round(nocadatatemp$scraper_bio*10), round(nocadatatemp$scraper_bio*10)~log(nocadatatemp$maxvol)+nocadatatemp$site*(log(nocadatatemp$k.scalar)+I(log(nocadatatemp$k.scalar)^2)), nocadatatemp)
aic.hydro(round(nocadatatemp$scraper_bio*10), round(nocadatatemp$scraper_bio*10)~log(nocadatatemp$maxvol)+nocadatatemp$site*(log(nocadatatemp$k.scalar)+I(log(nocadatatemp$k.scalar)^2)), family=negative.binomial(theta = 1.338),nocadatatemp)
# m19 (site*prop.driedout.days) m5  (long_dry)
aic.hydro.nb.best("scraper_bio", nocadata.noout,10)#site x propdriedoutdays
aic.hydro.nb.purex(round(nocadata.noout$scraper_bio*10), nocadata.noout)#same

bestmodsc<-glm.nb(round(scraper_bio*10)~log(maxvol)+site*prop.driedout.days,  data = nocadata.noout); Anova(bestmodsc, type=2)
visreg(bestmodsc, "prop.driedout.days", by="site", ylab="Scraper biomass", overlay=TRUE, partial=TRUE, band=FALSE)
nextmodsc<-glm.nb(round(scraper_bio*10)~log(maxvol)+site+long_dry,  data = nocadata.noout); Anova(nextmodsc, type=2)
visreg(nextmodsc, "long_dry", by="site", ylab="Scraper biomass", overlay=TRUE, partial=TRUE, band=FALSE)

#Actually the coverage of data much better for m19, and I notice Dsquared higher too
nocadata.noout$resid.sc.hydro<-resid(glm.nb(round(scraper_bio*10)~log(maxvol)+site,data=nocadata.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Scraper"]<-Dsquared(glm(resid.sc.hydro~site*prop.driedout.days,family=gaussian, nocadata.noout), adjust=TRUE) #was 0.0765, now 0.059...this is effect of correcting leaky?
aic.percent$dhydrofalse[aic.percent$response=="Scraper"]<-Dsquared(glm(resid.sc.hydro~site*prop.driedout.days,family=gaussian, nocadata.noout), adjust=FALSE)
aic.percent$hydrotype[aic.percent$response=="Scraper"]<-"contingent"
aic.percent$bestmodel[aic.percent$response=="Scraper"]<-"hydrology"


#gatherer - hydro -------------------

no67185datatemp<-filter(no67185data, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)
no67185data.noout<-filter(no67185data, site_brom.id%nin%"argentina_15")

aic.hydro.nb(round(no67185datatemp$gatherer_bio*10), round(no67185datatemp$gatherer_bio*10)~log(no67185datatemp$maxvol)+no67185datatemp$site*(I(log(no67185datatemp$k.scalar)^2)), no67185datatemp)
aic.hydro(round(no67185datatemp$gatherer_bio*10), round(no67185datatemp$gatherer_bio*10)~log(no67185datatemp$maxvol)+no67185datatemp$site*(I(log(no67185datatemp$k.scalar)^2)), family=negative.binomial(theta = 1.0889),no67185datatemp)
aic.hydro.nb.best("gatherer_bio", no67185data.noout,10)#last wet
aic.hydro.nb.purex(round(no67185data.noout$gatherer_bio*10), no67185data.noout)#same
#m0 only model...rainfall best!
bestmodga<-glm.nb(round(gatherer_bio*10)~log(maxvol)+site*(I(log(k.scalar)^2)), data = no67185data)
hydroga<-glm.nb(round(gatherer_bio*10)~log(maxvol)+site+cv.depth, data = no67185data.noout); Anova(hydroga, type=2)
hydroga2<-glm.nb(round(gatherer_bio*10)~log(maxvol)+site+last_wet, data = no67185data.noout); Anova(hydroga2, type=2)
#cv.depth general model is sig, just not as good as rainfall...
#aic.hydro.nb(round(no67185leakydatatemp$gatherer_bio*10), round(no67185leakydatatemp$gatherer_bio*10)~log(no67185leakydatatemp$maxvol)+no67185leakydatatemp$site*(I(log(no67185leakydatatemp$k.scalar)^2)), no67185leakydatatemp)
#even with leaky brom removed, best is not hydro
no67185data.noout$resid.ga.hydro<-resid(glm.nb(round(gatherer_bio*10)~log(maxvol)+site,data=no67185data.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Gatherer"]<-Dsquared(glm(resid.ga.hydro~site+last_wet,family=gaussian, no67185data.noout), adjust=TRUE)
aic.percent$dhydrofalse[aic.percent$response=="Gatherer"]<-Dsquared(glm(resid.ga.hydro~site+last_wet,family=gaussian, no67185data.noout), adjust=FALSE)
aic.percent$hydrotype[aic.percent$response=="Gatherer"]<-"general"
aic.percent$bestmodel[aic.percent$response=="Gatherer"]<-"rain"

#engulfer - hydro ------------------
noargco13datatemp<-filter(noargco13data, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)
noargco13data.noout<-filter(noargco13data, site_brom.id%nin%"argentina_15")
aic.hydro.nb(round(noargco13datatemp$engulfer_bio*100), round(noargco13datatemp$engulfer_bio*100)~log(noargco13datatemp$maxvol)+noargco13datatemp$site, noargco13datatemp)

aic.hydro.nb.best("engulfer_bio", noargco13data.noout,100)#m4
aic.hydro.nb.purex(round(noargco13data.noout$engulfer_bio*10), noargco13data.noout)#similar positive effect of mean depth when volume omitted
#m0 m4 (mean depth) m8 m24 m3
#so stick with best rainfall, which is null
bestmoden<-glm.nb(round(engulfer_bio*100)~log(maxvol)+site,  data = noargco13data.noout)
hydroen<-glm.nb(round(engulfer_bio*100)~log(maxvol)+site+mean.depth,  data = noargco13data.noout); Anova(hydroen, type=2)
#best hydro model is not sig
noargco13data.noout$resid.en.hydro<-resid(glm.nb(round(noargco13data.noout$engulfer_bio*100)~log(maxvol)+site,data=noargco13data.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Engulfer"]<-Dsquared(glm(resid.en.hydro~mean.depth,family=gaussian, noargco13data.noout), adjust=TRUE) #0
aic.percent$dhydrofalse[aic.percent$response=="Engulfer"]<-Dsquared(glm(resid.en.hydro~mean.depth,family=gaussian, noargco13data.noout), adjust=FALSE) #0
aic.percent$hydrotype[aic.percent$response=="Engulfer"]<-"ns"
aic.percent$bestmodel[aic.percent$response=="Engulfer"]<-"neither"

#piercer-hydro -----------------
nococrprdatatemp<-filter(nococrprdata, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)
nococrprdata.noout<-filter(nococrprdata, site_brom.id%nin%"argentina_15")
#nococrprleakydatatemp<-filter(nococrprdata, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
# filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
#filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)%>%
#filter(site_brom.id%nin%c("macae_B24", "macae_B22", "macae_B9", "macae_B2", "macae_B11", "macae_B41", "argentina_15"))

aic.hydro.nb(round(nococrprdatatemp$piercer_bio*10), round(nococrprdatatemp$piercer_bio*10)~log(nococrprdatatemp$maxvol)+nococrprdatatemp$site*(log(nococrprdatatemp$k.scalar)+I(log(k.scalar)^2)), nococrprdatatemp)
aic.hydro(round(nococrprdatatemp$piercer_bio*10), round(nococrprdatatemp$piercer_bio*10)~log(nococrprdatatemp$maxvol)+nococrprdatatemp$site*(log(nococrprdatatemp$k.scalar)+I(log(k.scalar)^2)), family=negative.binomial(theta = 0.3292),nococrprdatatemp)

aic.hydro.nb.best("piercer_bio", nococrprdata.noout,10)#overflow x site
aic.hydro.nb.purex(round(nococrprdata.noout$piercer_bio*10), nococrprdata.noout)
#m18 note only 3 sites
#aic.hydro.nb(round(nococrprleakydatatemp$piercer_bio*10), round(nococrprleakydatatemp$piercer_bio*10)~log(nococrprleakydatatemp$maxvol)+nococrprleakydatatemp$site*(log(nococrprleakydatatemp$k.scalar)), nococrprleakydatatemp)
#still m18 m6
bestmodpi<-glm.nb(round(piercer_bio*10)~log(maxvol)+site*prop.overflow.days, data = nococrprdata.noout); Anova(bestmodpi, type=3)
visreg(bestmodpi, "prop.overflow.days", by="site",ylab="Piercer biomass", overlay=TRUE, partial=TRUE, band=FALSE)
nococrprdata.noout$resid.pi.hydro<-resid(glm.nb(round(piercer_bio*10)~log(maxvol)+site,data=nococrprdata.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Piercer"]<-Dsquared(glm(resid.pi.hydro~site*prop.overflow.days,family=gaussian, nococrprdata.noout), adjust=TRUE) #0.0742
aic.percent$dhydrofalse[aic.percent$response=="Piercer"]<-Dsquared(glm(resid.pi.hydro~site*prop.overflow.days,family=gaussian, nococrprdata.noout), adjust=FALSE)
aic.percent$hydrotype[aic.percent$response=="Piercer"]<-"contingent"
aic.percent$bestmodel[aic.percent$response=="Piercer"]<-"hydrology"


#bacteria - hydro ---------------
noargco123datatemp<-filter(noargco123data, mean_temp%nin%NA)%>%filter(cv_mean_temp%nin%NA)%>%
  filter(cv.depth%nin%NA)%>%filter(long_dry%nin%NA)%>%filter(last_wet%nin%NA)%>%
  filter(prop.overflow.days%nin%NA)%>%filter(prop.driedout.days%nin%NA)
noargco123data.noout<-filter(noargco123data, site_brom.id%nin%"argentina_15")
aic.hydro.nb(round(noargco123datatemp$bacteria.per.nl.final*100), round(noargco123datatemp$bacteria.per.nl.final*100)~log(noargco123datatemp$maxvol)+noargco123datatemp$site+(log(noargco123datatemp$mu.scalar)), noargco123datatemp)
aic.hydro(round(noargco123datatemp$bacteria.per.nl.final*100), round(noargco123datatemp$bacteria.per.nl.final*100)~log(noargco123datatemp$maxvol)+noargco123datatemp$site+(log(noargco123datatemp$mu.scalar)), family=negative.binomial(theta = 2.5895),noargco123datatemp)

aic.hydro.nb.best("bacteria.per.nl.final", noargco123data.noout,100)#last wet
#m0 m6 (last_wet)
bestmodba<-glm.nb(round(bacteria.per.nl.final*100)~log(maxvol)+site+log(mu.scalar),  data = noargco123data.noout)
hydroba<-glm.nb(round(bacteria.per.nl.final*100)~log(maxvol)+site+last_wet,  data = noargco123data.noout); Anova(hydroba, type=2)
# so best hydro model is not sig
#therefore, use null model for best hydro when competing with best rain, refer to bacteria rain aic tables
noargco123data.noout$resid.ba.hydro<-resid(glm.nb(round(piercer_bio*10)~log(maxvol)+site,data=noargco123data.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Bacterial density"]<-Dsquared(glm(resid.ba.hydro~last_wet,family=gaussian, noargco123data.noout), adjust=TRUE) #0.0742
aic.percent$dhydrofalse[aic.percent$response=="Bacterial density"]<-Dsquared(glm(resid.ba.hydro~last_wet,family=gaussian, noargco123data.noout), adjust=FALSE) #0.0742
aic.percent$hydrotype[aic.percent$response=="Bacterial density"]<-"ns"
aic.percent$bestmodel[aic.percent$response=="Bacterial density"]<-"rain"

#totalbio - hydro
fulldatatemp<-temphydrotruedata
fulldatatemp.noout<-filter(fulldatatemp, site_brom.id%nin%"argentina_15")
fulldata.noout<-filter(fulldata, site_brom.id%nin%"argentina_15")
aic.hydro.nb(round(fulldatatemp.noout$totalbio*10), round(fulldatatemp.noout$totalbio*10)~log(fulldatatemp.noout$maxvol)+fulldatatemp.noout$site, fulldatatemp.noout)
aic.hydro(round(fulldatatemp.noout$totalbio*10), round(fulldatatemp.noout$totalbio*10)~log(fulldatatemp.noout$maxvol)+fulldatatemp.noout$site, family=negative.binomial(theta = 2.53), fulldatatemp.noout)
#m0 m8 m18 m5 m6 m1 m3
#+change mean   *prop.overflow
aic.hydro.nb.best("totalbio", nocadata.noout,10)#prop overflow  ...because the best rain model was the null, we used this model selection for rain vs hydro
besttb<-glm.nb(round(totalbio*10)~log(maxvol)+site,  data = nocadata.noout)
bestmodtb<-glm.nb(round(totalbio*10)~log(maxvol)+site,  data = nocadata.noout)
hydrotb<-glm.nb(round(totalbio*10)~log(maxvol)+site*prop.overflow.days,  data = nocadata.noout); Anova(hydrotb, type=2) #sig
visreg(hydrotb, "prop.overflow.days", by="site",ylab="Total biomass", overlay=TRUE, partial=TRUE, band=FALSE)
nocadata.noout$resid.tb.hydro<-resid(glm.nb(round(totalbio*10)~log(maxvol)+site,data=nocadata.noout, na.action=na.exclude))
aic.percent$dhydrotrue[aic.percent$response=="Total Invertebrates"]<-Dsquared(glm(resid.tb.hydro~prop.overflow.days,family=gaussian, nocadata.noout), adjust=TRUE) #0
aic.percent$dhydrofalse[aic.percent$response=="Total Invertebrates"]<-Dsquared(glm(resid.tb.hydro~prop.overflow.days,family=gaussian, nocadata.noout), adjust=FALSE)
aic.percent$hydrotype[aic.percent$response=="Total Invertebrates"]<-"contingent"
aic.percent$bestmodel[aic.percent$response=="Total Invertebrates"]<-"hydrology"

#models for final line graphs

#plots ---------------
fig_N_rain<-visreg(bestnit, "k.scalar", by="site", ylab="Nitrogen uptake", overlay=TRUE, partial=TRUE, band=FALSE)
fig_sh_rain<-visreg(bestsh, "k.scalar", by="site", ylab="Shredder biomass", overlay=TRUE, partial=TRUE, band=FALSE)
fig_ff_rain<-visreg(bestff, "mu.scalar", by="site", ylab="Filter feeder biomass",overlay=TRUE, partial=TRUE, band=FALSE)
fig_sc_rain<-visreg(bestsc, "k.scalar", by="site", ylab="Scraper biomass",overlay=TRUE, partial=TRUE, band=FALSE)
fig_ga_rain<-visreg(bestga, "k.scalar", by="site", ylab="Gatherer biomass",overlay=TRUE, partial=TRUE, band=FALSE)
fig_pi_rain<-visreg(bestpi, "mu.scalar", by="site", ylab="Piercer biomass",overlay=TRUE, partial=TRUE, band=FALSE)
fig_ba_rain<-visreg(bestba, "mu.scalar", by="site", ylab="Bacterial biomass",overlay=TRUE, partial=TRUE, band=FALSE)

fig_decomp_hydro<-visreg(bestmoddecomp, "cv.depth", by="site", ylab="Decomposition", overlay=TRUE, partial=TRUE, band=FALSE)
fig_N_hydro<-visreg(bestmodN2, "prop.driedout.days", by="site",ylab="Nitrogen uptake", overlay=TRUE, partial=TRUE, band=FALSE)
fig_ff_hydro<-visreg(bestmodff, "cv.depth", by="site", ylab="Filter feeder biomass", overlay=TRUE, partial=TRUE, band=FALSE)
fig_sh_hydro<-visreg(besthydrosh, "prop.overflow.days", by="site", ylab="Shredder biomass", overlay=TRUE, partial=TRUE, band=FALSE)
fig_sc_hydro<-visreg(bestmodsc, "prop.driedout.days", by="site", ylab="Scraper biomass", overlay=TRUE, partial=TRUE, band=FALSE)
fig_ga_hydro<-visreg(hydroga2, "last_wet", by="site", ylab="Gatherer biomass", overlay=TRUE, partial=TRUE, band=FALSE)
fig_pi_hydro<-visreg(bestmodpi, "prop.overflow.days", by="site",ylab="Piercer biomass", overlay=TRUE, partial=TRUE, band=FALSE)

aic.percent$response <- factor(aic.percent$response, levels = c(response<-as.vector(c("CO2 flux","Decomposition", "Nitrogen uptake", "Bacterial density", "Total Invertebrates","Engulfer","Shredder", "Piercer","Gatherer","Scraper", "Filter feeder"))))


require(readr)

figure04_summary_list <- list(
  fig_N_rain = fig_N_rain,
  fig_sh_rain = fig_sh_rain,
  fig_ff_rain = fig_ff_rain,
  fig_sc_rain = fig_sc_rain,
  fig_ga_rain = fig_ga_rain,
  fig_pi_rain = fig_pi_rain,
  fig_ba_rain = fig_ba_rain,
  fig_decomp_hydro = fig_decomp_hydro,
  fig_N_hydro = fig_N_hydro,
  fig_ff_hydro = fig_ff_hydro,
  fig_sh_hydro = fig_sh_hydro,
  fig_sc_hydro = fig_sc_hydro,
  fig_ga_hydro = fig_ga_hydro,
  fig_pi_hydro = fig_pi_hydro,
  aic.percent = aic.percent
)

saveRDS(figure04_summary_list, "figure/data/figure04_summary_list.rds")



abs_k_Nfig<-visreg(absnit, "intended.k", by="site", ylab="Nitrogen uptake", xtrans=log, overlay=FALSE, partial=TRUE, rug=TRUE, band=FALSE)

Appendix7_summary_list<- list(
  abs_k_Nfig = abs_k_Nfig
)
saveRDS(Appendix7_summary_list, "figure/data/Appendix7_summary_list.rds")


#how should we deal with ns models? just box plots?

#overall figure

dodge <- position_dodge(width = 0.9)


ggplot(data = aic.percent, aes(x = response, y = draintrue))+
  geom_bar(data=aic.percent, aes(y=draintrue,x=response, fill=raintype), stat="identity", width = 0.75)+
  labs(title = "Effect of rainfall on ecosystem responses", 
       y = "Proportion residual deviance explained by rainfall (adj. for number of parameters)", x = "Response type")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

ggplot(data = aic.percent, aes(x = response, y = dhydrotrue))+
  geom_bar(data=aic.percent, aes(y=dhydrotrue,x=response, fill=hydrotype), stat="identity", width = 0.75)+
  labs(title = "Effect of hydrology on ecosystem responses", 
       y = "Proportion residual deviance explained by hydrology (adj. for number of parameters)", x = "Response type")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
