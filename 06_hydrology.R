source("01_datasets_for_paper1.R")
source("02_functions.R")

#============hydrology==best aic models
Measure<-as.vector(c("Proportion days bromeliad dry", "Proportion days with overflow", "CV water depth", "Mean water depth", "Longest dry period", "Days since last wet", "Mean water temperature", "CV water temperature", "pH water"))
hydro.aic.percent<-as.data.frame(Measure)

aic.lmx((nocadata$prop.driedout.days)^0.33,  gaussian, nocadata)#bestmodel is m11
bestmod.propdried<-glm((prop.driedout.days)^0.33~log(maxvol)+site+log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)), family=gaussian, data = nocadata); Anova(bestmod.propdried,type=2)
#non sig interaction term removed
bestmod.propdried2<-glm((prop.driedout.days)^0.33~log(maxvol)+site+log(k.scalar)+(log(mu.scalar)+I(log(mu.scalar)^2)), family=gaussian, data = nocadata); Anova(bestmod.propdried2,type=2)
par(mfrow=c(1,1));visreg(bestmod.propdried, "mu.scalar", by="k.scalar", ylab="Proportion dried out", overlay=TRUE, partial=TRUE, band=FALSE)
par(mfrow=c(2,2)); plot(bestmod.propdried)

nocadata$resid.propdry<-resid(glm((prop.driedout.days)^0.33~log(maxvol)+site,data=nocadata, na.action=na.exclude))
bestpropdry.r<-glm(resid.propdry~log(k.scalar), family=gaussian, data = nocadata)
hydro.aic.percent$dsqtrue[hydro.aic.percent$Measure=="Proportion days bromeliad dry"]<-Dsquared(bestpropdry.r, adjust=TRUE)#0
hydro.aic.percent$raintype[hydro.aic.percent$Measure=="Proportion days bromeliad dry"]<-"general"


#costa rica seems to be very different - this site has very few rows of data so hydro variable
cr.propdried<-glm(sqrt(prop.driedout.days)~log(maxvol)+log(mu.scalar), family=gaussian, data = crdata)
outlierTest(cr.propdried)
par(mfrow=c(1,1));visreg(cr.propdried, "mu.scalar", ylab="Proportion dried out")#little trend and very variable up if anything
par(mfrow=c(2,2)); plot(cr.propdried)
aic.site(sqrt(crdata$prop.driedout.days), gaussian, crdata)#mo m2 (more affected by k than anything)
cr.propdriedk<-glm(sqrt(prop.driedout.days)~log(maxvol)+log(k.scalar), family=gaussian, data = crdata)
par(mfrow=c(1,1));visreg(cr.propdriedk, "k.scalar", ylab="Proportion dried out")#negative relationship, consistent with other sites

aic.lmx(sqrt(nocadata$prop.overflow.days),  gaussian, nocadata)#m0 and m4, m1, m2
bestmod.propoverflow<-glm(sqrt(prop.overflow.days)~log(maxvol)+site+log(k.scalar)+I(log(k.scalar)^2), family=gaussian, data = nocadata); Anova(bestmod.propoverflow, type=2)
par(mfrow=c(1,1));visreg(bestmod.propoverflow, "k.scalar", by="site",ylab="Proportion overflow", overlay=TRUE, partial=FALSE, band=FALSE)
par(mfrow=c(2,2)); plot(bestmod.propoverflow)#perfect diagostics

nocadata.15<-filter(nocadata, site_brom.id%nin%"argentina_15")
aic.lmx(log(nocadata.15$cv.depth),  gaussian, nocadata.15)
bestmod.cvdepth<-glm(log(cv.depth)~log(maxvol)+site*log(mu.scalar), family=gaussian, data = nocadata.15); Anova(bestmod.cvdepth, type=2)
par(mfrow=c(2,2)); plot(bestmod.cvdepth)
par(mfrow=c(1,1));visreg(bestmod.cvdepth, "mu.scalar", by="site",ylab="cv depth",trans=exp, overlay=TRUE, partial=TRUE, band=FALSE)
#ar ma and to a lesser degree pr have strongest decline in cv depth with mu; cr fg and co have no or variable change in cv with mu

nocadata.15<-filter(nocadata, site_brom.id%nin%"argentina_15")
aic.lmx(log(nocadata.15$cv.depth),  gaussian, nocadata.15)
bestmod.cvdepth<-glm(log(cv.depth)~log(maxvol)+site*log(mu.scalar), family=gaussian, data = nocadata.15); Anova(bestmod.cvdepth, type=2)
par(mfrow=c(2,2)); plot(bestmod.cvdepth)
par(mfrow=c(1,1));visreg(bestmod.cvdepth, "mu.scalar", by="site",ylab="cv depth",trans=exp, overlay=TRUE, partial=TRUE, band=FALSE)
#ar ma and to a lesser degree pr have strongest decline in cv depth with mu; cr fg and co have no or variable change in cv with mu

aic.lmx(log(nocadata.15$mean.depth),  gaussian, nocadata.15)# mu x site
bestmod.meandepth<-glm(log(mean.depth)~log(maxvol)+site*log(mu.scalar), family=gaussian, data = nocadata.15); Anova(bestmod.meandepth, type=2)
par(mfrow=c(2,2)); plot(bestmod.meandepth)

aic.lmx((nocadata$long_dry)^0.1,  gaussian, nocadata)#m9
bestmod.longdry<-glm((long_dry)^0.1~log(maxvol)+site+log(mu.scalar)*log(k.scalar), family=gaussian, data = nocadata); Anova(bestmod.longdry, type=3)
#but interaction not sig, removed:
bestmod.longdry2<-glm((long_dry)^0.1~log(maxvol)+site+log(mu.scalar)+log(k.scalar), family=gaussian, data = nocadata); Anova(bestmod.longdry2, type=2)
par(mfrow=c(1,1));visreg(bestmod.longdry, "mu.scalar", by="site",ylab="Longest dry period", overlay=TRUE, partial=FALSE, band=FALSE)
par(mfrow=c(1,1));visreg(bestmod.longdry, "k.scalar", by="site",ylab="Longest dry period", overlay=TRUE, partial=FALSE, band=FALSE)
visreg2d(bestmod.longdry, "mu.scalar", "k.scalar", zlab="Longest dry period",plot.type="persp")
par(mfrow=c(2,2)); plot(bestmod.longdry)

aic.lmx(nocadata$last_wet,  gaussian, nocadata)

aic.lmx(fulldata$mean_temp,  gaussian, fulldata)#m0 and m4, m1, m2
#bestmodel is m7
bestmod.meantemp<-glm(mean_temp~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data = fulldata); Anova(bestmod.meantemp, type=2)
visreg(bestmod.meantemp, "mu.scalar", by="site",ylab="Mean temp", overlay=FALSE, partial=FALSE, band=FALSE)

aic.lmx(fulldata$change_mean_temp,  gaussian, fulldata) #m7mu2 x site
bestmod.meantemp2<-glm(change_mean_temp~log(maxvol)+site*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data = fulldata); Anova(bestmod.meantemp2, type=2)
visreg(bestmod.meantemp2, "mu.scalar", by="site",ylab="Mean temp", overlay=FALSE, partial=TRUE, band=FALSE)

aic.lmx(fulldata$cv_mean_temp,  gaussian, fulldata)#m2 log k
bestmod.cvtemp<-glm(cv_mean_temp~log(maxvol)+site+log(k.scalar),family=gaussian, data = fulldata);Anova(bestmod.cvtemp, type=2)
aic.lmx(fulldata$change_cv_temp,  gaussian, fulldata) #m2 site+k
bestmod.cvtemp2<-glm(change_cv_temp~log(maxvol)+site+log(k.scalar),family=gaussian, data = fulldata)
visreg(bestmod.cvtemp, "k.scalar", by="site",ylab="cv temp", overlay=FALSE, partial=TRUE, band=FALSE)

aic.lmx(fulldata$ph.final,gaussian, fulldata)#no effect but residuals wonky
aic.lmx((fulldata$change_ph)^0.5,gaussian, fulldata) #needed to center ph on site first to get residuals to behave
bestmod.ph<-glm((change_ph)^0.5~log(maxvol)+site+log(mu.scalar),family=gaussian, data = fulldata); Anova(bestmod.ph, type=2)
par(mfrow=c(2,2)); plot(bestmod.ph) #truly ns

aic.lmx(fulldata$turbidity.final,gaussian, fulldata)#m3 site +(mu+mu2)
bestmod.turbid<-glm(log(turbidity.final)~log(maxvol)+site+log(mu.scalar)+I(log(mu.scalar)^2),family=gaussian, data = fulldata)
visreg(bestmod.turbid, "mu.scalar", by="site",ylab="Turbidity", overlay=FALSE, partial=TRUE, band=FALSE); #mu reduces turbid

aic.lmx(fulldata$oxygen.percent.final,gaussian, fulldata) #m5 log(mu)x site
bestmod.oxy<-glm(log(oxygen.percent.final)~log(maxvol)+site*log(mu.scalar),family=gaussian, data = fulldata)
visreg(bestmod.oxy, "mu.scalar", by="site",ylab="Oxygen", overlay=FALSE, partial=TRUE, band=FALSE); #mu reduces turbid
tapply(fulldata$oxygen.percent.final, fulldata$site, datacheck)#neither conc or percent in PR and CA
#goes down with mu in Arg and Mac, but up in other three countries..?
absmod.oxy<-glm(log(oxygen.percent.final)~log(maxvol)+site+log(intended.mu),family=gaussian, data = fulldata)
Anova(absmod.oxy, type=2)#no, not an absolute pattern

aic.lmx(waterchem.df$wc1,gaussian, waterchem.df) #mo so no effect of water chemistry
aic.lmx(waterchem.df$wc2,gaussian, waterchem.df) #log(mu)x site
bestmod.wc2<-glm(wc2~log(maxvol)+site*log(mu.scalar),family=gaussian, data = waterchem.df); Anova(bestmod.wc2, type=2)
#here site x mu sig effect on wc2 axis

tapply(fulldata$oxygen.percent.final, fulldata$site, datacheck) #none in cardoso or puerto rico
tapply(fulldata$ph.final, fulldata$site, datacheck)#good
tapply(fulldata$turbidity.final, fulldata$site, datacheck) #not in fg


#=================site level hydrology models

d.dried<-as.vector(rep(NA,7))
sensitive<-as.data.frame(d.dried)
sensitive$d.over<-rep(NA,7)
sensitive$d.cvdepth<-rep(NA,7)
sensitive$site<-c("argentina", "cardoso", "colombia", "costarica","frenchguiana", "macae", "puertorico")


aic.site((ardata$prop.driedout.days)^0.33,  gaussian, ardata)#m9
ardata$resid.prop.dried<-resid(glm((prop.driedout.days)^0.33~log(maxvol), data=ardata, na.action=na.exclude))
Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=ardata), adjust=TRUE) #0.467
sensitive$d.dried[sensitive$site=="argentina"]<-Dsquared(glm(resid.prop.dried~log(k.scalar)*log(mu.scalar),family=gaussian, data=ardata), adjust=TRUE) #0.417
Anova(glm(resid.prop.dried~log(k.scalar)*log(mu.scalar),family=gaussian, data=ardata), type=2) #sig

aic.site((codata$prop.driedout.days)^0.33,  gaussian, codata)#m0
codata$resid.prop.dried<-resid(glm((prop.driedout.days)^0.33~log(maxvol), data=codata, na.action=na.exclude))
sensitive$d.dried[sensitive$site=="colombia"]<-Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=codata), adjust=TRUE) #-0.15
Anova(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=codata), type=2) #ns

aic.site((crdata$prop.driedout.days)^0.33,  gaussian, crdata)#m0
crdata$resid.prop.dried<-resid(glm((prop.driedout.days)^0.33~log(maxvol), data=crdata, na.action=na.exclude))
sensitive$d.dried[sensitive$site=="costarica"]<-Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=crdata), adjust=TRUE) #-0.19
Anova(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=crdata), type=2) #ns

aic.site((fgdata$prop.driedout.days)^0.33,  gaussian, fgdata)#m3
fgdata$resid.prop.dried<-resid(glm((prop.driedout.days)^0.33~log(maxvol), data=fgdata, na.action=na.exclude))
Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), adjust=TRUE) #+0.14
sensitive$d.dried[sensitive$site=="frenchguiana"]<-Dsquared(glm(resid.prop.dried~(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), adjust=TRUE) #+0.20 with m3
Anova(glm(resid.prop.dried~(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), type=2) #sig as type 1, marg ns type 2

aic.site((madata$prop.driedout.days)^0.33,  gaussian, madata)#m3
madata$resid.prop.dried<-resid(glm((prop.driedout.days)^0.33~log(maxvol), data=madata, na.action=na.exclude))
Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=madata), adjust=TRUE) #+0.50
sensitive$d.dried[sensitive$site=="macae"]<-Dsquared(glm(resid.prop.dried~log(mu.scalar)+I(log(mu.scalar)^2),family=gaussian, data=madata), adjust=TRUE) #+0.56
Anova(glm(resid.prop.dried~log(mu.scalar)+I(log(mu.scalar)^2),family=gaussian, data=madata), type=2) #sig

aic.site((prdata$prop.driedout.days)^0.33,  gaussian, prdata)#m4
prdata$resid.prop.dried<-resid(glm((prop.driedout.days)^0.33~log(maxvol), data=prdata, na.action=na.exclude))
Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=prdata), adjust=TRUE) #0.001
sensitive$d.dried[sensitive$site=="puertorico"]<-Dsquared(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, data=prdata), adjust=TRUE) #0.11
Anova(glm(resid.prop.dried~(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, data=prdata), type=2) #sig

aic.site(log(ardata$cv.depth),  gaussian, ardata)#m9
ardata$resid.cv.depth<-resid(glm(log(cv.depth)~log(maxvol), data=ardata, na.action=na.exclude))
Dsquared(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=ardata), adjust=TRUE) #0.38
sensitive$d.cvdepth[sensitive$site=="argentina"]<-Dsquared(glm(resid.cv.depth~log(k.scalar)*log(mu.scalar),family=gaussian, data=ardata), adjust=TRUE) #0.42
Anova(glm(resid.cv.depth~log(k.scalar)*log(mu.scalar),family=gaussian, data=ardata), type=2) #sig

aic.site(log(codata$cv.depth),  gaussian, codata)#m0
codata$resid.cv.depth<-resid(glm(log(cv.depth)~log(maxvol), data=codata, na.action=na.exclude))
sensitive$d.cvdepth[sensitive$site=="colombia"]<-Dsquared(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=codata), adjust=TRUE) #-0.20
Anova(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=codata), type=2) #ns

aic.site(log(crdata$cv.depth),  gaussian, crdata)#m0
crdata$resid.cv.depth<-resid(glm(log(cv.depth)~log(maxvol), data=crdata, na.action=na.exclude))
sensitive$d.cvdepth[sensitive$site=="costarica"]<-Dsquared(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=crdata), adjust=TRUE) #-0.21
Anova(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=crdata), type=2) #ns

aic.site(log(fgdata$cv.depth),  gaussian, fgdata)#m11
fgdata$resid.cv.depth<-resid(glm(log(cv.depth)~log(maxvol), data=fgdata, na.action=na.exclude))
Dsquared(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), adjust=TRUE) #0.36
sensitive$d.cvdepth[sensitive$site=="frenchguiana"]<-Dsquared(glm(resid.cv.depth~log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), adjust=TRUE) #0.26
Anova(glm(resid.cv.depth~log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), type=3) #sig

aic.site(log(madata$cv.depth),  gaussian, madata)#m1
madata$resid.cv.depth<-resid(glm(log(cv.depth)~log(maxvol), data=madata, na.action=na.exclude))
Dsquared(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=madata), adjust=TRUE) #0.54
sensitive$d.cvdepth[sensitive$site=="macae"]<-Dsquared(glm(resid.cv.depth~log(mu.scalar),family=gaussian, data=madata), adjust=TRUE) #0.55
Anova(glm(resid.cv.depth~log(mu.scalar),family=gaussian, data=madata), type=2) #sig

aic.site(log(prdata$cv.depth),  gaussian, prdata)#m4
prdata$resid.cv.depth<-resid(glm(log(cv.depth)~log(maxvol), data=prdata, na.action=na.exclude))
Dsquared(glm(resid.cv.depth~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=prdata), adjust=TRUE) #0.04
sensitive$d.cvdepth[sensitive$site=="puertorico"]<-Dsquared(glm(resid.cv.depth~log(k.scalar)+I(log(k.scalar)^2),family=gaussian, data=prdata), adjust=TRUE) #0.11
Anova(glm(resid.cv.depth~log(k.scalar)+I(log(k.scalar)^2),family=gaussian, data=prdata), type=2) #sig

aic.site(sqrt(ardata$prop.overflow.days),  gaussian, ardata)#m2
ardata$resid.prop.over<-resid(glm(sqrt(prop.overflow.days)~log(maxvol), data=ardata, na.action=na.exclude))
Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=ardata), adjust=TRUE) #0.24
sensitive$d.over[sensitive$site=="argentina"]<-Dsquared(glm(resid.prop.over~log(k.scalar),family=gaussian, data=ardata), adjust=TRUE) #0.31
Anova(glm(resid.prop.over~log(k.scalar),family=gaussian, data=ardata), type=2) #sig

aic.site(sqrt(codata$prop.overflow.days),  gaussian, codata)#m0
codata$resid.prop.over<-resid(glm(sqrt(prop.overflow.days)~log(maxvol), data=codata, na.action=na.exclude))
sensitive$d.over[sensitive$site=="colombia"]<-Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=codata), adjust=TRUE) #-0.20
Anova(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=codata), type=2) #sig

aic.site(sqrt(crdata$prop.overflow.days),  gaussian, crdata)#m0
crdata$resid.prop.over<-resid(glm(sqrt(prop.overflow.days)~log(maxvol), data=crdata, na.action=na.exclude))
sensitive$d.over[sensitive$site=="costarica"]<-Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=crdata), adjust=TRUE) #-0.08
Anova(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=crdata), type=2) #sig

aic.site(sqrt(fgdata$prop.overflow.days),  gaussian, fgdata)#m0
fgdata$resid.prop.over<-resid(glm(sqrt(prop.overflow.days)~log(maxvol), data=fgdata, na.action=na.exclude))
sensitive$d.over[sensitive$site=="frenchguiana"]<-Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), adjust=TRUE) #-0.07
Anova(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=fgdata), type=3) #ns (one term is 0.09)

aic.site(sqrt(madata$prop.overflow.days),  gaussian, madata)#m9
madata$resid.prop.over<-resid(glm(sqrt(prop.overflow.days)~log(maxvol), data=madata, na.action=na.exclude))
Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=madata), adjust=TRUE) #0.37
sensitive$d.over[sensitive$site=="macae"]<-Dsquared(glm(resid.prop.over~log(k.scalar)*log(mu.scalar),family=gaussian, data=madata), adjust=TRUE) #0.14
Anova(glm(resid.prop.over~log(k.scalar)*log(mu.scalar),family=gaussian, data=madata), type=3) #sig

aic.site(sqrt(prdata$prop.overflow.days),  gaussian, prdata)#m4
prdata$resid.prop.over<-resid(glm(sqrt(prop.overflow.days)~log(maxvol), data=prdata, na.action=na.exclude))
Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2))*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, data=prdata), adjust=TRUE) #0.12
sensitive$d.over[sensitive$site=="puertorico"]<-Dsquared(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, data=prdata), adjust=TRUE) #0.22
Anova(glm(resid.prop.over~(log(k.scalar)+I(log(k.scalar)^2)),family=gaussian, data=prdata), type=2) #sig

#====how hydrologically variable are ambient conditions in sites?
#lets average over mu 0.6 to 1.5
ambient<-filter(nocadata, mu.scalar >= 0.6 & mu.scalar<1.5&k.scalar==1)
ambient.driedout<-as.vector(tapply(ambient$prop.driedout.days, ambient$site,mean))
ambient.overflow<-as.vector(tapply(ambient$prop.overflow.days, ambient$site,mean))
ambient.cvdepth<-as.vector(tapply(ambient$cv.depth, ambient$site,mean))

write.csv(sensitive, "Data/sensitive.csv")


#====hydrology correlations

nocadata.noout.hydro<-filter(nocadata, site_brom.id%nin%"argentina_15")%>%
 select(prop.driedout.days, prop.overflow.days, cv.depth, long_dry, last_wet, mean.depth)
cor(nocadata.noout.hydro)



#===hydrology plots


fig_cvdepth_rain<-visreg(bestmod.cvdepth, "mu.scalar", by="site", ylab="CV water depth", overlay=TRUE, partial=TRUE, band=FALSE)

require(readr)

figure05_hydro_list <- list(
  fig_cvdepth_rain = fig_cvdepth_rain
)

saveRDS(figure05_hydro_list, "figure/data/figure05_hydro_list.rds")
