source("01_datasets_for_paper1.R")
source("02_functions.R")
source("06_hydrology.R")

#===========repeat family analysis with aic-c model selection

sensitive <- read_csv("Data/sensitive.csv")

taxa<-as.vector(c("Ceratopogonidae_bio", "Chironomidae_bio","Coenagrionidae_bio", "Corethrellidae_bio", "Culicidae_bio", 
                  "Empididae_bio", "Limnocytheridae_bio", "Psychodidae_bio", "Scirtidae_bio","Syrphidae_bio", "Tabanidae_bio", 
                  "TipulidaeLimoniidae_bio", "Oligochaeta_bio", "Hirudinea_bio", "Calamoceratidae_bio","Cecidomyiidae_bio", "Candonidae_bio",
                  "Hydrophilidae_bio", "Pseudostigmatidae_bio","Stratiomyidae_bio","Ephydridae_bio", "Dytiscidae_bio"))
taxa.percent<-as.data.frame(taxa)

aa<-"Ceratopogonidae_bio"; bb<-noargprdata; cc<-filter(noargcaprdata,site_brom.id%nin%"argentina_15"); dd<-100
aic.lmx.nb.best(aa, bb, dd)#m9 (mu x k)m11
aic.hydro.nb.best(aa, cc, dd)#m4 (mean.depth) m6 m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(mu.scalar)*log(k.scalar), family=gaussian, bb), adjust=TRUE) #0.04
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~mean.depth,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(nofgdata$Chironomidae_bio*100), nofgdata)#m5 (site * mu) m3
aa<-"Chironomidae_bio"; bb<-nofgdata; cc<-filter(nocafgdata,site_brom.id%nin%"argentina_15"); dd<-100
aic.lmx.nb.best(aa, bb, dd)#m5 site x mu
aic.hydro.nb.best(aa, cc, dd)#m22 site x last wet
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~site*log(mu.scalar), family=gaussian, bb), adjust=TRUE) #0.04
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~site*last_wet,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(camadata$Coenagrionidae_bio*2),  camadata)# best thetas
aic.lmxnb.x(round(camadata$Coenagrionidae_bio*2),  camadata)#if remove volume, this is m24: site*(k+k2)
aa<-"Coenagrionidae_bio"; bb<-camadata; cc<-filter(madata,site_brom.id%nin%"argentina_15"); dd<-2
aic.lmx.nb.best(aa, bb, dd)# m8 (k+k2)
aic.sitehydro.nb.best(aa, cc, dd)# m3 prop.driedout.days 
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~1,data=cc, na.action=na.exclude))#not sig: plot(madata$Coenagrionidae_bio~madata$maxvol)
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~site*(log(k.scalar)+(I(log(k.scalar)^2))), family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~prop.driedout.days,family=gaussian,cc), adjust=TRUE) 

aic.lmxnb(round(cafgmadata$Corethrellidae_bio*100),  cafgmadata)#m9 (site+mu x k)
aa<-"Corethrellidae_bio"; bb<-cafgmadata; cc<-filter(fgmadata,site_brom.id%nin%"argentina_15"); dd<-100
aic.lmx.nb.best(aa, bb, dd)#m9 (mu x k)
aic.hydro.nb.best(aa, cc, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(mu.scalar)*log(k.scalar), family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(nocodata$Culicidae_bio*100),  nocodata)#m7 (site*mu+mu2)
aa<-"Culicidae_bio"; bb<-nocodata; cc<-filter(nocacodata,site_brom.id%nin%"argentina_15"); dd<-100
aic.lmx.nb.best(aa, bb, dd)#(site*mu+mu2)
aic.hydro.nb.best(aa, cc, dd)#m19 site*prop,driedout.days
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~site*(log(mu.scalar)+(I(log(mu.scalar)^2))), family=gaussian, bb), adjust=TRUE) #0.04
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~site*prop.driedout.days,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(camadata$Empididae_bio*100),  camadata)#m7 (site*mu+mu2) m5
aa<-"Empididae_bio"; bb<-camadata; cc<-madata; dd<-100
aic.lmx.nb.best(aa, bb, dd)#m7(site*mu+mu2)
aic.sitehydro.nb.best(aa, cc, dd)#m3 prop.driedout.days
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol),data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~site*(log(mu.scalar)+(I(log(mu.scalar)^2))), family=gaussian, bb), adjust=TRUE) #0.04
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~prop.driedout.days,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(noargcocrdata$Limnocytheridae_bio*100),  noargcocrdata)#m6 (site x k)
aa<-"Limnocytheridae_bio"; bb<-noargcocrdata; cc<-fgmaprdata; dd<-100
aic.lmx.nb.best(aa, bb, dd)#m6 (site x k)
aic.hydro.nb.best(aa, cc, dd)#m21(site x long_dry)
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~site*log(k.scalar), family=gaussian, bb), adjust=TRUE) #0.04
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~site*long_dry,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(noargcodata$Naididae_bio*100),  noargcodata) #m0 (no rainfall effect)
#see oligochaete

aic.lmxnb(round(argcacodata$Psychodidae_bio*10),  argcacodata)# gives m12 m0 m10 m6 m2
aic.lmx(round(argcacodata$Psychodidae_bio*100), family=negative.binomial(theta = 0.2222), argcacodata)#m0 (no rainfall effect)
aa<-"Psychodidae_bio"; bb<-argcacodata; cc<-filter(arcodata,site_brom.id%nin%"argentina_15"); dd<-10
aic.lmx.nb.best(aa, bb, dd)# m12 (k=k2)(mu+mu2)
aic.hydro.nb.best(aa, cc, dd)# m0 m1 m3 m2 m5
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~(log(k.scalar)+(I(log(k.scalar)^2)))*(log(mu.scalar)+(I(log(mu.scalar)^2))), family=gaussian, bb), adjust=TRUE) 
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian, cc), adjust=TRUE) 


aic.lmxnb(round(fulldata$Scirtidae_bio*100),  fulldata)#m8 m0 m2 m6 (no rainfall effect)
aa<-"Scirtidae_bio"; bb<-fulldata; cc<-filter(nocadata,site_brom.id%nin%"argentina_15"); dd<-100
aic.lmx.nb.best(aa, bb, dd)#m8 sitex(k+k2)
aic.hydro.nb.best(aa, cc, dd)#
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~site*(log(k.scalar)+(I(log(k.scalar)^2))), family=gaussian, bb), adjust=TRUE) 
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~long_dry,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(cocrdata$Syrphidae_bio*1),  cocrdata)#m2 (k) m0 m6 m1 m4 #conclude no rainfall effect
aic.lmx(round(cocrdata$Syrphidae_bio*1), poisson, cocrdata) #not as good a model, m16
aa<-"Syrphidae_bio"; bb<-cocrdata; cc<-filter(cocrdata,site_brom.id%nin%"argentina_15"); dd<-1
aic.lmx.nb.best(aa, bb, dd)# m2 (k) 
aic.hydro.nb.best(aa, cc, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(k.scalar), family=gaussian, bb), adjust=TRUE) 
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian, cc), adjust=TRUE) 


aic.lmxnb(round(argcrdata$Tabanidae_bio*1),  argcrdata)#m2, but note that theta blows up for half of models inc m0
tapply(fulldata$Tabanidae_bio, fulldata$site, datacheck)#argentina now down to 13 brom, after fixed
aa<-"Tabanidae_bio"; bb<-argcrdata; cc<-filter(argcrdata,site_brom.id%nin%"argentina_15"); dd<-1
aic.lmx.nb.best(aa, bb, dd)# m0
aic.hydro.nb.best(aa, cc, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~site,data=bb, na.action=na.exclude)) #not removed maxvol as destabilizing models
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1, family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(fulldata$TipulidaeLimoniidae_bio*10),  fulldata)#m4 (site +k+k2) m2
aa<-"TipulidaeLimoniidae_bio"; bb<-fulldata; cc<-filter(nocadata,site_brom.id%nin%"argentina_15"); dd<-10
aic.lmx.nb.best(aa, bb, dd)#m4 k+k2
aic.hydro.nb.best(aa, cc, dd)#m22 site x last_wet
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(k.scalar)+(I(log(k.scalar)^2)), family=gaussian, bb), adjust=TRUE) 
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~site*last_wet,family=gaussian, cc), adjust=TRUE) 

aic.lmxnb(round(noargcodata$Oligochaeta_bio*1000),  noargcodata)#m0 m1 (no rainfall effect)
aa<-"Oligochaeta_bio"; bb<-noargcodata; cc<-noargcacodata; dd<-1000
aic.lmx.nb.best(aa, bb, dd)#m0
aic.hydro.nb.best(aa, cc, dd)#m6 last_wet
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1, family=gaussian, bb), adjust=TRUE) 
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~last_wet,family=gaussian, cc), adjust=TRUE) 



#these are for interest, not included in d-squared table for now
aic.lmxnb(round(camaprdata$Tanypodinae_bio*1000),  camaprdata)#m1 m0 m3 (no rainfall effect)
aa<-"Tanypodinae_bio"; bb<-camaprdata; cc<-maprdata; dd<-200
aic.lmx.nb.best(aa, bb, dd)#m1 negative effect of mu
aic.hydro.nb.best(aa, cc, dd)#m22
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol)+site,data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol)+site,data=cc, na.action=na.exclude))
#taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(mu.scalar), family=gaussian, bb), adjust=TRUE)#-0.01 
#taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~site*last_wet,family=gaussian, cc), adjust=TRUE) #0.22
Dsquared(glm(cc[,"resid.taxa"]~site*prop.overflow.days,family=gaussian, cc), adjust=TRUE)#0.12


aic.lmxnb(round(cafgdata$predCerato_bio*1000),  cafgdata)#m0 m1 (no rainfall effect)
aic.lmxnb(round(cacrmadata$detCerato_bio*1000),  cacrmadata) #m9 m0 m2 m1 (no rainfall effect)
aic.lmxnb(round(nocodata$Culex*100),  nocodata) #m7 (site*mu+mu2) m13
aic.lmxnb(round(noargcacodata$Wyeomyia*100),  noargcacodata)#m7 (site*mu+mu2)
aic.lmxnb(round(nocodata$filtCulicidae*100),  nocodata)#m7 (site*mu+mu2)
aic.lmxnb(round(nofgdata$detChiron_bio*100),  nofgdata)#m5 (site * mu)

#single site families
aic.sitenb(round(cadata$Calamoceratidae_bio*4),  cadata)#m0
aa<-"Calamoceratidae_bio"; bb<-cadata; cc<-filter(cadata,site_brom.id%nin%"argentina_15"); dd<-4 #note no hydro data
aic.siterain.nb.best(aa, bb, dd)# m0
taxa.percent$draintrue[taxa.percent$taxa==aa]<-0 #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-NA


aic.sitenb(round(fgdata$Cecidomyiidae_bio*100),  fgdata) #m3 m0 m2 (no rainfall effect)
aa<-"Cecidomyiidae_bio"; bb<-fgdata; cc<-fgdata; dd<-100
aic.siterain.nb.best(aa, bb, dd)# m3 (mu + mu2)
aic.sitehydro.nb.best(aa, cc, dd)# m0 (null) m2 m1
bb[,"resid.taxa"]<-resid(glm(round(bb[,aa]*dd)~log(maxvol),family=negative.binomial(theta=0.2344),data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm(round(cc[,aa]*dd)~log(maxvol),family=negative.binomial(theta=0.2344),data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~(log(mu.scalar)+(I(log(mu.scalar)^2))), family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian,cc), adjust=TRUE) 

aic.sitenb(round(prdata$Candonidae_bio*100),  prdata)#m2 m4 m0 (no rainfall effect)
aa<-"Candonidae_bio"; bb<-prdata; cc<-prdata; dd<-100
aic.siterain.nb.best(aa, bb, dd)# m2: k
aic.sitehydro.nb.best(aa, cc, dd)# m3 prop.driedout.days
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol),data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(k.scalar),family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~prop.driedout.days,family=gaussian,cc), adjust=TRUE) 

aic.sitenb(round(prdata$Enchytraeoidae_bio*100),  prdata) #m0 m1 (no rainfall effect)
aa<-"Enchytraeoidae_bio"; bb<-prdata; cc<-prdata; dd<-100
#aic.siterain.nb.best(aa, bb, dd)# m0 from above
aic.sitehydro.nb.best(aa, cc, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol),data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1,family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian,cc), adjust=TRUE) 

aic.sitenb(round(crdata$Hydrophilidae_bio*1),  crdata) #m0 (no rainfall effect)
aa<-"Hydrophilidae_bio"; bb<-crdata; cc<-crdata; dd<-1
aic.siterain.nb.best(aa, bb, dd)#m0
aic.sitehydro.nb.best(aa, cc, dd)#m1 cv.depth
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol),data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1,family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~cv.depth,family=gaussian,cc), adjust=TRUE) 

aic.sitenb(round(crdata$Pseudostigmatidae_bio*3),  crdata) #m0 m4
aa<-"Pseudostigmatidae_bio"; bb<-crdata; cc<-crdata; dd<-3
aic.siterain.nb.best(aa, bb, dd)#m0 (null) m4
aic.sitehydro.nb.best(aa, cc, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
cc[,"resid.taxa"]<-resid(glm.nb(round(cc[,aa]*dd)~log(maxvol),data=cc, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1,family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-Dsquared(glm(cc[,"resid.taxa"]~1,family=gaussian,cc), adjust=TRUE) 

aic.sitenb(round(cadata$Stratiomyidae_bio*10),  cadata)#m0
aa<-"Stratiomyidae_bio"; bb<-cadata; cc<-cadata; dd<-10
aic.siterain.nb.best(aa, bb, dd)#m0 (null) m1 m2 m3
#aic.sitehydro.nb.best(aa, cc, dd)#m0  undefined as cardodo no hydro
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1,family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-NA 


aic.sitenb(round(cadata$Ephydridae_bio*100),  cadata)#m9 mu * k) m10
aa<-"Ephydridae_bio"; bb<-cadata;  dd<-100
aic.siterain.nb.best(aa, bb, dd)#m9 (mu x k) m10
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(k.scalar)*log(mu.scalar),family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-NA


aic.sitenb(round(cadata$Dytiscidae_bio*20),  cadata)#m0 (no rainfall effect) m3
aa<-"Dytiscidae_bio"; bb<-cadata;  dd<-20
aic.siterain.nb.best(aa, bb, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~1,family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-NA



aic.sitenb(round(cadata$Anopheles*100),  cadata) #m3 (mu = mu2) m4

aic.sitenb(round(cadata$Hirudinea*40),  cadata)#
aa<-"Hirudinea_bio"; bb<-cadata;  dd<-40
aic.siterain.nb.best(aa, bb, dd)#m0
bb[,"resid.taxa"]<-resid(glm.nb(round(bb[,aa]*dd)~log(maxvol),data=bb, na.action=na.exclude))
taxa.percent$draintrue[taxa.percent$taxa==aa]<-Dsquared(glm(bb[,"resid.taxa"]~log(k.scalar)*(log(mu.scalar)+I(log(mu.scalar)^2)),family=gaussian, bb), adjust=TRUE) #
taxa.percent$dhydrotrue[taxa.percent$taxa==aa]<-NA

##calculating biomass-weighted sensitivity index of site species pools

#shortcut to taxa.percent next line
#taxa.percent<-read.csv("C:/Users/Diane/Dropbox/BWG Drought Experiment/Paper 1_thresholds/Data/family_sensitivity.csv")

ar.poolmaker<-function(yvar){sum(ardata[,yvar], na.rm=TRUE)}
ca.poolmaker<-function(yvar){sum(cadata[,yvar], na.rm=TRUE)}
co.poolmaker<-function(yvar){sum(codata[,yvar], na.rm=TRUE)}
cr.poolmaker<-function(yvar){sum(crdata[,yvar], na.rm=TRUE)}
fg.poolmaker<-function(yvar){sum(fgdata[,yvar], na.rm=TRUE)}
ma.poolmaker<-function(yvar){sum(madata[,yvar], na.rm=TRUE)}
pr.poolmaker<-function(yvar){sum(prdata[,yvar], na.rm=TRUE)}

taxavector<-as.vector(taxa.percent$taxa)
taxa.percent$ar.pool.biomass<-as.vector(sapply(taxavector, ar.poolmaker))
taxa.percent$ca.pool.biomass<-as.vector(sapply(taxavector, ca.poolmaker))
taxa.percent$co.pool.biomass<-as.vector(sapply(taxavector, co.poolmaker))
taxa.percent$cr.pool.biomass<-as.vector(sapply(taxavector, cr.poolmaker))
taxa.percent$fg.pool.biomass<-as.vector(sapply(taxavector, fg.poolmaker))
taxa.percent$ma.pool.biomass<-as.vector(sapply(taxavector, ma.poolmaker))
taxa.percent$pr.pool.biomass<-as.vector(sapply(taxavector, pr.poolmaker))

site<-as.vector(c("argentina", "cardoso","colombia", "costarica", "frenchguiana","macae", "puertorico" ))

sens.index<-as.vector(c(
  sum(taxa.percent$ar.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$ar.pool.biomass),
  sum(taxa.percent$ca.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$ca.pool.biomass),
  sum(taxa.percent$co.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$co.pool.biomass),
  sum(taxa.percent$cr.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$cr.pool.biomass),
  sum(taxa.percent$fg.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$fg.pool.biomass),
  sum(taxa.percent$ma.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$ma.pool.biomass),
  sum(taxa.percent$pr.pool.biomass*taxa.percent$draintrue)/sum(taxa.percent$pr.pool.biomass)
))

sens.index.hydro<-as.vector(c(
  sum(taxa.percent$ar.pool.biomass*taxa.percent$dhydrotrue, na.rm=TRUE)/sum(taxa.percent$ar.pool.biomass),
  NA,
  sum(taxa.percent$co.pool.biomass*taxa.percent$dhydrotrue, na.rm=TRUE)/sum(taxa.percent$co.pool.biomass),
  sum(taxa.percent$cr.pool.biomass*taxa.percent$dhydrotrue, na.rm=TRUE)/sum(taxa.percent$cr.pool.biomass),
  sum(taxa.percent$fg.pool.biomass*taxa.percent$dhydrotrue, na.rm=TRUE)/sum(taxa.percent$fg.pool.biomass),
  sum(taxa.percent$ma.pool.biomass*taxa.percent$dhydrotrue, na.rm=TRUE)/sum(taxa.percent$ma.pool.biomass),
  sum(taxa.percent$pr.pool.biomass*taxa.percent$dhydrotrue, na.rm=TRUE)/sum(taxa.percent$pr.pool.biomass)
))


sensitive$sens.index<-sens.index
sensitive$sens.index.hydro<-sens.index.hydro
sensitive$ambient.driedout<-ambient.driedout
sensitive$ambient.overflow<-ambient.overflow
sensitive$ambient.cvdepth<-ambient.cvdepth

cor.test(sensitive$sens.index,sensitive$ambient.driedout) 
#t = -4.2189, df = 4, p-value = 0.01349  r = -0.90...so sites that dry out freq unde ambient conditions have species pools that are pre-adapted to be insensitive

cor.test(sensitive$ambient.driedout,sensitive$d.dried)
# p =0.39, r = -0.43 ...but sites that dry out more under ambient conditions are not more sensitive to rainfall change

write.csv(sensitive, "Data/family_sensitivity.csv")


ggplot(sensitive, aes(x=ambient.driedout, y=sens.index)) + geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  labs(x="Ambient proportion of days without water", y = "Species pool sensitivity to rain")+
  theme_bw()
