source("01_datasets_for_paper1.R")
source("02_functions.R")


#====correlation tests
#==note filter low in co, cr; engul in arg, co, piercer in co, cr,pr===

corr.fngrps<-function (dataset)
{out<-1:15
r.shredder<-resid(glm(shredder_bio~log(maxvol),family=poisson, data=dataset))
r.filter.feeder<-resid(glm(filter.feeder_bio~log(maxvol),family=poisson, data=dataset))
r.scraper<-resid(glm(scraper_bio~log(maxvol),family=poisson, data=dataset))
r.gatherer<-resid(glm(gatherer_bio~log(maxvol),family=poisson, data=dataset))
r.engulfer<-resid(glm(engulfer_bio~log(maxvol),family=poisson, data=dataset))
r.piercer<-resid(glm(piercer_bio~log(maxvol), family=poisson, data=dataset))
out[1]<-cor.test(r.shredder,r.filter.feeder, na.rm=TRUE)$estimate
out[2]<-cor.test(r.shredder,r.scraper, na.rm=TRUE)$estimate
out[3]<-cor.test(r.scraper,r.filter.feeder, na.rm=TRUE)$estimate
out[4]<-cor.test(r.shredder,r.gatherer, na.rm=TRUE)$estimate
out[5]<-cor.test(r.filter.feeder,r.gatherer, na.rm=TRUE)$estimate
out[6]<-cor.test(r.scraper,r.gatherer, na.rm=TRUE)$estimate
out[7]<-cor.test(r.engulfer,r.shredder, na.rm=TRUE)$estimate
out[8]<-cor.test(r.engulfer,r.filter.feeder, na.rm=TRUE)$estimate
out[9]<-cor.test(r.engulfer,r.scraper, na.rm=TRUE)$estimate
out[10]<-cor.test(r.engulfer,r.gatherer, na.rm=TRUE)$estimate
out[11]<-cor.test(r.piercer,r.shredder, na.rm=TRUE)$estimate
out[12]<-cor.test(r.piercer,r.filter.feeder, na.rm=TRUE)$estimate
out[13]<-cor.test(r.piercer,r.scraper, na.rm=TRUE)$estimate
out[14]<-cor.test(r.piercer,r.gatherer, na.rm=TRUE)$estimate
out[15]<-cor.test(r.piercer,r.engulfer, na.rm=TRUE)$estimate
return(out)
}


site<-c(rep("macae",15),rep("costarica", 15),rep("cardoso",15), rep("frenchguiana",15), rep("puertorico",15), rep("argentina",15), rep("colombia",15))
corrfg<-as.data.frame(site)
corrfg$fngrp1<-rep(c("filter.feeder", "scraper", "scraper", "gatherer", "gatherer", "gatherer", rep("engulfer",4), rep("piercer",5)),7)
corrfg$fngrp2<-rep(c("shredder", "shredder", "filter.feeder", "shredder", "filter.feeder", "scraper","shredder", "filter.feeder", "scraper", "gatherer","shredder", "filter.feeder", "scraper", "gatherer", "engulfer"),7)
corrfg$trophic<-as.factor(rep(c(rep("yes",6), rep("no",8), rep("no",1)),7))
corrfg$corr<-1:105
corrfg$corr[1:15]<-corr.fngrps(madata)
corrfg$corr[16:30]<-corr.fngrps(crdata)
corrfg$corr[31:45]<-corr.fngrps(cadata)
corrfg$corr[46:60]<-corr.fngrps(fgdata)
corrfg$corr[61:75]<-corr.fngrps(prdata)
corrfg$corr[76:90]<-corr.fngrps(ardata)
corrfg$corr[91:105]<-corr.fngrps(codata)

write.csv(corrfg, "Data/corrfg.csv")

mean(corrfg$corr)#r = 0.026
1.96*sd(corrfg$corr)/(105^0.5) # =/- 0.037

#site level
mean(tapply(corrfg$corr, corrfg$site, mean))#0.029
sd(tapply(corrfg$corr, corrfg$site, mean))#0.052
(sd(tapply(corrfg$corr, corrfg$site, mean)))/sqrt(7)#se is 0.020, twice this is 95% ci=0.038

mean(tapply(corrfg$corr[corrfg$trophic=="yes"], corrfg$site[corrfg$trophic=="yes"], mean))#0.077
sd(tapply(corrfg$corr[corrfg$trophic=="yes"], corrfg$site[corrfg$trophic=="yes"], mean))#0.036
(sd(tapply(corrfg$corr[corrfg$trophic=="yes"], corrfg$site[corrfg$trophic=="yes"], mean)))/sqrt(7)#se is 0.0137, twice this is 95% ci=0.027

mean(tapply(corrfg$corr[corrfg$trophic=="no"], corrfg$site[corrfg$trophic=="no"], mean))#-0.003
sd(tapply(corrfg$corr[corrfg$trophic=="no"], corrfg$site[corrfg$trophic=="no"], mean))#0.076
(sd(tapply(corrfg$corr[corrfg$trophic=="no"], corrfg$site[corrfg$trophic=="no"], mean)))/sqrt(7)#se is 0.028, twice this is 95% ci=0.056

mean(corrfg$corr[corrfg$trophic=="yes"&corrfg$site!="puertorico"])
mean(corrfg$corr[corrfg$trophic=="yes"])#r=0.080, or 0.078 without PR
1.96*sd(corrfg$corr[corrfg$trophic=="yes"])#0.35
mean(corrfg$corr[corrfg$trophic=="no"&corrfg$site!="puertorico"])#r=-0.0098, or -0.024 without PR
mean(corrfg$corr[corrfg$trophic=="no"])
1.96*sd(corrfg$corr[corrfg$trophic=="no"])#0.38
anova(lm(corr~trophic, data=corrfg)) #p=0.01882
Anova(lm(corr~trophic*site, data=corrfg), type=2) #p=0.020, diff overall between fn grps in same vs. diff trophic level
Anova(lm(corr~trophic*site, data=subset(corrfg, site!="puertorico")), type=2)#p=0.012

pval<-rep(1,1000)
for (i in 1:1000) {
  corrfg$ptrophic<-sample(corrfg$trophic, replace=FALSE)
  pval[i]<-Anova(lm(corr~ptrophic*site, data=subset(corrfg, site!="puertorico")), type=2)$Pr[1]#more sig
}

pval<-rep(1,1000)
for (i in 1:1000) {
  corrfg$ptrophic<-sample(corrfg$trophic, replace=FALSE)
  pval[i]<-Anova(lm(corr~ptrophic*site, data=corrfg), type=2)$Pr[1]#more sig
}
sort(pval)[25]#0.020...so our pvalue is just sig, but wouldn't be without excl. PR..remove flatworms?
sort(pval)[975]#0.9667

plot(corrfg$corr~corrfg$trophic, na.rm=TRUE)

subset(corrfg, site!="puertorico")
mean(corr.fngrps(madata)[1:6])#0.068
1.96*(sd(corr.fngrps(madata)[1:6]))/(6^0.5)#0.2
mean(corr.fngrps(madata)[7:14])#-0.09
1.96*(sd(corr.fngrps(madata)[7:14]))/(8^0.5)# 0.11

mean(corr.fngrps(crdata)[1:6])#0.12
1.96*(sd(corr.fngrps(crdata)[1:6]))/(6^0.5)#0.18
mean(corr.fngrps(crdata)[7:14])#-0.02
1.96*(sd(corr.fngrps(crdata)[7:14]))/(8^0.5)# 0.13

mean(corr.fngrps(fgdata)[1:6])#0.04
1.96*(sd(corr.fngrps(fgdata)[1:6]))/(6^0.5)#0.09
mean(corr.fngrps(fgdata)[7:14])#-0.09
1.96*(sd(corr.fngrps(fgdata)[7:14]))/(8^0.5)# 0.08

mean(corr.fngrps(prdata)[1:6])#0.09
1.96*(sd(corr.fngrps(prdata)[1:6]))/(6^0.5)#0.136
mean(corr.fngrps(prdata)[7:14])#0.10
1.96*(sd(corr.fngrps(prdata)[7:14]))/(8^0.5)# 0.178

#===plots

meancorr<-function(fn1,fn2)
{corrfg%>%filter(fngrp1==fn1 & fngrp2==fn2)%>%select(corr)%>%apply(2,mean.na)}

meancorr("filter.feeder","shredder")

fgpairs<-filter(corrfg, site=="macae")%>%select(fngrp1, fngrp2)

j<-1
for (j in 1:15)
{fgpairs[j,3]<-meancorr(fgpairs[j,1], fgpairs[j,2])}

View(fgpairs)

table.corr.fg<-tapply(corrfg$corr, corrfg$site, mean)%>%as.data.frame()
write.csv(table.corr.fg, "Data/tablecorrfg.csv")

write.csv(fgpairs, "Data/fgpairs.csv")

#======species interactions

#here I first considered single site models. Not much signal of rainfall in terms of prey biomass...
#guess this reflects the statistical averaging of diff detritivore effects

#this script wont run unless you convert "_bio" to "_biomass"!

aic.sitenb(round(ardata[,"prey_biomass"]*10), ardata) #m4, k+k2
aic.site(round(ardata[,"prey_biomass"]*10), negative.binomial(theta=1.8),ardata) 

aic.sitenb(round(cadata[,"prey_biomass"]*10), cadata) #m0  null
aic.site(round(cadata[,"prey_biomass"]*10), negative.binomial(theta=7.2),cadata)

aic.sitenb(round(codata[,"prey_biomass"]*10), codata) #m0  null
aic.site(round(codata[,"prey_biomass"]*10), negative.binomial(theta=1.37),codata)

aic.sitenb(round(crdata[,"prey_biomass"]*10), crdata) #m1 but not sig: mu
aic.site(round(crdata[,"prey_biomass"]*10), negative.binomial(theta=2.66),crdata)#still m1, m3 but m0 in delta<2

aic.sitenb(round(fgdata[,"prey_biomass"]*10), fgdata) #m0  null
aic.site(round(fgdata[,"prey_biomass"]*10), negative.binomial(theta=4.46),fgdata)

aic.sitenb(round(madata[,"prey_biomass"]*10), madata) #m0  null
aic.site(round(madata[,"prey_biomass"]*10), negative.binomial(theta=8.57),madata)

aic.sitenb(round(prdata[,"prey_biomass"]*10), prdata) #m0  null
aic.site(round(prdata[,"prey_biomass"]*10), negative.binomial(theta=3.64),prdata)

b1<-glm(scraper_bio~maxvol+engulfer_bio, data =madata, family = poisson)
summary(b1)
Anova(b1, type=2) #engulfer_bio sig and negative effect, but looks like its drien by two bromeliads
visreg(b1, "engulfer_bio")

b1<-glm(shredder_bio~log(maxvol)+predator_biomass, data =fgdata, family = poisson); anova(b1, test="Chi")#sig
b1<-glm(prey_biomass~log(maxvol)+predator_biomass, data =fgdata, family = poisson); anova(b1, test="Chi")#ns
b1<-glm(prey_biomass~log(maxvol)+(k.scalar+I(k.scalar^2))+predator_biomass, data =fgdata, family = poisson); anova(b1, test="Chi")

b1<-glm(prey_biomass~log(maxvol)+(k.scalar+I(k.scalar^2))+predator_biomass, data =ardata, family = poisson); anova(b1, test="Chi")#but POSITIVE effect of pred!
b1<-glm(prey_biomass~log(maxvol)+predator_biomass, data =ardata, family = poisson); anova(b1, test="Chi") #super sig

b1<-glm(shredder_bio~log(maxvol)+predator_biomass, data =crdata, family = poisson); anova(b1, test="Chi")#this is just sig, plot looks convincing
b1<-glm(shredder_bio~log(maxvol)+(mu.scalar)+predator_biomass, data =crdata, family = poisson); anova(b1, test="Chi")#this is just sig, plot looks convincing
b1<-glm(prey_biomass~log(maxvol)+predator_biomass, data =crdata, family = poisson); anova(b1, test="Chi")#this is just sig, plot looks convincing
b1<-glm(prey_biomass~log(maxvol)+(mu.scalar)+predator_biomass, data =crdata, family = poisson); anova(b1, test="Chi")#both mu and pred biomass sig

b1<-glm(prey_biomass~log(maxvol)+predator_biomass, data =madata, family = poisson); anova(b1, test="Chi") #p=0.07
b1<-glm(prey_biomass~log(maxvol)+Odonata_bio, data =madata, family = poisson); anova(b1, test="Chi") #ns
b1<-glm(prey_biomass~log(maxvol)+(k.scalar+I(k.scalar^2))+predator_biomass, data =madata, family = poisson); anova(b1, test="Chi")#only pred sig, not k

b1<-glm(prey_biomass~log(maxvol)+predator_biomass, data =cadata, family = poisson); anova(b1, test="Chi")#ns

b3<-glm(prey_biomassmass~log(maxvol)+predator_biomassmass*site, data =nocacoprdata, family = poisson); Anova(b1, type=2)
padj<-b3$deviance/b3$df.residual
b4<-glm(prey_biomassmass/padj~log(maxvol)+predator_biomassmass*site, data =nocoprdata, family = poisson)
Anova(b4, type=2)

b5<-glm.nb(round(prey_biomass*10)~log(maxvol)+predator_biomass*site, data =nocacoprdata); Anova(b5, type=2)#sig!
b6<-glm.nb(round(prey_biomass*10)~log(maxvol)+predator_biomass*site, data =filter(fulldata,site%in%c("costarica", "macae", "cardoso"))); Anova(b6, type=2) #no
b6<-glm.nb(round(prey_biomass*10)~log(maxvol)+site+predator_biomass, data =filter(fulldata,site%in%c("costarica", "macae"))); Anova(b6, type=2) #no

#note that these 3 sites have lots of odonates, but odonate biomass not the only driver (see below)
b6<-glm.nb(round(prey_biomassmass*10)~log(maxvol)+predator_biomassmass*site, data =filter(fulldata,site%nin%c("puertorico"))); Anova(b6, type=2)
b6<-glm.nb(round(prey_biomassmass*10)~log(maxvol)+predator_biomassmass*site, data =fulldata); Anova(b6, type=2) #hey this is sig in full data!
summary(b6)

b6<-glm.nb(round(prey_biomass*10)~log(maxvol)+site*(k.scalar+I(k.scalar^2))+predator_biomass*site, data =fulldata); Anova(b6, type=3) #
summary(b6)

#so, although predator biomass is correlated, this is driven by their common response to rainfall, not an effect of rainfall on predators than is transmitted to prey
b6a<-glm.nb(round(prey_biomass*10)~log(maxvol)+site*(k.scalar+I(k.scalar^2))+predator_biomass*site, data =fulldata) ; Anova(b6a, type=3)
b6b<-glm.nb(round(prey_biomass*10)~log(maxvol)+site*(k.scalar+I(k.scalar^2)), data =fulldata) ; Anova(b6b, type=3)
anova(b6a, b6b)#LR(df=7) =3.86, p=0.80, changed to LR (df=7) = 5.83, p=0.56 reported in ms

c6a<-glm.nb(round(shredder_bio*10)~log(maxvol)+site*(k.scalar+I(k.scalar^2))+predator_biomass*site, data =fulldata) ; Anova(c6a, type=3)
c6b<-glm.nb(round(shredder_bio*10)~log(maxvol)+site*(k.scalar+I(k.scalar^2)), data =fulldata) ; Anova(c6b, type=3)
anova(c6a, c6b)#LR(df=7) =3.86, p=0.80, changed to LR (df=7) = 5.83, p=0.56 reported in ms




tapply(fulldata$Odonata_bio, fulldata$site, datacheck)
b6<-glm.nb(round(prey_biomass*10)~log(maxvol)+Odonata_bio*site, data =filter(fulldata,site%in%c("costarica", "macae", "cardoso"))); Anova(b6, type=2) #no

b2<-rq((shredder_bio)^0.5~maxvol+predator_biomass,  tau = 0.05, data =madata)
summary(b2, se = "boot")
fullplot(crdata$predator_biomass, poisson, crdata)
fullsum(cacrmadata$Odonata_bio, poisson, cacrmadata)
m25a<-glm(predator_biomass~maxvol+(mu.scalar+I(mu.scalar^2))*(k.scalar+I(k.scalar^2)), family=poisson, data=madata);
padj<-m25a$deviance/m25a$df.residual
m25b<-glm(predator_biomass/padj~maxvol+(mu.scalar+I(mu.scalar^2))*(k.scalar+I(k.scalar^2)), family=poisson, data=madata);
plot((cadata$detritivore_bio)^0.5~(cadata$predator_biomass))
Anova(m25b, type=2, test="LR")

m25a<-glm(detritivore_bio~maxvol+(mu.scalar+I(mu.scalar^2))*(k.scalar+I(k.scalar^2)), family=poisson, data=ardata);
padj<-m25a$deviance/m25a$df.residual
m25b<-glm(detritivore_bio/padj~maxvol+(mu.scalar+I(mu.scalar^2))*(k.scalar+I(k.scalar^2)), family=poisson, data=ardata)
m25c<-glm(detritivore_bio/padj~maxvol, family=poisson, data=ardata)
anova(m25b, m25c, test= "LRT") # rainfall, largely k, explains 18.5 deviance
Anova(m25b, type = 2)
m25d<-glm(detritivore_bio/padj~maxvol+predator_biomass+(mu.scalar+I(mu.scalar^2))*(k.scalar+I(k.scalar^2)), family=poisson, data=ardata)
m25e<-glm(detritivore_bio/padj~maxvol+predator_biomass, family=poisson, data=ardata)
anova(m25d, m25e, test= "LRT")#rainfall explains 20.2 deviance??
anova(m25d, test= "LRT")#predator no longer sig!