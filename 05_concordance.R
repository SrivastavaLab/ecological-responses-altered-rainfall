source("01_datasets_for_paper1.R")
source("02_functions.R")

#===concordance===memory intensive script!==will take a while to run

library(Agreement)

preds<-c(12,2,3,5,4,5,6,7,3)
preds<-c(NA,NA,NA,NA,NA,NA,NA,NA)
a<-c(2,5,3,1,7,8,11,12)
obs<-log(round(a)+0.001)

preds<-rnorm(30, mean=10, sd=2)
obs<-rnorm(30, mean=10, sd=2)
predobs<-as.data.frame(cbind(preds,obs))
xx<-as.vector(predobs$preds)
yy<-as.vector(predobs$obs)
sc.agr<-agreement(x=xx, y=yy, error="const",TDI_a = 1, target="random")
sc.agr$Estimate$CCC
sc.agr$Conf_Limit$CCC
sc.agr$Estimate$Precision
sc.agr$Conf_Limit$Precision
sc.agr$Estimate$Accuracy
sc.agr$Conf_Limit$Accuracy
predobs$id<-1:30
long<-gather(predobs,type, value,-id)
estccc=cccvc(long,"value","id","type")
estccc

kendall.global(predobs)
W<-kendall.global(predobs)
Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman") 

out<-as.data.frame(cbind(0,0))
for (i in 1:10000)
 { preds<-rnorm(30, mean=10, sd=2)
obs<-rnorm(30, mean=10, sd=2)
predobs<-as.data.frame(cbind(preds,obs))
xx<-as.vector(predobs$preds)
yy<-as.vector(predobs$obs)
sc.agr<-agreement(x=xx, y=yy, error="const",TDI_a = 1, target="random")
out[i,1]<-sc.agr$Estimate$CCC
out[i,2]<-sc.agr$Conf_Limit$CCC}
View(arrange(out, V2)) # lower confidence limit is zero at CCC=0.29


no126data$scaled.n15.bromeliad.final<-(no126data$n15.bromeliad.final+4)^0.125
fulldata$sqrt.decomp<-sqrt(fulldata$decomp)
fulldata$log.co2.final<-log(fulldata$co2.final)

concord.out<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())

singlesite.concord.hydro("filter.feeder_bio", fulldata, "puertorico", "argentina", 100, 2, "gaussian")
sites<-c("puertorico","argentina", "macae", "frenchguiana")
multisite.concord.hydro("filter.feeder_bio", fulldata, sites, "argentina", 100, 2, "gaussian")
multisite.concord("filter.feeder_bio", fulldata, sites, "argentina", 100, 2, "gaussian")

singlesite.concord<-function(a, dataset, reference, target, scalar, deltalimit, family) 
{
  b<-filter(dataset, site%in%reference)
  c<-filter(dataset, site%in%target)
  summ<-cbind(a, reference, target, "rain", 1 ,"NA", "NA", "NA", "NA", "NA", "NA")
  
  try({
    if (family=="nb") {
      correction<-mean(log(round(c[,a]*scalar)+0.001))-mean(log(round(b[,a]*scalar)+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.sitenb(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)]) #what is theta of best nb model? theta differs between nb models so can't average
      allsc<-aic.site(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else{
      correction<-mean(c[,a])-mean(b[,a])
      obs<-c[,a]   
      allsc<-aic.site((b[,a]),gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0); preds<- as.vector(sapply(sc.ests, predict, newdata = c))+correction#default to top model if only one
    try({
      sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#try because wont run if only one model within deltalimit
      preds<-as.vector(predict(sc.ests,c, full = TRUE))+correction
    })  
    sc.agr<-agreement(x=preds, y=obs, error="const",TDI_a = 1, target="random") #tdi is bogus, just to get CCC output 
    predobs<-as.matrix(cbind(preds,obs))
    W<-kendall.global(predobs)
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman")
    summ<-cbind(a, reference, target, "rain", 1 ,sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  
  return(summ)
}


singlesite.concord.hydro<-function(a, dataset, reference, target, scalar, deltalimit, family) 
{
  b<-filter(dataset, site%in%reference)
  c<-filter(dataset, site%in%target)
  summ<-cbind(a, reference, target, "hydro", 1 ,"NA", "NA", "NA", "NA", "NA", "NA")
  try({
    if (family=="nb") {
      correction<-mean(log(c[,a]*scalar+0.001))-mean(log(b[,a]*scalar+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.site.hydro.nb(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)]) #what is theta of best nb model? theta differs between nb models so can't average
      allsc<-aic.site.hydro(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else{
      correction<-mean(c[,a])-mean(b[,a])
      obs<-c[,a]   
      allsc<-aic.site.hydro(b[,a],gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0); preds<- as.vector(sapply(sc.ests, predict, newdata = c))+correction#default to top model if only one
    try({
      sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#try because wont run if only one model within deltalimit
      preds<-as.vector(predict(sc.ests,c, full = TRUE))+correction
    })  
    sc.agr<-agreement(x=preds, y=log(round(c[,a]*scalar)+0.001), error="const",TDI_a = 1, target="random") #tdi is bogus, just to get CCC output 
    predobs<-as.matrix(cbind(preds, log(round(c[,a]*scalar)+0.001)))
    W<-kendall.global(predobs)
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman")
    summ<-cbind(a, reference, target, "hydro", 1 ,sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  return(summ)
}


multisite.concord<-function(a, dataset, sites, target, scalar, deltalimit, family) #a is variable in quotes, target is site name in quotes
{
  b<-filter(dataset, site%in%sites)%>%filter(site%nin%target)
  c<-filter(dataset, site%in%target)
  nsites<-nrow(distinct(select(b, site)))
  predmatrix<-cbind(c["maxvol"], c["k.scalar"],c["mu.scalar"])
  predmatrix[,a]<-c[,a]
  nr<-nrow(c)
  predmatrix$site<-rep(b[,"site"][1], nr)
  summ<-cbind(a, deparse(substitute(dataset)), target, "rain", nsites ,"NA", "NA", "NA", "NA", "NA", "NA")
  try({
    if (family=="nb") 
    {
      correction<-mean(log(c[,a]*scalar+0.001))-mean(log(b[1:nr,a]*scalar+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.lmxnb.add(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
      allsc<-aic.lmx.add(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else
    {
      correction<-mean(c[,a])-mean(b[1:nr,a])
      obs<-c[,a]   
      allsc<-aic.lmx.add(b[,a],gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0)
    preds<- as.vector(sapply(sc.ests, predict, newdata = predmatrix))
    try({
      sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#now get average parameter estimates from plausible set
      preds<-as.vector(predict(sc.ests,predmatrix, full = TRUE))+correction
    }) 
    predobs<-as.matrix(cbind(preds, obs))
    sc.agr$Estimate$CCC[1]<-"NA"
    sc.agr$Estimate$Precision[1]<- cor(predobs, use="pairwise.complete.obs", method="pearson")[2]
    sc.agr$Estimate$Accuracy[1]<-"NA"
    try({sc.agr<-agreement(x=preds, y=obs, error="const",TDI_a = 1, target="random")
    })
    W<-kendall.global(predobs)
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman") 
    summ<-cbind(a, deparse(substitute(dataset)), target, "rain", nsites, sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  return(summ)
}


multisite.concord.hydro<-function(a, dataset, sites, target, scalar, deltalimit, family) #a is variable in quotes, target is site name in quotes
{
  b<-filter(dataset, site%in%sites)%>%filter(site%nin%target)
  c<-filter(dataset, site%in%target)
  nsites<-nrow(distinct(select(b, site)))
  predmatrix<-cbind(c["maxvol"], c["cv.depth"],c["prop.driedout.days"],c["prop.overflow.days"],c["mean.depth"], c["long_dry"], c["last_wet"])
  predmatrix[,a]<-c[,a]
  nr<-nrow(c)
  predmatrix$site<-rep(b[,"site"][1], nr)
  summ<-cbind(a, deparse(substitute(dataset)), target, "hydro", nsites ,"NA", "NA", "NA", "NA", "NA", "NA")
  try({
    if (family=="nb") 
    {
      correction<-mean(log(c[,a]*scalar+0.001))-mean(log(b[1:nr,a]*scalar+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.hydro.nb.add(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
      allsc<-aic.hydro.add(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else
    {
      correction<-mean(c[,a])-mean(b[1:nr,a])
      obs<-c[,a]   
      allsc<-aic.hydro.add(b[,a],gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0)
    preds<- as.vector(sapply(sc.ests, predict, newdata = predmatrix))+correction
    try({    
      sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#now get average parameter estimates from plausible set
      preds<-as.vector(predict(sc.ests,predmatrix, full = TRUE))+correction
    })
    predobs<-as.matrix(cbind(preds, obs))
    sc.agr$Estimate$CCC[1]<-"NA"
    sc.agr$Estimate$Precision[1]<- cor(predobs, use="pairwise.complete.obs", method="pearson")[2]
    sc.agr$Estimate$Accuracy[1]<-"NA"
    try({sc.agr<-agreement(x=preds, y=obs, error="const",TDI_a = 1, target="random")
    })
    W<-kendall.global(predobs)
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman") 
    summ<-cbind(a, deparse(substitute(dataset)), target, "hydro", nsites, sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  return(summ)
}

multicombo<-function(sites1, response, dataset, scalar, deltalimit, family){
  comb1<-combn(sites1, 1)
  comb2<-ncol(comb1)
  outr<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
  outh<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
  
  for (i in 1:comb2) 
  {
    outr[i,]<- multisite.concord(response, dataset, sites1, comb1[1,i],  scalar, deltalimit, family)
    outh[i,]<- multisite.concord.hydro(response, dataset, sites1, comb1[1,i],  scalar, deltalimit, family)
  }
  out<-rbind(outr, outh)
  return(out)
}


multimachine<-function(response, dataset, sites, scalar, deltalimit, family)
{
  for (f in 3:length(sites))
  {
    k<-nrow(concord.out)
    siter<-combn(sites,f)#must be higherthan 3
    for (g in 1:ncol(siter))
    {
      output<- multicombo(siter[,g], response, dataset, scalar, deltalimit, family)
      concord.out<-rbind(concord.out, output)
    }
  }
  return(concord.out)}


concord.magic<-function(sites, response, dataset, scalar, deltalimit, family)
{
  concord.out<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
  combin<-combn(sites, 2)
  combos<-ncol(combin)
  i<-0
  for (i in 1:combos) 
  {
    concord.out[i,]<- singlesite.concord(response, dataset, combin[1,i], combin[2,i], scalar, deltalimit, family)
  }
  j<-nrow(concord.out); i<-0
  for (i in 1:combos) 
  {
    concord.out[j+i,]<- singlesite.concord.hydro(response, dataset, combin[1,i], combin[2,i], scalar, deltalimit, family)
  }
  multiout<-multimachine(response, dataset, sites, scalar, deltalimit, family)
  concord.out<-rbind(concord.out, multiout)
  concord.out$Reference<-as.factor(concord.out$Reference)
  concord.out$Target<-as.factor(concord.out$Target)
  concord.out$Model<-as.factor(concord.out$Model)
  concord.out$Sites<-as.numeric(concord.out$Sites)
  concord.out$CCC<-as.numeric(concord.out$CCC)
  concord.out$Kendall<-as.numeric(concord.out$Kendall)
  concord.out$Spearman<-as.numeric(concord.out$Spearman)
  concord.out$Precision<-as.numeric(concord.out$Precision)
  concord.out$Accuracy<-as.numeric(concord.out$Accuracy)
  return(concord.out)
}
#

sites<-c("puertorico","argentina", "macae", "frenchguiana")
concord.out1<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
concord.out1<-concord.magic(sites, "filter.feeder_bio", nocadata, 100, 2, "nb")

sites<-c("puertorico","argentina", "macae", "frenchguiana", "costarica", "colombia") 
concord.out2<-concord.magic(sites, "scraper_bio", nocadata, 100, 2, "nb")#works fine now save 2 models

sites<-c("puertorico","argentina", "macae", "frenchguiana", "costarica", "colombia") 
concord.out3<-concord.magic(sites, "shredder_bio", nocadata, 100, 2, "nb")#works fine now save 2 models

sites<-c("puertorico", "macae", "frenchguiana", "costarica") 
concord.out4<-concord.magic(sites, "engulfer_bio", noargco13data, 10, 2, "nb")#all models work, but nost corr are neg

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina") 
concord.out5<-concord.magic(sites, "gatherer_bio", no67185data, 100, 2, "nb")#a handful of models fail, not too bad

sites<-c("macae", "frenchguiana", "argentina") 
concord.out6<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
concord.out6<-concord.magic(sites, "piercer_bio", nococrprdata, 100, 2, "nb")#all but one model works

sites<-c("puertorico", "macae", "frenchguiana", "costarica")
noargco123cleandata<-filter(noargco123data, bacteria.per.nl.final%nin%NA)#original fails as puertorico has NA values
concord.out7<-concord.magic(sites, "bacteria.per.nl.final", noargco123cleandata, 100, 2, "nb")#all models work!

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina") 
concord.out8<-concord.magic(sites, "totalbio", noargco123data, 10, 2, "nb")#a mess!

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina") 
no126data$sqrt.n15.bromeliad.final<-(no126data$n15.bromeliad.final+4)^0.5
concord.out9<-concord.magic(sites, "sqrt.n15.bromeliad.final", no126data, 10, 2, "gaussian")#this one was used in the end
concord.out9b<-concord.magic(sites, "scaled.n15.bromeliad.final", no126data, 10, 2, "gaussian")

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina")
concord.out11<-concord.magic(sites, "sqrt.decomp", fulldata, 10, 2, "gaussian")#rain models hardly ever converge

sites<-c("macae", "frenchguiana", "costarica", "colombia", "argentina")#many NA models
nocaprdata$log.co2.final<-log(nocaprdata$co2.final)
concord.out10<-concord.magic(sites, "log.co2.final", nocaprdata, 10, 2, "gaussian")#rain models hardly ever converge

concord.out.all<-rbind(concord.out1, concord.out2, concord.out3, concord.out4, concord.out5, concord.out6, concord.out7, concord.out9, concord.out9b, concord.out10)

write.csv(concord.out.all, "Data/concord_out_all.csv")


concord.out.all%>%filter(Sites==1)%>%select(Response, Reference, Target, Model, Kendall)%>%filter(Kendall>0.65)
#both scrapers and gatherers in colombbia are well predicted by other sites....?why?

concord.out.all%>%select(Response, Sites, Reference, Target, Model, Kendall)%>%filter(Kendall>0.65)
#gatherers..ar, cr, co all predicetd well, never macae or frencguiana
#shredders: costa rica and puertorico predicetd well

tapply(fulldata$maxvol, fulldata$site, min)
tapply(fulldata$maxvol, fulldata$site, max)
codata$maxvol
codata$gatherer_bio
plot(log(codata$gatherer_bio)~log(codata$maxvol))
plot(log(codata$scraper_bio)~log(codata$maxvol))
cor.test(prdata$gatherer_bio, prdata$maxvol, method="kendall")
cor.test(fgdata$scraper_bio, fgdata$maxvol, method="kendall")



#consider cohen.kappa(cbind(x,y))

mean.concord<-concord.out.all %>%
  group_by(Response,Model, Sites) %>%
  summarise(mCCC = mean(CCC, na.rm=TRUE), mSpearman=mean(Spearman, na.rm=TRUE), mKendall=mean(Kendall, na.rm=TRUE), mPearson=mean(Precision, na.rm=TRUE))

ggplot(mean.concord, aes(Sites, mSpearman, linetype=Model, colour=Response)) +
  geom_point(size = 3) +
  geom_smooth(span=1.2,se=FALSE)

ggplot(mean.concord, aes(Sites, mKendall, linetype=Model, colour=Response)) +
  geom_point(size = 3) +
  geom_smooth(span=1.2,se=FALSE)

ggplot(mean.concord, aes(Sites, mPearson, linetype=Model, colour=Response)) +
  geom_point(size = 3) +
  geom_smooth(span=1.2,se=FALSE)

ggplot(concord.out.all, aes(Sites, CCC, linetype=Model, colour=Response)) +
  geom_jitter(width=0.25, size = 1, alpha = 1/3) +
  geom_smooth(span=1.2, se=FALSE)

ggplot(concord.out.all, aes(Sites, Precision, linetype=Model, colour=Response)) +
  geom_jitter(width=0.25, size = 1, alpha = 1/3) +
  geom_smooth(span=1.2, se=FALSE)

ggplot(concord.out.all, aes(Sites, Accuracy, linetype=Model, colour=Response)) +
  geom_jitter(width=0.25, size = 1, alpha = 1/3) +
  geom_smooth(span=1.2, se=FALSE)

ggplot(concord.out.all, aes(Sites, Kendall, linetype=Model, colour=Response)) +library(cccrm)
  geom_jitter(width=0.25, size = 1, alpha = 1/3) +
  geom_smooth(span=1.2, se=TRUE, alpha = 1/10)

ggplot(concord.out.all, aes(Sites, Spearman, linetype=Model, colour=Response)) +
  geom_jitter(width=0.25, size = 1, alpha = 1/3) +
  geom_smooth(span=1.2, se=TRUE, alpha = 1/10)

#=======reduced set of hydrologically responsive sites
# these are ar and ma, and then pr and then fg

sites<-c("puertorico","argentina", "macae", "frenchguiana")
concord.outr1<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
concord.outr1<-concord.magic(sites, "filter.feeder_bio", fulldata, 100, 2, "nb")

sites<-c("puertorico","argentina", "macae", "frenchguiana") 
concord.outr2<-concord.magic(sites, "scraper_bio", nocadata, 100, 2, "nb")#works fine now save 2 models

sites<-c("puertorico","argentina", "macae", "frenchguiana") 
concord.outr3<-concord.magic(sites, "shredder_bio", nocadata, 100, 2, "nb")#works fine now save 2 models

sites<-c("puertorico", "macae", "frenchguiana") 
concord.outr4<-concord.magic(sites, "engulfer_bio", noargco13data, 10, 2, "nb")#all models work, but nost corr are neg

sites<-c("puertorico", "macae", "frenchguiana", "argentina") 
concord.outr5<-concord.magic(sites, "gatherer_bio", no67185data, 100, 2, "nb")#a handful of models fail, not too bad

sites<-c("macae", "frenchguiana", "argentina") 
concord.outr6<-concord.magic(sites, "piercer_bio", nococrprdata, 100, 2, "nb")#all but one model works

sites<-c("puertorico", "macae", "frenchguiana")
noargco123cleandata<-filter(noargco123data, bacteria.per.nl.final%nin%NA)#original fails as puertorico has NA values
concord.outr7<-concord.magic(sites, "bacteria.per.nl.final", noargco123cleandata, 100, 2, "nb")#all models work!

sites<-c("puertorico", "macae", "frenchguiana", "argentina") 
concord.outr8<-concord.magic(sites, "totalbio", noargco123data, 10, 2, "nb")#a mess!

sites<-c("puertorico", "macae", "frenchguiana",  "argentina") 
no126data$sqrt.n15.bromeliad.final<-(no126data$n15.bromeliad.final+4)^0.5
concord.outr9<-concord.magic(sites, "sqrt.n15.bromeliad.final", no126data, 10, 2, "gaussian")#this one was used in the end

sites<-c("puertorico", "macae", "frenchguiana", "argentina")
concord.outr11<-concord.magic(sites, "sqrt.decomp", fulldata, 10, 2, "gaussian")#rain models hardly ever converge

sites<-c("macae", "frenchguiana", "argentina")#many NA models
nocaprdata$log.co2.final<-log(nocaprdata$co2.final)
concord.outr10<-concord.magic(sites, "log.co2.final", nocaprdata, 10, 2, "gaussian")#rain models hardly ever converge

concordr.out.all<-rbind(concord.outr1, concord.outr2, concord.outr3, concord.outr4, concord.outr5, concord.outr6, concord.outr7, concord.outr9, concord.outr10)

mean.concordr<-concordr.out.all %>%
  group_by(Response,Model, Sites) %>%
  summarise(mCCC = mean(CCC, na.rm=TRUE), mSpearman=mean(Spearman, na.rm=TRUE), mKendall=mean(Kendall, na.rm=TRUE), mPearson=mean(Precision, na.rm=TRUE))

ggplot(mean.concordr, aes(Sites, mKendall, linetype=Model, colour=Response)) +
  geom_point(size = 3) +
  geom_smooth(span=1.2,se=FALSE)

## remove NA models from concordance

singlesite.concord2<-function(a, dataset, reference, target, scalar, deltalimit, family) 
{
  b<-filter(dataset, site%in%reference)
  c<-filter(dataset, site%in%target)
  summ<-cbind(a, reference, target, "rain", 1 ,"NA", "NA", "NA", "NA", "NA", "NA")
  
  try({
    if (family=="nb") {
      correction<-mean(log(round(c[,a]*scalar)+0.001))-mean(log(round(b[,a]*scalar)+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.sitenb(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)]) #what is theta of best nb model? theta differs between nb models so can't average
      allsc<-aic.site(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else{
      correction<-mean(c[,a])-mean(b[,a])
      obs<-c[,a]   
      allsc<-aic.site((b[,a]),gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0); preds<- as.vector(sapply(sc.ests, predict, newdata = c))+correction#default to top model if only one
    modlist<-c(rownames(subset(allsc, delta < 2)))
    print(modlist)
    if("m0"%in%modlist){preds<-as.vector(rep(NA,nrow(c)))} else
    {
      try({
        sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#try because wont run if only one model within deltalimit
        preds<-as.vector(predict(sc.ests,c, full = TRUE))+correction
      })  
    }
    sc.agr<-agreement(x=preds, y=obs, error="const",TDI_a = 1, target="random") #tdi is bogus, just to get CCC output 
    predobs<-as.matrix(cbind(preds,obs))
    if("m0"%in%modlist){W$Concordance_analysis[1]<-"NA"} else {W<-kendall.global(predobs)}
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman")
    summ<-cbind(a, reference, target, "rain", 1 ,sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  
  return(summ)
}



singlesite.concord.hydro2<-function(a, dataset, reference, target, scalar, deltalimit, family) 
{
  b<-filter(dataset, site%in%reference)
  c<-filter(dataset, site%in%target)
  summ<-cbind(a, reference, target, "hydro", 1 ,"NA", "NA", "NA", "NA", "NA", "NA")
  try({
    if (family=="nb") {
      correction<-mean(log(c[,a]*scalar+0.001))-mean(log(b[,a]*scalar+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.site.hydro.nb(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)]) #what is theta of best nb model? theta differs between nb models so can't average
      allsc<-aic.site.hydro(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else{
      correction<-mean(c[,a])-mean(b[,a])
      obs<-c[,a]   
      allsc<-aic.site.hydro(b[,a],gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0); preds<- as.vector(sapply(sc.ests, predict, newdata = c))+correction#default to top model if only one
    modlist<-c(rownames(subset(allsc, delta < 2)))
    print(modlist)
    if("m0"%in%modlist){preds<-as.vector(rep(NA,nrow(c)))} else
    {  
      try({
        sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#try because wont run if only one model within deltalimit
        preds<-as.vector(predict(sc.ests,c, full = TRUE))+correction
      }) 
    }
    sc.agr<-agreement(x=preds, y=log(round(c[,a]*scalar)+0.001), error="const",TDI_a = 1, target="random") #tdi is bogus, just to get CCC output 
    predobs<-as.matrix(cbind(preds, log(round(c[,a]*scalar)+0.001)))
    if("m0"%in%modlist){W$Concordance_analysis[1]<-"NA"} else {W<-kendall.global(predobs)}
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman")
    summ<-cbind(a, reference, target, "hydro", 1 ,sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  return(summ)
}

multisite.concord2<-function(a, dataset, sites, target, scalar, deltalimit, family) #a is variable in quotes, target is site name in quotes
{
  b<-filter(dataset, site%in%sites)%>%filter(site%nin%target)
  c<-filter(dataset, site%in%target)
  nsites<-nrow(distinct(select(b, site)))
  predmatrix<-cbind(c["maxvol"], c["k.scalar"],c["mu.scalar"])
  predmatrix[,a]<-c[,a]
  nr<-nrow(c)
  predmatrix$site<-rep(b[,"site"][1], nr)
  summ<-cbind(a, deparse(substitute(dataset)), target, "rain", nsites ,"NA", "NA", "NA", "NA", "NA", "NA")
  try({
    if (family=="nb") 
    {
      correction<-mean(log(c[,a]*scalar+0.001))-mean(log(b[1:nr,a]*scalar+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.lmxnb.add(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
      allsc<-aic.lmx.add(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else
    {
      correction<-mean(c[,a])-mean(b[1:nr,a])
      obs<-c[,a]   
      allsc<-aic.lmx.add(b[,a],gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0)
    preds<- as.vector(sapply(sc.ests, predict, newdata = predmatrix))
    modlist<-c(rownames(subset(allsc, delta < 2)))
    print(modlist)
    if("m0"%in%modlist){preds<-as.vector(rep(NA,nrow(predmatrix)))} else
    {
      try({
        sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#now get average parameter estimates from plausible set
        preds<-as.vector(predict(sc.ests,predmatrix, full = TRUE))+correction
      }) 
    }
    predobs<-as.matrix(cbind(preds, obs))
    sc.agr$Estimate$CCC[1]<-"NA"
    sc.agr$Estimate$Precision[1]<- cor(predobs, use="pairwise.complete.obs", method="pearson")[2]
    sc.agr$Estimate$Accuracy[1]<-"NA"
    try({sc.agr<-agreement(x=preds, y=obs, error="const",TDI_a = 1, target="random")
    })
    if("m0"%in%modlist){W$Concordance_analysis[1]<-"NA"} else {W<-kendall.global(predobs)}
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman") 
    summ<-cbind(a, deparse(substitute(dataset)), target, "rain", nsites, sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  return(summ)
}


multisite.concord.hydro2<-function(a, dataset, sites, target, scalar, deltalimit, family) #a is variable in quotes, target is site name in quotes
{
  b<-filter(dataset, site%in%sites)%>%filter(site%nin%target)
  c<-filter(dataset, site%in%target)
  nsites<-nrow(distinct(select(b, site)))
  predmatrix<-cbind(c["maxvol"], c["cv.depth"],c["prop.driedout.days"],c["prop.overflow.days"],c["mean.depth"], c["long_dry"], c["last_wet"])
  predmatrix[,a]<-c[,a]
  nr<-nrow(c)
  predmatrix$site<-rep(b[,"site"][1], nr)
  summ<-cbind(a, deparse(substitute(dataset)), target, "hydro", nsites ,"NA", "NA", "NA", "NA", "NA", "NA")
  try({
    if (family=="nb") 
    {
      correction<-mean(log(c[,a]*scalar+0.001))-mean(log(b[1:nr,a]*scalar+0.001))
      obs<-log(round(c[,a]*scalar)+0.001)
      nbset<-aic.hydro.nb.add(round(b[,a]*scalar), b)  #set of nb models
      newtheta<-as.numeric(as.vector(levels(nbset[1]$init.theta))[as.numeric(nbset[1]$init.theta)])
      allsc<-aic.hydro.add(round(b[,a]*scalar),negative.binomial(theta=newtheta), b) #rerun AIC with fixed theta
    }
    else
    {
      correction<-mean(c[,a])-mean(b[1:nr,a])
      obs<-c[,a]   
      allsc<-aic.hydro.add(b[,a],gaussian, b)
    }
    sc.ests<-get.models(allsc, subset= delta == 0)
    preds<- as.vector(sapply(sc.ests, predict, newdata = predmatrix))+correction
    modlist<-c(rownames(subset(allsc, delta < 2)))
    print(modlist)
    if("m0"%in%modlist){preds<-as.vector(rep(NA,nrow(predmatrix)))} else
    {
      try({    
        sc.ests<-model.avg(allsc, subset= delta < deltalimit, revised.var = TRUE, fit = TRUE)#now get average parameter estimates from plausible set
        preds<-as.vector(predict(sc.ests,predmatrix, full = TRUE))+correction
      })
    }
    predobs<-as.matrix(cbind(preds, obs))
    sc.agr$Estimate$CCC[1]<-"NA"
    sc.agr$Estimate$Precision[1]<- cor(predobs, use="pairwise.complete.obs", method="pearson")[2]
    sc.agr$Estimate$Accuracy[1]<-"NA"
    try({sc.agr<-agreement(x=preds, y=obs, error="const",TDI_a = 1, target="random")
    })
    if("m0"%in%modlist){W$Concordance_analysis[1]<-"NA"} else {W<-kendall.global(predobs)}
    Sp<-cor(predobs, use="pairwise.complete.obs", method="spearman") 
    summ<-cbind(a, deparse(substitute(dataset)), target, "hydro", nsites, sc.agr$Estimate$CCC[1], sc.agr$Estimate$Precision[1], sc.agr$Estimate$Accuracy[1], W$Concordance_analysis[1], Sp[2])
  })
  return(summ)
}

multicombo2<-function(sites1, response, dataset, scalar, deltalimit, family){
  comb1<-combn(sites1, 1)
  comb2<-ncol(comb1)
  outr<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
  outh<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
  
  for (i in 1:comb2) 
  {
    outr[i,]<- multisite.concord2(response, dataset, sites1, comb1[1,i],  scalar, deltalimit, family)
    outh[i,]<- multisite.concord.hydro2(response, dataset, sites1, comb1[1,i],  scalar, deltalimit, family)
  }
  out<-rbind(outr, outh)
  return(out)
}


multimachine2<-function(response, dataset, sites, scalar, deltalimit, family)
{
  for (f in 3:length(sites))
  {
    k<-nrow(concord.out)
    siter<-combn(sites,f)#must be higherthan 3
    for (g in 1:ncol(siter))
    {
      output<- multicombo2(siter[,g], response, dataset, scalar, deltalimit, family)
      concord.out<-rbind(concord.out, output)
    }
  }
  return(concord.out)}


concord.magic2<-function(sites, response, dataset, scalar, deltalimit, family)
{
  concord.out<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
  combin<-combn(sites, 2)
  combos<-ncol(combin)
  i<-0
  for (i in 1:combos) 
  {
    concord.out[i,]<- singlesite.concord2(response, dataset, combin[1,i], combin[2,i], scalar, deltalimit, family)
  }
  j<-nrow(concord.out); i<-0
  for (i in 1:combos) 
  {
    concord.out[j+i,]<- singlesite.concord.hydro2(response, dataset, combin[1,i], combin[2,i], scalar, deltalimit, family)
  }
  multiout<-multimachine2(response, dataset, sites, scalar, deltalimit, family)
  concord.out<-rbind(concord.out, multiout)
  concord.out$Reference<-as.factor(concord.out$Reference)
  concord.out$Target<-as.factor(concord.out$Target)
  concord.out$Model<-as.factor(concord.out$Model)
  concord.out$Sites<-as.numeric(concord.out$Sites)
  concord.out$CCC<-as.numeric(concord.out$CCC)
  concord.out$Kendall<-as.numeric(concord.out$Kendall)
  concord.out$Spearman<-as.numeric(concord.out$Spearman)
  concord.out$Precision<-as.numeric(concord.out$Precision)
  concord.out$Accuracy<-as.numeric(concord.out$Accuracy)
  return(concord.out)
}

singlesite.concord2("filter.feeder_bio", fulldata, "cardoso", "argentina", 100, 2, "gaussian")
singlesite.concord.hydro2("filter.feeder_bio", fulldata, "frenchguiana", "argentina", 100, 2, "gaussian")

singlesite.concord.hydro("gatherer_bio", fulldata, "colombia", "costarica", 100, 2, "gaussian")#0.71
singlesite.concord.hydro("gatherer_bio", fulldata, "costarica", "colombia", 100, 2, "gaussian")#0.79
singlesite.concord.hydro("gatherer_bio", fulldata, "costarica", "frenchguiana", 100, 2, "gaussian")#0.41
singlesite.concord.hydro("gatherer_bio", fulldata, "frenchguiana","costarica",  100, 2, "gaussian")#0.48
singlesite.concord.hydro("gatherer_bio", fulldata, "frenchguiana","puertorico", 100, 2, "gaussian")#0.50
singlesite.concord.hydro("gatherer_bio", fulldata, "puertorico","frenchguiana", 100, 2, "gaussian")#0.43
singlesite.concord.hydro("gatherer_bio", fulldata, "argentina","macae", 100, 2, "gaussian")#0.48
singlesite.concord.hydro("gatherer_bio", fulldata, "macae","argentina", 100, 2, "gaussian")#0.62

singlesite.concord.hydro("shredder_bio", fulldata, "colombia", "costarica", 100, 2, "gaussian")#0.57
singlesite.concord.hydro("shredder_bio", fulldata, "costarica", "colombia", 100, 2, "gaussian")#0.45
singlesite.concord.hydro("shredder_bio", fulldata, "costarica", "frenchguiana", 100, 2, "gaussian")#0.37
singlesite.concord.hydro("shredder_bio", fulldata, "frenchguiana","costarica",  100, 2, "gaussian")#0.50
singlesite.concord.hydro("shredder_bio", fulldata, "frenchguiana","puertorico", 100, 2, "gaussian")#0.49
singlesite.concord.hydro("shredder_bio", fulldata, "puertorico","frenchguiana", 100, 2, "gaussian")#0.37
singlesite.concord.hydro("shredder_bio", fulldata, "argentina","macae", 100, 2, "gaussian")#0.56
singlesite.concord.hydro("shredder_bio", fulldata, "macae","argentina", 100, 2, "gaussian")#0.44

sites<-c("cardoso", "argentina")
multisite.concord2("filter.feeder_bio", fulldata, sites, "macae", 100, 2, "gaussian") #a is variable in quotes, target is site name in quotes
sites<-c("puertorico","argentina", "frenchguiana")
multisite.concord2("filter.feeder_bio", fulldata, sites, "macae", 100, 2, "gaussian") #a is variable in quotes, target is site name in quotes
multisite.concord.hydro2("filter.feeder_bio", fulldata, sites, "macae", 100, 2, "gaussian") #a is variable in quotes, target is site name in quotes
multisite.concord.hydro2("engulfer_bio", fulldata, sites, "macae", 100, 2, "gaussian") #a is variable in quotes, target is site name in quotes

sites<-c("puertorico","argentina", "macae", "frenchguiana")
concord.out1<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
concord.out1<-concord.magic2(sites, "filter.feeder_bio", fulldata, 100, 2, "nb")

sites<-c("puertorico","argentina", "macae", "frenchguiana", "costarica", "colombia") 
concord.out2<-concord.magic2(sites, "scraper_bio", nocadata, 100, 2, "nb")#works fine now save 2 models

sites<-c("puertorico","argentina", "macae", "frenchguiana", "costarica", "colombia") 
concord.out3<-concord.magic2(sites, "shredder_bio", nocadata, 100, 2, "nb")#works fine now save 2 models

sites<-c("puertorico", "macae", "frenchguiana", "costarica") 
concord.out4<-concord.magic2(sites, "engulfer_bio", noargco13data, 10, 2, "nb")#all models work, but nost corr are neg

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina") 
concord.out5<-concord.magic2(sites, "gatherer_bio", no67185data, 100, 2, "nb")#a handful of models fail, not too bad

sites<-c("macae", "frenchguiana", "argentina") 
concord.out6<-data.frame(Response=numeric(), Reference= numeric(), Target =numeric(), Model=numeric(), Sites=numeric(), CCC=numeric(),Precision=numeric(),Accuracy=numeric(), Kendall=numeric(), Spearman=numeric())
concord.out6<-concord.magic2(sites, "piercer_bio", nococrprdata, 100, 2, "nb")#all but one model works

sites<-c("puertorico", "macae", "frenchguiana", "costarica")
noargco123cleandata<-filter(noargco123data, bacteria.per.nl.final%nin%NA)#original fails as puertorico has NA values
concord.out7<-concord.magic2(sites, "bacteria.per.nl.final", noargco123cleandata, 100, 2, "nb")#all models work!

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina") 
concord.out8<-concord.magic2(sites, "totalbio", noargco123data, 10, 2, "nb")#a mess!

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina") 
no126data$sqrt.n15.bromeliad.final<-(no126data$n15.bromeliad.final+4)^0.5
concord.out9<-concord.magic2(sites, "sqrt.n15.bromeliad.final", no126data, 10, 2, "gaussian")#this one was used in the end

sites<-c("puertorico", "macae", "frenchguiana", "costarica", "colombia", "argentina")
concord.out11<-concord.magic2(sites, "sqrt.decomp", fulldata, 10, 2, "gaussian")#rain models hardly ever converge

sites<-c("macae", "frenchguiana", "costarica", "colombia", "argentina")#many NA models
nocaprdata$log.co2.final<-log(nocaprdata$co2.final)
concord.out10<-concord.magic2(sites, "log.co2.final", nocaprdata, 10, 2, "gaussian")#rain models hardly ever converge

concord.out.all2<-rbind(concord.out1, concord.out2, concord.out3, concord.out4, concord.out5, concord.out6, concord.out7, concord.out9, concord.out10)

mean.concord2<-concord.out.all2 %>%
  group_by(Response,Model, Sites) %>%
  summarise(mCCC = mean(CCC, na.rm=TRUE), mSpearman=mean(Spearman, na.rm=TRUE), mKendall=mean(Kendall, na.rm=TRUE), mPearson=mean(Precision, na.rm=TRUE))

ggplot(mean.concord2, aes(Sites, mKendall, linetype=Model, colour=Response)) +
  geom_point(size = 3) +
  geom_smooth(span=1.2,se=FALSE)


