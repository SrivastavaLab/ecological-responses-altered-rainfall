library(ade4)
library(lme4)
library(lmerTest)
library(MASS)
library(vcd)
library(fitdistrplus)
library(car)
library(knitr)
library(visreg)
library(vegan)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(Hmisc)
library(quantreg)
library(MuMIn)
library(readr)
library(tibble)
library(Agreement)

fulldata<-read.csv("../00_BWGrainfall_data/Data/BWG_wide_functional_groups_ibuttons.csv")

sim<-read.csv("Data/simblank.csv", row.names=1)

pitfam<-read.csv("Data/CostaRica_overflow_families.csv")

#===creation of custom datasets===================

fulldata$site<-as.factor(fulldata$site)
fulldata$mu.scalar.2<-fulldata$mu.scalar^2
fulldata$k.scalar.2<-fulldata$k.scalar^2
fulldata$intended.mu.2<-fulldata$intended.mu^2
fulldata$intended.k.2<-fulldata$intended.k^2
fulldata$change.k<-abs(log(fulldata$k.scalar))
fulldata$change.mu<-abs(log(fulldata$mu.scalar))

fulldata$bacteria.per.nl.final<-(fulldata$bacteria.per.ml.final/1000000)

# predator_biomass is the sum of predator_engulfer_total_biomass and
# predator_piercer_total_biomass. Check:
with(fulldata,
     stopifnot(all.equal(predator_biomass, 
                         predator_engulfer_total_biomass + 
                           predator_piercer_total_biomass)))

# prey_biomass is the sum of prey_scraper_total_biomass,
# prey_gatherer_total_biomass, prey_shredder_total_biomass and
# prey_piercer_total_biomass and prey_filter.feeder_total_biomass
with(fulldata,
     # Stop everything if prey_biomass does not equal the sum of these columns
     stopifnot(all.equal(prey_biomass, 
                         prey_scraper_total_biomass + 
                           prey_gatherer_total_biomass +
                           prey_shredder_total_biomass +
                           prey_piercer_total_biomass + 
                           prey_filter.feeder_total_biomass)))


# detritivore_total_biomass is a separate thing from total_prey_biomass. This is
# ONLY the scrapers, gatherers and shredders
fulldata$detritivore_total_biomass <- with(fulldata, prey_scraper_total_biomass + 
                                             prey_gatherer_total_biomass +
                                             prey_shredder_total_biomass)


fulldata$PDratio <- with(fulldata, predator_biomass / detritivore_total_biomass)

# This next line seems unnecessary now -- we have this number as prey_biomass.
# This includes also prey which are piercing (Cecidomyiidae) 

# fulldata$totalbio <- fulldata$predator_bio + fulldata$detritivore_bio + fulldata$filter.feeder_bio

# fulldata <- 

full_data_temp_trophic <- fulldata %>% 
  # Separating taxonomic groups by their trophic position
  mutate(detChiron_bio           = Chironomidae_bio - Tanypodinae_bio,
         predCerato_bio          = Bezzia_bio + Sphaeromias_bio + Stilobezzia_bio + Culicoides_bio,
         detCerato_bio           = Ceratopogonidae_bio - predCerato_bio,
         TipulidaeLimoniidae_bio = Tipulidae_bio + Limoniidae_bio,
         filtCulicidae           = Culex_bio + Wyeomyia_bio + Anopheles_bio) %>% 
  # deriving environmental data
  mutate(exposure = prop.driedout.days * mean_temp,
         # How does temperature vary from a the site average? This is equivalent
         # to subtracting the site mean from mean_temp
         change_mean_temp = resid(glm(mean_temp~site, family=gaussian, data=fulldata, na.action=na.exclude)),
         change_cv_temp = resid(glm(cv_mean_temp~site, family=gaussian, data=fulldata, na.action=na.exclude)))

#damn, once again we have lost the cardoso temperature data...why does this keep happening? 
#and the NA buttons in FG are registering as mean_temp = 0...need to change this too!

# dropping NAs in mean_temp would drop.. 30 bromeliads?! yup, this would be cardoso...instead we want to filter the 16 
# nonfuntional ibuttons from FG, and a few nonfunctional ones in cardoso

full_data_temp_trophic %>% 
  filter(is.na(mean_temp))

full_data_temp_trophic %>% 
  group_by(site) %>% 
  count


# renaming fulldata to match downstream analysis ---------------------------------------------------

## from here, rename things just to permit the script to run!
fulldata <- full_data_temp_trophic %>% ungroup %>% as.data.frame

# creating new columns that match names to old ones
fulldata <- fulldata %>% 
  mutate(shredder_bio       = prey_shredder_total_biomass,
         filter.feeder_bio  = prey_filter.feeder_total_biomass,
         scraper_bio        = prey_scraper_total_biomass,
         gatherer_bio       = prey_gatherer_total_biomass,
         engulfer_bio       = predator_engulfer_total_biomass,
         piercer_bio        = predator_piercer_total_biomass,
         totalbio           = prey_biomass + predator_biomass) # Is that the right equation?? does it leave out NA (ie animals whose trophic role is unknown)
#answer to Andrew's question - true, but very few unknown trophic, might be better to just sum everything eventually
#DS orrected error: piercer_bio != prey_piercer_total_biomass

#making a water chemistry axis

waterchem<-select(fulldata, oxygen.percent.final, ph.final, turbidity.final)%>%
  filter(oxygen.percent.final%nin%NA)%>%
  filter(ph.final%nin%NA)%>%
  filter(turbidity.final%nin%NA)
pcawc<- dudi.pca(waterchem, center = TRUE, scale = TRUE, scan = FALSE)
waterchem.df<-fulldata%>%
  filter(oxygen.percent.final%nin%NA)%>%
  filter(ph.final%nin%NA)%>%
  filter(turbidity.final%nin%NA)%>%
  select(oxygen.percent.final, ph.final, turbidity.final, k.scalar, mu.scalar, site, maxvol)%>%
  mutate(wc1=pcawc$li$Axis1, wc2=pcawc$li$Axis2)

#axis 1 58% inertia; Axis 2 33% inertia 
#                         CS1        CS2
#oxygen.percent.final  0.7035388 0.03277309
#ph.final             -0.1990401 0.96804271
#turbidity.final       0.6822142 0.24863475

fulldata$change_ph<-resid(glm(ph.final~site, family=gaussian, data=fulldata, na.action=na.exclude))

# Making smaller customized datasets in terms of filters ----------------------------------------------------



temptruedata<-filter(fulldata, mean_temp%nin%NA)
temphydrotruedata<-filter(temptruedata, cv.depth%nin%NA)
noleakydatatemp<-filter(temptruedata, site_brom.id%nin%c("macae_B24", "macae_B22", "macae_B9", "macae_B2", "macae_B11", "macae_B41", "argentina_15"))
noleakydata<-filter(fulldata, site_brom.id%nin%c("macae_B24", "macae_B22", "macae_B9", "macae_B2", "macae_B11", "macae_B41", "argentina_15"))


nocadata<-subset(fulldata,site!="cardoso")
nocadata$resid.driedout<-resid(glm(sqrt(prop.driedout.days)~maxvol+site*log(mu.scalar)*log(k.scalar), family=gaussian, data=nocadata))
nocadata$resid.overflow<-resid(glm(sqrt(prop.overflow.days)~maxvol+site*log(mu.scalar)*log(k.scalar), family=gaussian, data=nocadata))
nocadata$resid.cvdepth<-resid(glm(log(cv.depth)~maxvol+site*log(mu.scalar)*log(k.scalar), family=gaussian, data=nocadata))
nocadatatemp<-filter(nocadata, mean_temp%nin%NA)

ardata<-subset(fulldata,site=="argentina")
cadata<-subset(fulldata,site=="cardoso")
codata<-subset(fulldata,site=="colombia")
crdata<-subset(fulldata,site=="costarica")
fgdata<-subset(fulldata,site=="frenchguiana")
madata<-subset(fulldata,site=="macae")
prdata<-subset(fulldata,site=="puertorico")


noargdata<-subset(fulldata,site!="argentina")
noargmacdata<-subset(noargdata,site!="macae")
noargcodata<-subset(noargdata,site!="colombia" )
noargcadata<-subset(noargdata,site!="cardoso" )
noargprdata<-subset(noargdata,site!="puertorico" )
noargcaprdata<-subset(noargcadata,site!="puertorico" )
noargcocrdata<-subset(noargcodata,site!="costarica" )
nocodata<-subset(fulldata,site!="colombia")
nocoprdata<-subset(nocodata,site!="puertorico")
nofgdata<-subset(fulldata,site!="frenchguiana")
nocafgdata<-subset(nocadata,site!="frenchguiana")
nocacodata<-subset(nocadata,site!="colombia")
nocaprdata<-subset(nocadata,site!="puertorico")
nocacoprdata<-subset(nocacodata,site!="puertorico")
nococrdata<-subset(nocodata,site!="costarica")
nococrprdata<-subset(nococrdata,site!="puertorico")
nocacocrdata<-subset(nocacodata,site!="costarica")
noprdata<-subset(fulldata,site!="puertorico")
nomadata<-subset(fulldata,site!="macae")
noargcacodata<-subset(nocacodata,site!="argentina" )

camadata<-rbind(cadata, madata)
arcodata<-rbind(ardata, codata)
cafgmadata<-rbind(cadata, fgdata, madata)
camaprdata<-rbind(cadata, madata, prdata)
argcacodata<-rbind(ardata, cadata, codata)
cocrdata<-rbind(codata, crdata)
argcrdata<-rbind(ardata, crdata)
cafgdata<-rbind(cadata, fgdata)
cacrmadata<-rbind(cadata, crdata, madata)
maprdata<-rbind(madata, prdata)
fgmadata<-rbind(fgdata, madata)
fgmaprdata<-rbind(fgdata, madata, prdata)

nococr140data<-nococrdata[-140,]#working, corrects filter feeders
no67185data<-fulldata[-(c(67,185)),] #working
noargco13data<-noargcodata[-c(13, 127, 129),]
no126data<-fulldata[-126,]#working
noargco123data<-noargcodata[-123,]
noargco123datatemp<-filter(noargco123data, mean_temp%nin%NA)

nocaco170173<-nocacodata[-c(113,110),]
target<-c("67","61","66","68","90","79","61")
nocafgdata<-rownames_to_column(nocafgdata, "row")
nocafgclean<-filter(nocafgdata, row%nin%target)


