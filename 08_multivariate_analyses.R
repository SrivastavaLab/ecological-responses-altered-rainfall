source("01_datasets_for_paper1.R")
source("02_functions.R")

#=====multivariate analysis-------------------

#+ multivariate

fngrp<-select(fulldata, engulfer_bio, filter.feeder_bio, gatherer_bio, piercer_bio,scraper_bio,shredder_bio)
brom<-select(fulldata, site, change.mu, change.k, mu.scalar, k.scalar, intended.mu, intended.k, leaf.number,mean.diam,catchment.area, maxvol, cv.depth, prop.overflow.days,prop.driedout.days)
brom.min<-select(brom, site, mu.scalar, k.scalar, maxvol)
brom.med<-select(brom, site, mu.scalar, k.scalar, intended.mu, intended.k, maxvol)
brom.temptrue<-select(temphydrotruedata, site, change.mu, change.k, mu.scalar, k.scalar, intended.mu, intended.k, leaf.number,mean.diam,catchment.area, maxvol)
brom.noca<-select(nocadata, site, change.mu, change.k, mu.scalar, k.scalar, intended.mu, intended.k, leaf.number,mean.diam,catchment.area, maxvol)

hydro<-select(temphydrotruedata, cv.depth, prop.overflow.days,prop.driedout.days, mean.depth,long_dry,last_wet, change_cv_temp,change_mean_temp)
purehydro<-select(nocadata, cv.depth, prop.overflow.days,prop.driedout.days, mean.depth,long_dry,last_wet)
hydro.names<-names(purehydro)

families<-select(fulldata, Calamoceratidae_bio,Candonidae_bio,Cecidomyiidae_bio,Ceratopogonidae_bio,Chironomidae_bio,
                 Coenagrionidae_bio,Corethrellidae_bio,Culicidae_bio, Dytiscidae_bio,
                 Empididae_bio,Enchytraeoidae_bio,Ephydridae_bio,Hirudinea_bio,Hydrophilidae_bio,Lampyridae_bio,Limnocytheridae_bio, 
                 Naididae_bio, Pseudostigmatidae_bio,Psychodidae_bio,Ptilodactylidae_bio,
                 Scirtidae_bio,Stratiomyidae_bio,Syrphidae_bio,Tabanidae_bio,TipulidaeLimoniidae_bio)
#I am not including Curculionidae_bio...leafminers
apply(families, 2, datacheck)
#Excluded families with 5 or fewer bromeliads: Aeolosomatidae_bio (1), Anisopodidae_bio (5), Dolichopodidae_bio (5),
#Elateridae_bio (4),  Elmidae_bio (1), Muscidae_bio (4), Periscelididae_bio (2), Phoridae_bio (2), 
#Scatopsidae_bio (4), Sciaridae_bio (1), Sphaeroceridae_bio (5)

fn.names<-names(fngrp)
par(mfrow = c(1, 1))

fn.dist<-vegdist(fngrp)
fn.mds1<-metaMDS(fn.dist)
ordiplot(fn.mds1)
fnenv1<-envfit(fn.mds1, brom, permu = 999, na.rm=TRUE); fnenv1
ordiplot(fn.mds1); plot(fnenv1)
fnenv2<-envfit(fn.mds1, brom.med, permu = 999, na.rm=TRUE); fnenv2
ordiplot(fn.mds1); plot(fnenv2)
fnenv3<-envfit(fn.mds1, brom.min, permu = 999, na.rm=TRUE); fnenv3
ordiplot(fn.mds1); with(brom.min, ordiellipse(fn.mds1, site, kind = "se", conf = 0.95));plot(fnenv3)

fam.dist<-vegdist(families)
fam.mds1<-metaMDS(fam.dist)
ordiplot(fam.mds1)
famenv1<-envfit(fam.mds1, brom, permu = 999, na.rm=TRUE); famenv1
ordiplot(fam.mds1); plot(famenv1)
famenv2<-envfit(fam.mds1, brom.med, permu = 999, na.rm=TRUE); famenv2
ordiplot(fam.mds1); plot(famenv2)

hydro.dist<-vegdist(purehydro, method="euclidean")

famenv3<-envfit(fam.mds1, brom.min, permu = 999, na.rm=TRUE); famenv3
ordiplot(fam.mds1); with(brom.min, ordiellipse(fam.mds1, site, kind = "se", conf = 0.95)); plot(famenv3)

protest(fn.mds1, fam.mds1) #r = 0.48, procrustes sumsq = 0.765, p =0.001


#the following are anovas based on CCAs. These do not consider interactions between site and variables
#so just looking at general effects
fn.cca1 <- cca(fngrp ~ site+mu.scalar+k.scalar+maxvol+intended.k+intended.mu, brom); fn.cca1
fn.cca1b <- cca(fngrp ~ site+log(mu.scalar)+log(k.scalar)+maxvol+log(intended.k)+log(intended.mu), brom)
plot(fn.cca1); with(brom, ordiellipse(fn.cca1, site, kind = "se", conf = 0.95))
anova(fn.cca1, by = "term", step=200) #just site sig with type 1 ss
anova(fn.cca1, by = "margin", step=200)#site and k.scalar sig with type 3 ss
anova(fn.cca1b, by = "margin", step=200)#site, logk, log intendedk
plot(fn.cca1b)
orditorp(fn.cca1b,dis="sp",labels=fn.names)
plot(fn.cca1b, dis="sp")#mu and k cause shift from piercers to scrapers...

fam.cca1 <- cca(families ~ site+mu.scalar+k.scalar+maxvol+intended.k+intended.mu, brom); fam.cca1
plot(fam.cca1); with(brom, ordiellipse(fam.cca1, site, kind = "se", conf = 0.95))
anova(fam.cca1, by = "term", step=200) #site and intended mu sig with type 1 ss
anova(fam.cca1, by = "margin", step=200)#site, mu.scalar and intended.mu sig with type 3 ss
# so one story is that fn composition follows k, whereas taxonomy follows mu

#The following are CCAs with interactions, these cannot use type 3 (=marginal) tests
fn.cca0 <- cca(fngrp ~ 1, brom);fn.cca0
fn.cca0a <- cca(fngrp ~ log(maxvol)+site, brom);fn.cca0a #explains 44.85% of interia; 0.7588 of 1.3759 unexplained
fn.cca0b <- cca(fngrp ~ log(maxvol)+site+log(mu.scalar)+log(k.scalar)+change.k+change.mu+log(mu.scalar):log(k.scalar)+change.mu:change.k, brom); fn.cca0b
#without interactions, 46.69% of inertia explained, 0.7335 of 1.3759 unexplained
fn.cca2 <- cca(fngrp ~ mu.scalar+k.scalar+maxvol+site+site:k.scalar+site:mu.scalar+mu.scalar:k.scalar, brom); fn.cca2
fn.cca2b <- cca(fngrp ~ log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), brom); fn.cca2b
fn.cca2d <- cca(fngrp ~ log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar), brom); fn.cca2d
fn.cca2e <- cca(fngrp ~ log(mu.scalar)+log(k.scalar)+change.k+change.mu+log(maxvol)+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, brom); fn.cca2e
#so 0.5392 inertia unexplained...rainfall explains (0.7588-0.5392)/(0.7588) = 28.9% of non-site effect; of this rainfall effect (0.7335-0.5392)/(0.7588-0.5392)=88% is explained by contingency
plot(fn.cca2b); with(brom, ordiellipse(fn.cca2b, site, kind = "se", conf = 0.95))
anova(fn.cca2, by = "term", step=200)#only site and maxvol with type 1, type 3 ignores main effects
anova(fn.cca2b, by = "term", step=1000)#site, maxvol, logk p =0.045 (this is jumping around 0.05)
anova(fn.cca2d, by = "term", step=1000)#site, maxvol, log k p = 0.048 (this is jumping around 0.05)
anova(fn.cca2e, by = "term", step=1000)#site, maxvol, absolute change in k from 1 x site, explains 60% of inertia

ordistep(fn.cca2e, scope = formula(fn.cca2e), direction="both", perm.max = 200)#yields fn.cca2g

fn.cca2c <- cca(fngrp ~ log(k.scalar)+maxvol+site, brom); fn.cca2c
anova(fn.cca2c, by = "term", step=200)
plot(fn.cca2c); with(brom, ordiellipse(fn.cca2c, site, kind = "se", conf = 0.95)); text(fn.cca2c,dis="sp",labels=fn.names)
fn.cca2f <- cca(fngrp ~ log(maxvol)+site+site:change.k, brom); anova(fn.cca2f, by = "term", step=1000)#simplified model, all very sig
fn.cca2g <- cca(fngrp ~ site*change.k, brom); anova(fn.cca2g, by = "term", step=1000)#simplified model, all very sig


#saveRDS of fn.cca2f

require(readr)

figure06_cca_list <- list(
  fn.cca2g = fn.cca2g
)

saveRDS(figure06_cca_list, "figure/data/figure06_cca_list.rds")


#===functional group best plot
png(file = "FunctionalComp.png", units = "px", width = 7780, height=7380, res=1440)
par(mfrow=c(1,1))
plot(fn.cca2f); with(brom, ordiellipse(fn.cca2f, site, kind = "se", conf = 0.95)); text(fn.cca2f,dis="sp",labels=fn.names)

scl <- 3 ## scaling = 3
colvec <- c("red2", "green4", "mediumblue", "tan4", "plum4", "lawngreen", "cyan")
plot(fn.cca2f, type="n", scaling = scl)
with(brom, points(fn.cca2f, display = "sites", col = colvec[site],
                  scaling = scl, pch = 21, bg = colvec[site]))
points(fn.cca2f, display = "bp", scaling=scl)
with(brom, text(fn.cca2f,dis="sp",labels=fn.names, cex=0.5))
with(brom, legend("bottomleft", legend = levels(site), bty = "n",cex=0.5,
                  col = colvec, pch = 21, pt.bg = colvec))
dev.off()

#hydrology
brom.noca<-filter(brom, site%nin%"cardoso")
hydro.cca0 <- cca(purehydro ~ log(maxvol)+log(mu.scalar)+log(k.scalar)+change.k+change.mu+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, brom.noca) 
anova(hydro.cca0, by = "term", step=1000)#
hydro.cca1 <- cca(purehydro ~ log(maxvol)+log(mu.scalar)+log(k.scalar)+change.k+change.mu+site+site:log(k.scalar)+site:log(mu.scalar)+site:change.k+site:change.mu, brom.noca) 
anova(hydro.cca1, by = "term", step=1000)#

png(file = "HydroComp.png", units = "px", width = 7780, height=7380, res=1440)
par(mfrow=c(1,1))
plot(hydro.cca1); with(brom.noca, ordiellipse(hydro.cca1, site, kind = "se", conf = 0.95)); text(hydro.cca1,dis="sp",labels=hydro.names)

scl <- 3 ## scaling = 3
colvec <- c("red2", "green4", "mediumblue", "tan4", "plum4", "lawngreen", "cyan")
plot(hydro.cca1, type="n", scaling = scl)
with(brom.noca, points(hydro.cca1, display = "sites", col = colvec[site],
                       scaling = scl, pch = 21, bg = colvec[site]))
points(hydro.cca1, display = "bp", scaling=scl)
with(brom.noca, text(hydro.cca1,dis="sp",labels=hydro.names, cex=0.5))
with(brom.noca, legend("topleft", legend = levels(site), bty = "n",cex=0.5,
                       col = colvec, pch = 21, pt.bg = colvec))
dev.off()



fam.cca2 <- cca(families ~ mu.scalar+k.scalar+maxvol+site+site:k.scalar+site:mu.scalar+mu.scalar:k.scalar, brom); fam.cca2
fam.cca2b <- cca(families ~ log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), brom); fam.cca2b
plot(fam.cca2b); with(brom, ordiellipse(fam.cca2, site, kind = "se", conf = 0.95))
anova(fam.cca2, by = "term", step=200)#maxvol, site, mu.scalar x site
anova(fam.cca2b, by = "term", step=200)#maxvol, site, logmu.scalar x site
fam.cca2c <- cca(families ~ maxvol+site+log(mu.scalar):site, brom)
anova(fam.cca2c, by = "term", step=200)
fam.cca2d <- cca(families ~ log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar), brom); fam.cca2d
anova(fam.cca2d, by = "term", step=200)#maxvol, site, logmu.scalar x site
fam.cca2e <- cca(families ~ log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, brom); fam.cca2e
anova(fam.cca2e, by = "term", step=200)
fam.cca2f <- cca(families ~ maxvol+site+site:change.k+site:log(mu.scalar), brom);
anova(fam.cca2f, by = "term", step=200)#all sig

scl <- 3 ## scaling = 3
colvec <- c("red2", "green4", "mediumblue", "tan4", "plum4", "lawngreen", "cyan")
plot(fam.cca2c, type="n", scaling = scl)
with(brom, points(fam.cca2c, display = "sites", col = colvec[site],
                  scaling = scl, pch = 21, bg = colvec[site]))
points(fam.cca2c, display = "bp", scaling=scl)
with(brom, legend("bottomleft", legend = levels(site), bty = "n",cex=0.5,
                  col = colvec, pch = 21, pt.bg = colvec))

scl <- 3 ## scaling = 3
colvec <- c("red2", "green4", "mediumblue", "tan4", "plum4", "lawngreen", "cyan")
plot(fam.cca2f, type="n", scaling = scl)
with(brom, points(fam.cca2f, display = "sites", col = colvec[site],
                  scaling = scl, pch = 21, bg = colvec[site]))
points(fam.cca2f, display = "bp", scaling=scl)
with(brom, legend("bottomleft", legend = levels(site), bty = "n",cex=0.5,
                  col = colvec, pch = 21, pt.bg = colvec))

demo<-plot(fam.cca2c, display=c("cn","bp", "sites"), cex=0.5); with(brom, ordiellipse(fam.cca2c, site, kind = "se", conf = 0.95)); orditorp(fam.cca2c,dis="sp")
#identify(demo, "biplot")

#this loosely agrees with adonis that follows

#use adonis2 for marginal term testing but tosses out main effects, adonis uses sequential

fn.dist<-vegdist(fngrp, method="horn")

fn.adon1b<-adonis(fn.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", brom, perm = 200)
fn.adon1b #maxvol, site, log k is 0.02

fn.adon1c<-adonis(fn.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar), by="margin", brom, perm = 200)
fn.adon1c#maxvol, site, log k is 0.054

fn.adon1d<-adonis(fn.dist~log(k.scalar)+maxvol+site, by="margin", brom, perm = 200)
fn.adon1d#maxvol, site, logk is p =0.049, #so, like the cca, the effect of log k is bouncing around sig

fn.adon1e<-adonis(fn.dist~log(mu.scalar)+log(k.scalar)+log(maxvol)+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+change.k+change.mu+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", brom, perm = 200)
fn.adon1e#maxvol, site, log k is 0.05, oddly change.k not sig here unlike cca

fn.adon1e2<-adonis(fn.dist~log(k.scalar)+log(maxvol)+site+site:log(k.scalar)+change.k+site:change.k, by="margin", brom, perm = 1000)
fn.adon1e2#site and log k, (site x change k is at p=0.11)

fn.adon1f<-adonis(fn.dist~log(maxvol)+site, by="margin", brom, perm = 200)
fn.adon1f #so 29.79 of the 52.43 sunsq are unexplained, model fn.adon1e left unexplained 23.284, so (29.79-23.284)/ 29.79 = 21.8%

fn.adon1g<-adonis(fn.dist~log(maxvol)+site+log(mu.scalar)+log(k.scalar)+log(mu.scalar):log(k.scalar)+change.k+change.mu+change.mu:change.k, by="margin", brom, perm = 200)
fn.adon1g #here 28.875 unexplained, so general effects of rainfall explain (29.79 - 28.875)/(29.79-23.284) = 14%, and contingent explain (28.875-23.284)/(29.79-23.284) = 86% of this diff

fam.adon1<-adonis(fam.dist~mu.scalar+k.scalar+maxvol+site+site:k.scalar+site:mu.scalar+mu.scalar:k.scalar, by="margin", brom, perm = 200)
fam.adon1#maxvol, site, k*site. mu*site

fam.adon1b<-adonis(fam.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", brom, perm = 200)
fam.adon1b#maxvol, site, logk, logmu*site (logk*site is 0.055)

fam.adon1c<-adonis(fam.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar), by="margin", brom, perm = 200)
fam.adon1c#maxvol, site, logk, logmu*site (logk*site is 0.055) use this model

fam.adon2<-adonis(fam.dist~mu.scalar+k.scalar+maxvol, strata = brom$site, brom, by="margin", perm = 200)
fam.adon2#confirming above, if we restrict permutations to be within site, there is then no overall effect of mu or k (just maxvol)

#so another story is that contingency in taxonomic composition response doesnt extend to functional composition
#both adonis and the cca suggest that k marginally affects fn composition, and both mu and k additionally affect family comp but mu effects are contingent.

mantel(fn.dist, fam.dist, na.rm=TRUE) #r=0.86, p = 0.001

#hydrology multivariate

hydro.cca <- cca(hydro ~ log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, brom.temptrue)
anova(hydro.cca, by = "term", step=200) #overwhelmingly log k (not even site!)

hydro.cca2 <- cca(hydro ~ log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+site:change.k+site:change.mu+intended.mu+intended.k, brom.temptrue)
anova(hydro.cca2, by = "term", step=200) #overwhelmingly log k (not even site!)


#===hydrology on composition
fngrp2<-select(fulldata, site.x, engulfer_bio, filter.feeder_bio, gatherer_bio, piercer_bio,scraper_bio,shredder_bio)
bromnoca<-filter(brom, site!="cardoso")%>%select(-leaf.number)
fngrpnoca<-filter(fngrp2, site.x!="cardoso")%>%select(-site.x)
fn.cca2e <- cca(fngrpnoca ~ maxvol+site+log(mu.scalar)+log(k.scalar)+prop.driedout.days+prop.overflow.days+cv.depth, bromnoca); fn.cca2e
anova(fn.cca2e, by = "term", step=200)#amazing, no hydrology sig just maxvol, site, marg k (p=0.08)
plot(fn.cca2e);text(fn.cca2e,dis="sp",labels=fn.names)
#note that mu and k and overflow go in one direction, and driedout and cvdepth in the other, as expected
#here, increasing mu and k associated with scrapers, and away from piercers

#====separate fn grps===multivariate==
namedata<-read.csv( "C:/Users/Diane/Dropbox/BWG Drought experiment/Paper 1_thresholds/Data/bwgdb_names.csv")

eng_family<-as.data.frame(namedata$family[namedata$functional_group=="engulfer"])
eng_fam<-as.factor(unique(eng_family[,1]))
eng_bio<-select(fulldata,  Dytiscidae_bio, Corethrellidae_bio,
                Hydrophilidae_bio, Coenagrionidae_bio, Pseudostigmatidae_bio, Hirudinea_bio,Toxorhynchites_bio)

#*Axymyiidae removed as did not occur in drought exp,
#Dolichopodidae_bio, Periscelididae_bio removed as in 5 or fewer bromeliads
#Glossiphoniidae is removed as rarely identified as a family of Hirudinea (subclass)
#Muscidae removed as 6 records cannot confidently be attributed to engulfing pred

gath_family<-as.data.frame(namedata$family[namedata$functional_group=="gatherer"])
gath_fam<-as.factor(unique(gath_family[,1]))
gath_bio<-select(fulldata,  detCerato_bio,Ephydridae_bio,
                 Psychodidae_bio,  Stratiomyidae_bio, Syrphidae_bio, Enchytraeoidae_bio,
                 Naididae_bio, detChiron_bio)

#removed aff. Drosophilidae, Aulacigastridae (honduras), Lumbricidae, Entomobryidae as not in drought exp, and latter terrestrial
#removed Anisopodidae_bio,Sciaridae_bio, Phoridae_bio, Scatopsidae_bio, Sphaeroceridae_bio,
#consider using Oligochaeta_bio column (subclass) instead of Enchytraeoidae and Naididae as families within (as most not id to families)

pier_family<-as.data.frame(namedata$family[namedata$functional_group=="piercer"])
pier_fam<-as.factor(unique(pier_family[,1]))
pier_bio<-select(fulldata,  predCerato_bio, Empididae_bio,  Lampyridae_bio, Tabanidae_bio, Tanypodinae_bio)

#Cecidomyiidae_bio was removed at source as not predator
#Veliidae/Vellidae, Gerridae, Mesoveliidae, Enicocephalidae (dominica) removed as not in drought exp,
#changed chironomidae to tanypodinae
#Hydrophilidae_bio removed as piercer, as most classified as engulfers in database

#need to add Bezzia/Sphaeromias/Stilobezzia/Culicoides (genus) columns for pred ceratos.

filt_family<-as.data.frame(namedata$family[namedata$functional_group=="filter.feeder"])
filt_fam<-as.factor(unique(filt_family[,1]))
filt_bio<-select(fulldata, Anopheles_bio, Wyeomyia_bio, Culex_bio)
#Canthocamptidae, Daphnidae/Daphniidae, Cyclopidae removed as crustaceans, Dixidae as not in drought expt

shred_family<-as.data.frame(namedata$family[namedata$functional_group=="shredder"])
shred_fam<-as.factor(unique(shred_family[,1]))
shred_bio<-select(families, TipulidaeLimoniidae_bio, Ptilodactylidae_bio,Calamoceratidae_bio)

scrap_family<-as.data.frame(namedata$family[namedata$functional_group=="scraper"])
scrap_fam<-as.factor(unique(scrap_family[,1]))
scrap_bio<-select(families,  Scirtidae_bio,Candonidae_bio,Limnocytheridae_bio)
#removed Thaumaleidae (colombia), Cyprididae (cardoso) as not in drought expt
#Elmidae_bio is five or fewer bromeliads

##analysis==

eng_bio$eng_sum<-rowSums(eng_bio)
brom$eng_sum<-rowSums(eng_bio)
eng_bio_nozero<-filter(eng_bio,eng_sum>0)%>%select(-eng_sum)
eng.dist<-vegdist(eng_bio_nozero)
adonis(eng.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", filter(brom,eng_sum>0), perm = 1000)
#maxvol, site, log k*site, log mu *site, marginally log mu (shows contingency as in regressions)
mean(eng.dist)#0.803
engulfer.beta<-0.803



gath_bio$gath_sum<-rowSums(gath_bio)
brom$gath_sum<-rowSums(gath_bio)
gath_bio_nozero<-filter(gath_bio,gath_sum>0)%>%select(-gath_sum)
gath.dist<-vegdist(gath_bio_nozero)
adonis(gath.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", filter(brom,gath_sum>0), perm = 1000)
#maxvol, site, logmu x site  (now agrees with regr analysis, where both gatherer and chironmid biomass conting and responsive to rain)
mean(gath.dist)#0.793



pier_bio$pier_sum<-rowSums(pier_bio)
brom$pier_sum<-rowSums(pier_bio)
pier_bio_nozero<-filter(pier_bio,pier_sum>0)%>%select(-pier_sum)
pier.dist<-vegdist(pier_bio_nozero)
adonis(pier.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", filter(brom,pier_sum>0), perm = 1000)
#maxvol, site only (agrees with regr analysis, piercer biomass not contingent or responsive to rain)
mean(pier.dist)#0.842


filt_bio$filt_sum<-rowSums(filt_bio)
brom$filt_sum<-rowSums(filt_bio)
filt_bio_nozero<-filter(filt_bio,filt_sum>0)%>%select(-filt_sum)
filt.dist<-vegdist(filt_bio_nozero)
adonis(filt.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", filter(brom,filt_sum>0), perm = 1000)
#maxvol, site, log mu  (marg log mu x site, p =0.052)
mean(filt.dist)#0.666

scrap_bio$scrap_sum<-rowSums(scrap_bio)
brom$scrap_sum<-rowSums(scrap_bio)
scrap_bio_nozero<-filter(scrap_bio,scrap_sum>0)%>%select(-scrap_sum)
scrap.dist<-vegdist(scrap_bio_nozero)
adonis(scrap.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", filter(brom,scrap_sum>0), perm = 1000)
#maxvol, site, log k x site (agrees with regr analysis, both scraper and scirtid biomass contingent)
mean(scrap.dist)#0.650

shred_bio$shred_sum<-rowSums(shred_bio)
brom$shred_sum<-rowSums(shred_bio)
shred_bio_nozero<-filter(shred_bio,shred_sum>0)%>%select(-shred_sum)
shred.dist<-vegdist(shred_bio_nozero)
adonis(shred.dist~log(mu.scalar)+log(k.scalar)+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar), by="margin", filter(brom,shred_sum>0), perm = 1000)
#max vol, site only (agrees with regr analysis, neither shredder or tipulid biomass contingent nor that responsive to rain)
mean(shred.dist)#0.557

beta<-c(0.803,0.84, 0.79,0.67, 0.65, 0.56  )
fngrp<-c("engulfer","piercer","gatherer","filter.feeder", "scraper", "shredder" )
family.beta<-as.data.frame(beta)
family.beta$fngrp<-fngrp
family.beta$trophic<-c(1,1,0,0,0,0)

dodge <- position_dodge(width = 0.9)
family.beta$fngrp <- factor(family.beta$fngrp, levels = c(response<-as.vector(c("engulfer","piercer","gatherer","filter.feeder", "scraper", "shredder"))))


ggplot(data = family.beta, aes(x = fngrp, y = beta))+
  geom_bar(data=family.beta, aes(x = fngrp, y = beta, fill=trophic), stat="identity", width = 0.75)+
  labs( y = "Family beta diversity between sites", x = "Functional feeding group")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
