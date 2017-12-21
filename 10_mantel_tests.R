source("01_datasets_for_paper1.R")
source("02_functions.R")
source("08_multivariate_analyses.R")
sensitive<-read.csv ("Data/family_sensitivity.csv")
bioclim<-read.csv ("Data/bioclim.csv")

# ===mantel tests====

#species pool sizes
#CR 33; PR13; AR 13;FG 25; CO 21; MA 69; CA 82




simmat.fncomp<-vegdist(sim, method="bray")
pairsfc<-function(a,b,d){
  fngrp1<-fngrp
  fngrp1$site<-brom$site
  pairfn<-filter(fngrp1, site%in%c(a,b))%>%select(-site)
  pairx<-filter(brom, site%in%c(a,b))
  pair.dist<-vegdist(pairfn, method="jaccard")
  adon.comp<-adonis(pair.dist~log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", pairx, perm = 1000)
  conting3<-sum(adon.comp$aov.tab$SumsOfSqs[c(7:8, 10:11, 13:14)])/adon.comp$aov.tab$SumsOfSqs[16]
  simmat.fncomp[d]<-conting3
  return(simmat.fncomp)
}

simmat.famcomp<-vegdist(sim, method="bray")
pairsfam<-function(a,b,d){
  families1<-families
  families1$site<-brom$site
  pairfn<-filter(families1, site%in%c(a,b))%>%select(-site)
  pairx<-filter(brom, site%in%c(a,b))
  pair.dist<-vegdist(pairfn, method="jaccard")
  adon.comp<-adonis(pair.dist~log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", pairx, perm = 1000)
  conting3<-sum(adon.comp$aov.tab$SumsOfSqs[c(7:8, 10:11, 13:14)])/adon.comp$aov.tab$SumsOfSqs[16]
  simmat.famcomp[d]<-conting3
  return(simmat.famcomp)
}





pairs1<-function(b,c,d){
  a<-rbind(b,c)
  n<-1
  sum.outq<-data.frame(full.d2=numeric(n),noint.d2=numeric(n), intend.d2=numeric(n), site.d2=numeric(n), conting.p=numeric(n),general.p=numeric(n), siterain.p=numeric(n), rainfall.p=numeric(n), sensitive.p=numeric(n))
  sum.outq[1,]<-allmodel(2,a$gatherer_bio,poisson, a)
  conting<-(sum.outq[,1]-sum.outq[,2])
  simmatg[d]<-conting
  return(simmatg)
}#gatherer pairwise site effect

pairs2<-function(b,c,d){
  a<-rbind(b,c)
  n<-1
  sum.outq<-data.frame(full.d2=numeric(n),noint.d2=numeric(n), intend.d2=numeric(n), site.d2=numeric(n), conting.p=numeric(n),general.p=numeric(n), siterain.p=numeric(n), rainfall.p=numeric(n),sensitive.p=numeric(n))
  sum.outq[1,]<-allmodel(2,a$engulfer_bio,poisson, a)
  print(sum.outq)
  conting<-(sum.outq[,1]-sum.outq[,2])
  simmate[d]<-conting
  return(simmate)
}#engulfer pairwise site effect

pairs3<-function(b,c,d){
  a<-rbind(b,c)
  n<-1
  sum.outf<-data.frame(full.d2=numeric(n),noint.d2=numeric(n), intend.d2=numeric(n), site.d2=numeric(n), conting.p=numeric(n),general.p=numeric(n), siterain.p=numeric(n), rainfall.p=numeric(n),sensitive.p=numeric(n))
  sum.outf[1,]<-allmodel(2,a$filter.feeder_bio,poisson, a)
  print(sum.outf)
  conting<-(sum.outf[,1]-sum.outf[,2])
  simmatf[d]<-conting
  return(simmatf)
}#filterer pairwise site effect



#          ardata     cadata     codata     crdata     fgdata     madata
#cadata      1
#codata      2           7
#crdata      3           8          12
#fgdata      4           9          13        16
#madata      5          10          14        17        19
#prdata      6          11          15        18        20      21

i<-0
simnat.maker<-function(matrix.name, pairs.function)
{
  matrix.name[1]<-pairs.function("cardoso", "argentina") 
  matrix.name[2]<-pairs.function("colombia", "argentina")
  matrix.name[3]<-pairs.function("costarica", "argentina")
  matrix.name[4]<-pairs.function("frenchguiana", "argentina")
  matrix.name[5]<-pairs.function("macae", "argentina")
  matrix.name[6]<-pairs.function("puertorico", "argentina")
  matrix.name[7]<-pairs.function("colombia", "cardoso")
  matrix.name[8]<-pairs.function("costarica", "cardoso")
  matrix.name[9]<-pairs.function("frenchguiana", "cardoso")
  matrix.name[10]<-pairs.function("macae", "cardoso")
  matrix.name[11]<-pairs.function("puertorico", "cardoso")
  matrix.name[12]<-pairs.function("costarica", "colombia")
  matrix.name[13]<-pairs.function("frenchguiana", "colombia")
  matrix.name[14]<-pairs.function("macae", "colombia")
  matrix.name[15]<-pairs.function("puertorico", "colombia")
  matrix.name[16]<-pairs.function("frenchguiana", "costarica")
  matrix.name[17]<-pairs.function("macae", "costarica")
  matrix.name[18]<-pairs.function("puertorico", "costarica")
  matrix.name[19]<-pairs.function("macae", "frenchguiana")
  matrix.name[20]<-pairs.function("puertorico", "frenchguiana")
  matrix.name[21]<-pairs.function("puertorico", "macae")
  return(matrix.name)
}

pairs.sens<-function(b,c){
  a<-abs(sensitive$sens.index [sensitive$site==b]-sensitive$sens.index [sensitive$site==c])
  return(a)
}#sensitiity index pairwise site effect

pairs.sens.hydro<-function(b,c){
  a<-abs(sensitive$sens.index.hydro [sensitive$site==b]-sensitive$sens.index.hydro [sensitive$site==c])
  return(a)
}

simmat.sens<-vegdist(sim, method="bray")
simmat.sensout<-simnat.maker(simmat.sens, pairs.sens)
simmat.sensout.hydro<-simnat.maker(simmat.sens, pairs.sens.hydro)

simmat.fncomp<-pairsfc("cardoso", "argentina",1);simmat.fncomp<-pairsfc("colombia", "argentina",2);simmat.fncomp<-pairsfc("costarica", "argentina",3)
simmat.fncomp<-pairsfc("frenchguiana", "argentina",4);simmat.fncomp<-pairsfc("macae", "argentina",5);simmat.fncomp<-pairsfc("puertorico", "argentina",6)
simmat.fncomp<-pairsfc("colombia", "cardoso",7);simmat.fncomp<-pairsfc("costarica", "cardoso",8);simmat.fncomp<-pairsfc("frenchguiana", "cardoso",9)
simmat.fncomp<-pairsfc("macae", "cardoso",10);simmat.fncomp<-pairsfc("puertorico", "cardoso",11);simmat.fncomp<-pairsfc("costarica", "colombia",12)
simmat.fncomp<-pairsfc("frenchguiana", "colombia",13);simmat.fncomp<-pairsfc("macae", "colombia",14);simmat.fncomp<-pairsfc("puertorico", "colombia",15)
simmat.fncomp<-pairsfc("frenchguiana", "costarica",16);simmat.fncomp<-pairsfc("macae", "costarica",17);simmat.fncomp<-pairsfc("puertorico", "costarica",18)
simmat.fncomp<-pairsfc("macae", "frenchguiana",19);simmat.fncomp<-pairsfc("puertorico", "frenchguiana",20);simmat.fncomp<-pairsfc("puertorico", "macae",21)
simmat.fncomp

simmat.famcomp<-pairsfam("cardoso", "argentina",1);simmat.famcomp<-pairsfam("colombia", "argentina",2);simmat.famcomp<-pairsfam("costarica", "argentina",3)
simmat.famcomp<-pairsfam("frenchguiana", "argentina",4);simmat.famcomp<-pairsfam("macae", "argentina",5);simmat.famcomp<-pairsfam("puertorico", "argentina",6)
simmat.famcomp<-pairsfam("colombia", "cardoso",7);simmat.famcomp<-pairsfam("costarica", "cardoso",8);simmat.famcomp<-pairsfam("frenchguiana", "cardoso",9)
simmat.famcomp<-pairsfam("macae", "cardoso",10);simmat.famcomp<-pairsfam("puertorico", "cardoso",11);simmat.famcomp<-pairsfam("costarica", "colombia",12)
simmat.famcomp<-pairsfam("frenchguiana", "colombia",13);simmat.famcomp<-pairsfam("macae", "colombia",14);simmat.famcomp<-pairsfam("puertorico", "colombia",15)
simmat.famcomp<-pairsfam("frenchguiana", "costarica",16);simmat.famcomp<-pairsfam("macae", "costarica",17);simmat.famcomp<-pairsfam("puertorico", "costarica",18)
simmat.famcomp<-pairsfam("macae", "frenchguiana",19);simmat.famcomp<-pairsfam("puertorico", "frenchguiana",20);simmat.famcomp<-pairsfam("puertorico", "macae",21)
simmat.famcomp

seven.site<-matrix( data=NA, nrow=7, ncol=7, byrow = TRUE) 
colnames(seven.site)<-c("a","b", "c", "d", "e", "f", "g")
first.seven<-as.data.frame(seven.site)
second.seven<-as.data.frame(seven.site)

first.seven[lower.tri(first.seven)]<-simmat.fncomp
second.seven[lower.tri(second.seven)]<-simmat.famcomp

linear.first<-gather(first.seven, siteid, diss)
linear.second<-gather(second.seven, siteid, diss)
plot(linear.first$diss~linear.second$diss)



simmatg<-vegdist(sim, method="bray")
simmatg<-pairs1(cadata, ardata,1);simmatg<-pairs1(codata, ardata,2);simmatg<-pairs1(crdata, ardata,3)
simmatg<-pairs1(fgdata, ardata,4);simmatg<-pairs1(madata, ardata,5);simmatg<-pairs1(prdata, ardata,6)
simmatg<-pairs1(codata, cadata,7);simmatg<-pairs1(crdata, cadata,8);simmatg<-pairs1(fgdata, cadata,9)
simmatg<-pairs1(madata, cadata,10);simmatg<-pairs1(prdata, cadata,11);simmatg<-pairs1(crdata, codata,12)
simmatg<-pairs1(fgdata, codata,13);simmatg<-pairs1(madata, codata,14);simmatg<-pairs1(prdata, codata,15)
simmatg<-pairs1(fgdata, crdata,16);simmatg<-pairs1(madata, crdata,17);simmatg<-pairs1(prdata, crdata,18)
simmatg<-pairs1(madata, fgdata,19);simmatg<-pairs1(prdata, fgdata,20);simmatg<-pairs1(prdata, madata,21)
simmatg

#engulfers: no arg, co, pr

#          ca    cr    fg
#crdata      1
#fgdata      2     4
#madata      3     5     6

sim2<-as.data.frame(rep(1,4))
simmate<-vegdist(sim2, method="bray")
simmate<-pairs2(crdata, cadata,1);simmate<-pairs2(fgdata, cadata,2);simmate<-pairs2(madata, cadata,3)
simmate<-pairs2(fgdata, crdata,4);simmate<-pairs2(madata, crdata,5);simmate<-pairs2(fgdata, madata,6)
simmate

pool<-cbind(select(fulldata, site), families)
sppool<-group_by(pool, site)%>%
  summarise(
    a=sum(Calamoceratidae_bio),
    b=sum(Candonidae_bio),
    c=sum(Cecidomyiidae_bio),
    d=sum(Ceratopogonidae_bio),
    e=sum(Chironomidae_bio),
    f=sum(Coenagrionidae_bio),
    g=sum(Corethrellidae_bio),
    h=sum(Culicidae_bio),
    i=sum(Dytiscidae_bio),
    j=sum(Empididae_bio),
    k=sum(Enchytraeoidae_bio),
    l=sum(Ephydridae_bio),
    m=sum(Hydrophilidae_bio),
    n=sum(Lampyridae_bio),
    o=sum(Limnocytheridae_bio),
    p=sum(Naididae_bio),
    q=sum(Pseudostigmatidae_bio),
    r=sum(Psychodidae_bio),
    s=sum(Ptilodactylidae_bio),
    t=sum(Scirtidae_bio),
    u=sum(Stratiomyidae_bio),
    v=sum(Syrphidae_bio),
    w=sum(Tabanidae_bio),
    x=sum(TipulidaeLimoniidae_bio),
    y=sum(Hirudinea_bio))
sppool2<-select(sppool, -site)
row.names(sppool2)<-sppool$site
poolhell<-decostand(sppool2,method="hellinger")
simpool<-vegdist(poolhell, method="jaccard")

full.families<-select(fulldata, Calamoceratidae_bio,Candonidae_bio,Cecidomyiidae_bio,Ceratopogonidae_bio,Chironomidae_bio,
                      Coenagrionidae_bio,Corethrellidae_bio,Culicidae_bio, Dytiscidae_bio,
                      Empididae_bio,Enchytraeoidae_bio,Ephydridae_bio,Hirudinea_bio,Hydrophilidae_bio,Lampyridae_bio,Limnocytheridae_bio, 
                      Naididae_bio, Pseudostigmatidae_bio,Psychodidae_bio,Ptilodactylidae_bio,
                      Scirtidae_bio,Stratiomyidae_bio,Syrphidae_bio,Tabanidae_bio,TipulidaeLimoniidae_bio, Aeolosomatidae_bio, Anisopodidae_bio,
                      Dolichopodidae_bio,Elateridae_bio,  Elmidae_bio, Muscidae_bio, Periscelididae_bio, Phoridae_bio, Scatopsidae_bio, 
                      Sciaridae_bio, Sphaeroceridae_bio)
full.pool<-cbind(select(fulldata, site), full.families)
full.sppool<-group_by(full.pool, site)%>%
  summarise(
    a=sum(Calamoceratidae_bio),
    b=sum(Candonidae_bio),
    c=sum(Cecidomyiidae_bio),
    d=sum(Ceratopogonidae_bio),
    e=sum(Chironomidae_bio),
    f=sum(Coenagrionidae_bio),
    g=sum(Corethrellidae_bio),
    h=sum(Culicidae_bio),
    i=sum(Dytiscidae_bio),
    j=sum(Empididae_bio),
    k=sum(Enchytraeoidae_bio),
    l=sum(Ephydridae_bio),
    m=sum(Hydrophilidae_bio),
    n=sum(Lampyridae_bio),
    o=sum(Limnocytheridae_bio),
    p=sum(Naididae_bio),
    q=sum(Pseudostigmatidae_bio),
    r=sum(Psychodidae_bio),
    s=sum(Ptilodactylidae_bio),
    t=sum(Scirtidae_bio),
    u=sum(Stratiomyidae_bio),
    v=sum(Syrphidae_bio),
    w=sum(Tabanidae_bio),
    x=sum(TipulidaeLimoniidae_bio),
    y=sum(Hirudinea_bio), 
    z=sum(Aeolosomatidae_bio), 
    aa=sum(Anisopodidae_bio),
    bb=sum(Dolichopodidae_bio),
    cc=sum(Elateridae_bio),  
    dd=sum(Elmidae_bio), 
    ee=sum(Muscidae_bio), 
    ff=sum(Periscelididae_bio), 
    gg=sum(Phoridae_bio), 
    hh=sum(Scatopsidae_bio), 
    ii=sum(Sciaridae_bio), 
    jj=sum(Sphaeroceridae_bio))

full.sppool2<-select(full.sppool, -site)
row.names(full.sppool2)<-full.sppool$site
per.sim<-designdist(full.sppool2, method="J/((A+B)/2)", terms=c("binary"))

bioclim.nosite<-mutate(bioclim,sqrtppt=sqrt(ppt))%>%select(-site,-ppt)
simbioclim<-vegdist(bioclim.nosite, method="euclidean")
ppt.nosite<-select(bioclim,-site,-temp)
simppt<-vegdist(ppt.nosite, method="euclidean")

fncomp.pool<-mantel(simmat.fncomp, simpool, method="pearson", permutations=1000) #not sig r=-0.393, p =0.96
fncomp.fam<-mantel(simmat.fncomp, simmat.famcomp, method="pearson", permutations=1000) #r = 0.78, p = 0.002
#fncomp.sens<-mantel(simmat.fncomp, simmat.sensout, method="pearson", permutations=1000) # r = -0.16, p =0.76
fncomp.bioclim<-mantel(simmat.fncomp, simbioclim, method="pearson", permutations=1000) #r = 0.042 p = 0.417
fncomp.ppt<-mantel(simmat.fncomp, simppt, method="pearson", permutations=1000) #r = -0.0448, p = 0.577



fncomp.pool$statistic
fncomp.pool$signif

sppool.gatherer<-select(sppool2,d,e,l,u,r,v,p,k)
row.names(sppool.gatherer)<-sppool$site
poolhell.gath<-decostand(sppool.gatherer,method="hellinger")
simpool.gath<-vegdist(poolhell.gath, method="bray")
simpool.gath

mantel(simmatg, simpool.gath, method="pearson", permutations=5039) #not sig just for gatherers r = -0.41, p=0.92

simmatg3<-vegdist(sim, method="bray")
pairsg3<-function(a,b,d){
  gath_bio$site<-brom$site
  gath_bio_nozero<-filter(gath_bio,gath_sum>0)%>%filter(site%in%c(a,b))%>%select(-gath_sum, -site)
  xvar.gath<-filter(brom,gath_sum>0)%>%filter(site%in%c(a,b))
  gath.dist<-vegdist(gath_bio_nozero)
  adon.gath<-adonis(gath.dist~log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", xvar.gath, perm = 1000)
  conting3<-sum(adon.gath$aov.tab$SumsOfSqs[c(7:8, 10:11, 13:14)])/adon.gath$aov.tab$SumsOfSqs[16]
  simmatg3[d]<-conting3
  return(simmatg3)
}

simmatg3<-pairsg3("cardoso", "argentina",1);simmatg3<-pairsg3("colombia", "argentina",2);simmatg3<-pairsg3("costarica", "argentina",3)
simmatg3<-pairsg3("frenchguiana", "argentina",4);simmatg3<-pairsg3("macae", "argentina",5);simmatg3<-pairsg3("puertorico", "argentina",6)
simmatg3<-pairsg3("colombia", "cardoso",7);simmatg3<-pairsg3("costarica", "cardoso",8);simmatg3<-pairsg3("frenchguiana", "cardoso",9)
simmatg3<-pairsg3("macae", "cardoso",10);simmatg3<-pairsg3("puertorico", "cardoso",11);simmatg3<-pairsg3("costarica", "colombia",12)
simmatg3<-pairsg3("frenchguiana", "colombia",13);simmatg3<-pairsg3("macae", "colombia",14);simmatg3<-pairsg3("puertorico", "colombia",15)
simmatg3<-pairsg3("frenchguiana", "costarica",16);simmatg3<-pairsg3("macae", "costarica",17);simmatg3<-pairsg3("puertorico", "costarica",18)
simmatg3<-pairsg3("macae", "frenchguiana",19);simmatg3<-pairsg3("puertorico", "frenchguiana",20);simmatg3<-pairsg3("puertorico", "macae",21)
simmatg3

mantel(simmatg, simmatg3, method="pearson", permutations=999) #yeah! r=0.515, p = 0.037
mantel(simpool.gath, simmatg3, method="pearson", permutations=999) #taxonomic conting unrelated to similarity in pool
mantel.partial(simmatg, simmatg3, simpool.gath, method="pearson", permutations=999)#still sig after pulling out taxonomic diff in species pools
mantel.partial(simpool.gath,  simmatg3, simmatg, method="pearson", permutations=999)

exclude<-c("argentina", "colombia", "puertorico")
sppool.engulfer<-select(sppool,site,f,g,i,m,q,y,z)%>%filter(site%nin%exclude)%>%select(-site)
poolhell.eng<-decostand(sppool.engulfer,method="hellinger")
simpool.eng<-vegdist(poolhell.eng, method="bray")
simpool.eng
mantel(simmate, simpool.eng, method="pearson", permutations=999) #r = 0.442, p=0.33

sim2<-as.data.frame(rep(1,4))
simmate3<-vegdist(sim2, method="bray")

pairse3<-function(a,b,d){
  eng_bio$site<-brom$site
  eng_bio_nozero<-filter(eng_bio,eng_sum>0)%>%filter(site%in%c(a,b))%>%select(-eng_sum, -site)
  xvar.eng<-filter(brom,eng_sum>0)%>%filter(site%in%c(a,b))
  eng.dist<-vegdist(eng_bio_nozero)
  adon.eng<-adonis(eng.dist~log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", xvar.eng, perm = 1000)
  conting3<-sum(adon.eng$aov.tab$SumsOfSqs[c(7:8, 10:11, 13:14)])/adon.eng$aov.tab$SumsOfSqs[16]
  simmate3[d]<-conting3
  return(simmate3)
}

simmate3<-pairse3("costarica", "cardoso",1);simmate3<-pairse3("frenchguiana", "cardoso",2);simmate3<-pairse3("macae", "cardoso",3)
simmate3<-pairse3("frenchguiana", "costarica",4);simmate3<-pairse3("macae", "costarica",5);simmate3<-pairse3("frenchguiana", "macae",6)
simmate3

mantel(simmate, simmate3, method="pearson", permutations=23) #r=0.858, p =0.12 looks like an effect, but not enouh power


simdist<-vegdist(sim, method="euclidean")
simdist[1]<-1336;simdist[2]<-4078; simdist[3]<-5234; simdist[4]<-3920; simdist[5]<-2027; simdist[6]<-5342
simdist[7]<-4365;simdist[8]<-5711; simdist[9]<-3395; simdist[10]<-691; simdist[11]<-5195; simdist[12]<-1420
simdist[13]<-2340;simdist[14]<-4641; simdist[15]<-1759; simdist[16]<-3624; simdist[17]<-6035; simdist[18]<-2263
simdist[19]<-3290;simdist[20]<-2020; simdist[21]<-5224
simdist

fn.comp.dist<-mantel(simmat.fncomp, simdist, method="pearson", permutations=1000) # r = 0.08, p = 0.32
mantel(simmat, simdist, method="pearson", permutations=1000) #not sig
mantel(simmat, simpool, method="pearson", permutations=1000) #not sig
mantel.partial(simmat, simpool, simdist, method="pearson", permutations=1000) #not sig

#hydrologic sensitivity


pairsh<-function(a,b,d){
  hydro1<-purehydro
  hydro1$site<-brom.noca$site
  pairhydro<-filter(hydro1, site%in%c(a,b))%>%select(-site)
  pairx<-filter(brom.noca, site%in%c(a,b))
  pair.h.dist<-vegdist(pairhydro, method="jaccard")
  adon.comp<-adonis(pair.h.dist~log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", pairx, perm = 1000)
  conting3<-sum(adon.comp$aov.tab$SumsOfSqs[c(7:8, 10:11, 13:14)])/adon.comp$aov.tab$SumsOfSqs[16]
  simmath[d]<-conting3
  return(simmath)
}

sim3<-subset(sim,rownames(sim)!="cadata")

simmat.h.fncomp<-vegdist(sim3, method="bray")
pairshfc<-function(a,b,d){
  fngrp1<-fngrp
  fngrp1$site<-brom$site
  pairfn<-filter(fngrp1, site%in%c(a,b))%>%select(-site)
  pairx<-filter(brom, site%in%c(a,b))
  pair.dist<-vegdist(pairfn, method="jaccard")
  adon.comp<-adonis(pair.dist~log(mu.scalar)+log(k.scalar)+change.k+change.mu+maxvol+site+site:log(k.scalar)+site:log(mu.scalar)+log(mu.scalar):log(k.scalar)+site:log(mu.scalar):log(k.scalar)+site:change.k+site:change.mu+change.mu:change.k+site:change.mu:change.k, by="margin", pairx, perm = 1000)
  conting3<-sum(adon.comp$aov.tab$SumsOfSqs[c(7:8, 10:11, 13:14)])/adon.comp$aov.tab$SumsOfSqs[16]
  simmat.h.fncomp[d]<-conting3
  return(simmat.h.fncomp)
}

#           ardata codata crdata fgdata madata
#codata      1
#crdata      2      6
#fgdata      3      7      10
#madata      4      8      11      13
#prdata      5      9      12      14      15

sim3<-subset(sim,rownames(sim)!="cadata")
simmath<-vegdist(sim3, method="bray")
simmath<-pairsh("colombia", "argentina",1);simmath<-pairsh("costarica", "argentina",2);simmath<-pairsh("frenchguiana", "argentina",3)
simmath<-pairsh("macae", "argentina",4);simmath<-pairsh("puertorico", "argentina",5);simmath<-pairsh("colombia", "costarica",6)
simmath<-pairsh("colombia", "frenchguiana",7);simmath<-pairsh("macae", "colombia",8);simmath<-pairsh("puertorico", "colombia",9)
simmath<-pairsh("frenchguiana", "costarica",10);simmath<-pairsh("macae", "costarica",11);simmath<-pairsh("costarica", "puertorico",12)
simmath<-pairsh("frenchguiana", "macae",13);simmath<-pairsh("frenchguiana", "puertorico",14);simmath<-pairsh("puertorico","macae",15)
simmath

simmat.h.fncomp<-vegdist(sim3, method="bray")
simmat.h.fncomp<-pairshfc("colombia", "argentina",1);simmat.h.fncomp<-pairshfc("costarica", "argentina",2);simmat.h.fncomp<-pairshfc("frenchguiana", "argentina",3)
simmat.h.fncomp<-pairshfc("macae", "argentina",4);simmat.h.fncomp<-pairshfc("puertorico", "argentina",5);simmat.h.fncomp<-pairshfc("colombia", "costarica",6)
simmat.h.fncomp<-pairshfc("colombia", "frenchguiana",7);simmat.h.fncomp<-pairshfc("macae", "colombia",8);simmat.h.fncomp<-pairshfc("puertorico", "colombia",9)
simmat.h.fncomp<-pairshfc("frenchguiana", "costarica",10);simmat.h.fncomp<-pairshfc("macae", "costarica",11);simmat.h.fncomp<-pairshfc("costarica", "puertorico",12)
simmat.h.fncomp<-pairshfc("frenchguiana", "macae",13);simmat.h.fncomp<-pairshfc("frenchguiana", "puertorico",14);simmat.h.fncomp<-pairshfc("puertorico","macae",15)
simmat.h.fncomp

fncomp.hydro<-mantel(simmath, simmat.h.fncomp, method="pearson", permutations=1000) # r=0.59, p = 0.026 (worth pursuing!)jaccard+jaccard

six.site<-matrix( data=NA, nrow=6, ncol=6, byrow = TRUE) 
colnames(six.site)<-c("argentina","colombia", "costarica", "frenchguiana", "macae", "puertorico")
first.six<-as.data.frame(six.site)
second.six<-as.data.frame(six.site)

first.six[lower.tri(first.six)]<-simmath
second.six[lower.tri(second.six)]<-simmat.h.fncomp

linear.first<-gather(first.six, siteid, diss1)
linear.second<-gather(second.six, siteid, diss2)
linear.mat<-cbind(linear.first, linear.second)
plot(linear.first$diss~linear.second$diss)
m1<-lm(linear.first$diss~linear.second$diss)
ggplot(linear.mat, aes(x=diss1, y=diss2)) + geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  labs(x="Contingency in hydrological sensitivity to rain", y = "Contingency in functional composition sensitivity to rain")+
  theme_bw()




pairsht<-function(b,c,d){
  simmatax[d]<-abs(sensitive$sens.index [sensitive$site==b]-sensitive$sens.index [sensitive$site==c])
  return(simmatax)
}

sim3a<-subset(sim,rownames(sim)!="cadata")
simmatax<-vegdist(sim3a, method="bray")
simmatax<-pairsht("colombia", "argentina",1);simmatax<-pairsht("costarica", "argentina",2);simmatax<-pairsht("frenchguiana", "argentina",3)
simmatax<-pairsht("macae", "argentina",4);simmatax<-pairsht("puertorico", "argentina",5);simmatax<-pairsht("colombia", "costarica",6)
simmatax<-pairsht("colombia", "frenchguiana",7);simmatax<-pairsht("macae", "colombia",8);simmatax<-pairsht("puertorico", "colombia",9)
simmatax<-pairsht("frenchguiana", "costarica",10);simmatax<-pairsht("macae", "costarica",11);simmatax<-pairsht("costarica", "puertorico",12)
simmatax<-pairsht("frenchguiana", "macae",13);simmatax<-pairsht("frenchguiana", "puertorico",14);simmatax<-pairsht("puertorico","macae",15)
simmatax

pairsht2<-function(b,c,d){
  simmatax[d]<-abs(sensitive$sens.index.hydro [sensitive$site==b]-sensitive$sens.index.hydro [sensitive$site==c])
  return(simmatax)
}

sim3a<-subset(sim,rownames(sim)!="cadata")
simmatax<-vegdist(sim3a, method="bray")
simmatax<-pairsht("colombia", "argentina",1);simmatax<-pairsht("costarica", "argentina",2);simmatax<-pairsht("frenchguiana", "argentina",3)
simmatax<-pairsht("macae", "argentina",4);simmatax<-pairsht("puertorico", "argentina",5);simmatax<-pairsht("colombia", "costarica",6)
simmatax<-pairsht("colombia", "frenchguiana",7);simmatax<-pairsht("macae", "colombia",8);simmatax<-pairsht("puertorico", "colombia",9)
simmatax<-pairsht("frenchguiana", "costarica",10);simmatax<-pairsht("macae", "costarica",11);simmatax<-pairsht("costarica", "puertorico",12)
simmatax<-pairsht("frenchguiana", "macae",13);simmatax<-pairsht("frenchguiana", "puertorico",14);simmatax<-pairsht("puertorico","macae",15)
simmatax

mantel(simmatax, simmat.h.fncomp, method="pearson", permutations=1000) # ns, either with rain or hydro How fn grps affected by rain not related to hydro sensitivity of taxa
mantel(simmatax, simmath, method="pearson", permutations=1000) 
mantel.partial(simmath, simmat.h.fncomp, simmatax, method="pearson", permutations=1000) #still sig
mantel.partial(simmat.h.fncomp, simmatax, simmath, method="pearson", permutations=1000) #ns rain / hydro even after partialling out rain->hydro


simmatf<-vegdist(sim3, method="bray")
simmatf<-pairs3(codata, ardata,1);simmatf<-pairs3(crdata, ardata,2);simmatf<-pairs3(fgdata, ardata,3)
simmatf<-pairs3(madata, ardata,4);simmatf<-pairs3(prdata, ardata,5);simmatf<-pairs3(codata, crdata,6)
simmatf<-pairs3(codata, fgdata,7);simmatf<-pairs3(madata, codata,8);simmatf<-pairs3(prdata, codata,9)
simmatf<-pairs3(fgdata, crdata,10);simmatf<-pairs3(madata, crdata,11);simmatf<-pairs3(crdata, prdata,12)
simmatf<-pairs3(fgdata, madata,13);simmatf<-pairs3(fgdata, prdata,14);simmatf<-pairs3(prdata, madata,15)
simmatf

mantel(simmath, simmatf, method="pearson", permutations=5000)#blast r=0.28, p =0.16 even filter feeder contin not related to propdried days conting
#also tried for gatherers, nonsig
