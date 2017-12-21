source("01_datasets_for_paper1.R")
source("02_functions.R")


#===rate of invertebrate flushing====

mean(flu$total) #2.079365
sd(flu$total)/sqrt(length(flu$total))#se = 0.276

min(flu$water)#225.575
max(flu$water) #9142.877

hist(flu$total, breaks=15)

#therefore, on average 2.1 +/- 0.28 invertebrate individuals were flushed from bromeliads per overflow event (225-9142 ml of water)

#===composition of flushed species

names(flu)
flu2<-mutate(flu, Chironomidae=orthoclad+chiron,
             Tipulidae=tipulid,
             Pseudostigmatidae=damsel,
             Psychodidae=psyco,
             Scirtidae=scirtid,
             Ceratopogonidae= cerato,
             Naididae = oligo,
             Culicidae = mossie,
             Syrphidae = syrphid)

propflush<-group_by(flu2,bromeliad.id)%>%summarise(propchiron=sum(Chironomidae)/sum(total),
                                        proptip=sum(Tipulidae)/sum(total),
                                        proppseudo=sum(Pseudostigmatidae)/sum(total),
                                        proppsyc=sum(Psychodidae)/sum(total),
                                        propsci=sum(Scirtidae)/sum(total),
                                        propcerato=sum(Ceratopogonidae)/sum(total),
                                        propoligo= sum(Naididae)/sum(total),
                                        propcul=sum(Culicidae)/sum(total),
                                        propsyr=sum(Syrphidae)/sum(total))


pitfam.sum<-pitfam%>%group_by(site_brom.id, family)%>%summarise(abund=sum(abundance))%>%spread(key=family, value=abund, fill=0)
pitfam.sum$total<-rowSums(pitfam.sum[,2:17])

pitfam.sum<-pitfam.sum%>%mutate(propchiron=sum(Chironomidae)/(total),
                                                         proptip=sum(Tipulidae)/(total),
                                                         proppseudo=sum(Pseudostigmatidae)/(total),
                                                         proppsyc=sum(Psychodidae)/(total),
                                                         propsci=sum(Scirtidae)/(total),
                                                         propcerato=sum(Ceratopogonidae)/(total),
                                                         propoligo= sum(Naididae)/(total),
                                                         propcul=sum(Culicidae)/(total),
                                                         propsyr=sum(Syrphidae)/(total))

m1<-propflush%>%summarise(mean(propchiron),
                      mean(proptip),
                      mean(proppseudo),
                      mean(proppsyc),
                      mean(propsci),
                      mean(propcerato),
                      mean(propoligo),
                      mean(propcul),
                      mean(propsyr))

m2<-pitfam.sum%>%ungroup%>%summarise(mean(propchiron),
                       mean(proptip),
                       mean(proppseudo),
                       mean(proppsyc),
                       mean(propsci),
                       mean(propcerato),
                       mean(propoligo),
                       mean(propcul),
                       mean(propsyr))

#==chiron, oligo abundant in flush as abundant in bromeliads, tipulids and psychodids are susceptible to flushing, other taxa not flushed

allinvert<-rbind(m1,m2)
rownames(allinvert)<-c("Flushed", "Retained")

write.csv(allinvert, "Data/FlushedRetained")
