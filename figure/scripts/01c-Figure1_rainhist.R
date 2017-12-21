# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)
library(scales)

# loading data ------------------------------------------------------------------------------------------------------------------------

arrain<-read.csv("../rainfall.schedules/Argentinaschedule.csv")
arrain1<-filter(arrain, trt.name%in%"mu1k1")%>%select(day.3:day.63)

crrain<-read.csv("../rainfall.schedules/CostaRicaschedule.csv")
crrain1<-filter(crrain, trt.name%in%"mu1k1")%>%select(X6:X66)

fgrain<-read.csv("../rainfall.schedules/FrenchGuianaschedule.csv")
fgrain1<-filter(fgrain, trt.name%in%"mu1k1")%>%select(X6:X66)

# tidying data ------------------------------------------------------------------------------------------------------------------------

longar <- arrain1 %>% 
  gather(key = rain, value = mm, day.3:day.63)%>%mutate(mm=as.numeric(mm))

longcr <- crrain1 %>% 
  gather(key = rain, value = mm, X6:X66)%>%mutate(mm=as.numeric(mm))

longfg <- fgrain1 %>% 
  gather(key = rain, value = mm, X6:X66)%>%mutate(mm=as.numeric(mm))

arhist<-hist(longar$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
     include.lowest = TRUE, right = TRUE,
     col = "black", border = NULL,
     main=NULL, xlim=c(0,120), ylim=c(0,60), xlab=NULL, ylab=NULL,
     axes = TRUE, plot = TRUE, labels = FALSE)

crhist<-hist(longcr$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
     include.lowest = TRUE, right = TRUE,
     col = "black", border = NULL,
     main=NULL, xlim=c(0,120),ylim=c(0,60), xlab=NULL, ylab=NULL,
     axes = TRUE, plot = TRUE, labels = FALSE)

fghist<-hist(longfg$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
     include.lowest = TRUE, right = TRUE,
     col = "black", border = NULL,
     main=NULL, xlim=c(0,120),ylim=c(0,60), xlab=NULL, ylab=NULL,
     axes = TRUE, plot = TRUE, labels = FALSE)



# saving figure -----------------------------------------------------------------------------------------------------------------------

png(filename = "figure/figures/Fig1_ar_hist.png",
    width = 100, height = 100, pointsize=12, res=300, units = "mm")
arhist<-hist(longar$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
             include.lowest = TRUE, right = TRUE,
             col = "black", border = NULL,
             main=NULL, xlim=c(0,120), ylim=c(0,60), xlab=NULL, ylab=NULL,
             axes = TRUE, plot = TRUE, labels = FALSE)
arhist
dev.off()


png(filename = "figure/figures/Fig1_cr_hist.png",
    width = 100, height = 100, pointsize=12, res=300, units = "mm")
crhist<-hist(longcr$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
             include.lowest = TRUE, right = TRUE,
             col = "black", border = NULL,
             main=NULL, xlim=c(0,120),ylim=c(0,60), xlab=NULL, ylab=NULL,
             axes = TRUE, plot = TRUE, labels = FALSE)
crhist
dev.off()
