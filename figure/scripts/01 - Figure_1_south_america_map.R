# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(fwdata)
library(maps)
library(grid)
library(gridExtra)
library(scales)
library(ggsn)
library(maptools)
library(multipanelfigure)


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

ar<-ggplot(longar, aes(mm))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size=1),
        axis.text.x= element_text(face="bold", colour="black", size=12),
        axis.text.y= element_text(face="bold", colour="black", size=12),
        axis.title.x = element_text(face="bold", colour="black", size=12,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", colour="black", size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  geom_histogram(breaks=seq(0,120, by=10),fill="black")+
  scale_x_continuous(expand=c(0,0),breaks=c(0,60,120), limits=c(0, 125))+
  scale_y_continuous(expand=c(0,0), breaks=c(0,30,60), limits=c(0, 62))+
  labs(x="Daily rain (mm)",y="Frequency", size=3, colour="black")
ar
ar.grob <- ggplotGrob(ar)

crhist<-hist(longcr$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
             include.lowest = TRUE, right = TRUE,
             col = "black", border = NULL,
             main=NULL, xlim=c(0,120),ylim=c(0,60), xlab=NULL, ylab=NULL,
             axes = TRUE, plot = TRUE, labels = FALSE)

cr<-ggplot(longcr, aes(mm))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size=1),
        axis.text.x= element_text(face="bold", colour="black", size=12),
        axis.text.y= element_text(face="bold", colour="black", size=12),
        axis.title.x = element_text(face="bold", colour="black", size=12,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", colour="black", size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  geom_histogram(breaks=seq(0,120, by=10),fill="black")+
  scale_x_continuous(expand=c(0,0),breaks=c(0,60,120), limits=c(0, 125))+
  scale_y_continuous(expand=c(0,0), breaks=c(0,30,60), limits=c(0, 62))+
  labs(x="Daily rain (mm)",y="Frequency", size=3, colour="black")

cr.grob <- ggplotGrob(cr)

fghist<-hist(longfg$mm, breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120),
             include.lowest = TRUE, right = TRUE,
             col = "black", border = NULL,
             main=NULL, xlim=c(0,120),ylim=c(0,60), xlab=NULL, ylab=NULL,
             axes = TRUE, plot = TRUE, labels = FALSE)

fg<-ggplot(longfg, aes(mm))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size=1),
        axis.text.x= element_text(face="bold", colour="black", size=10),
        axis.text.y= element_text(face="bold", colour="black", size=10),
        axis.title.x = element_text(face="bold", colour="black", size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", colour="black", size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  geom_histogram(breaks=seq(0,120, by=10),fill="black")+
  scale_x_continuous(expand=c(0,0),breaks=c(0,60,120), limits=c(0, 125))+
  scale_y_continuous(expand=c(0,0), breaks=c(0,30,60), limits=c(0, 62))+
  labs(x="Daily rain (mm)",y="Frequency", size=3, colour="black")
fg





# load coordinates --------------------------------------------------------------------------------------------------------------------

site_data <- fw_data("0.7.7")$datasets %>% 
  filter(name %in% c("Pitilla2010", "RioBlanco2014", "Cardoso2011", "LasGamas2013", 
                     "Macae2008", "PetitSaut2014", "ElVerde2010"))
site_data

# colors ------------------------------------------------------------------------------------------------------------------------------

# c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")
# ARGENTINA, CARDOSO, COLOMBIA, COSTA RICA, FRENCH GUIANA, MACAE, PUERTO RICO


# draw map ----------------------------------------------------------------------------------------------------------------------------
df <- data.frame(x1 = -105, x2 = -85-2, y1 = 0, y2 = 11-1,x3 = -105, x4 = -60-2, y3 = -40, y4 = -29-1)
label.df<-data.frame(lat=c(10.9830+1, 5.0170, 18.3210, -22.2075+1, -25.0720-1,-29.4330+2, 4.0830), 
                     lng=c(-85.4330-11, -73.7000-11, -65.8170+12, -41.4847+10, -47.9230+12, -60.4670-10, -52.6830+13),
                     country=c("Costa Rica", "Colombia", "Puerto Rico", "Macae BR", "Cardoso BR", "Argentina", "French Guiana"))


our_map <- ggplot() +
  borders(database = "world", fill = "white", colour = "black", size = 0.5) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(), #element_rect(fill = "grey97"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  coord_quickmap(xlim = c(-110, -25), ylim = c(-55, 30)) +
  geom_point(data = site_data, mapping = aes(x = lng, y = lat, fill = location), 
             shape = 21, size = 4, stroke = 0.5) +
  scale_fill_manual(values = c("darkorange2", "firebrick1", "mediumpurple4", "lawngreen", "mediumblue", "deepskyblue1", "forestgreen"))+
  geom_text(data = label.df, mapping = aes(x = lng, y = lat, label=country, fontface = "bold", size =4)) +
  geom_segment(aes(x=x2, xend=x1, y=y2, yend=y1), data=df,size=1,
                 arrow = arrow(length = unit(0.03, "npc")))+
  geom_segment(aes(x=x4, xend=x3, y=y4, yend=y3), data=df,size=1,
               arrow = arrow(length = unit(0.03, "npc")))+
  annotation_custom(
    grob = cr.grob, xmin = -125, xmax = -95, ymin = -15, ymax = 15)+
  annotation_custom(
    grob = ar.grob, xmin = -125, xmax = -95, ymin = -55, ymax = -25)

our_map



# saving figure -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_1_map.pdf",
       plot = our_map, width = 644/3.5, height = 537/3.5, dpi = 300, units = "mm")

ggsave(filename = "figure/figures/Figure_1_map.png",
       plot = our_map, width = 644/3.5, height = 537/3.5, dpi = 600, units = "mm")

# looks like all three panels could be combined with multi_panel_figure command from multipanelfigure package
 
