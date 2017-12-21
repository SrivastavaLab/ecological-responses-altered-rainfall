
# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)

# loading data ------------------------------------------------------------------------------------------------------------------------

figure06 <- readRDS("figure/data/figure06_cca_list.rds")
figure06

# which data in which? ----------------------------------------------------------------------------------------------------------------

## cca arrows
figure06$fn.cca2g$CCA$biplot

## species
figure06$fn.cca2g$CCA$v

## individual bromeliads
figure06$fn.cca2g$CCA$wa

## eigenvalues for each axis
figure06$fn.cca2g$CCA$eig

# tidying data before plotting --------------------------------------------------------------------------------------------------------

## cca arrows
cca_arrows <- data.frame(figure06$fn.cca2g$CCA$biplot) %>% 
  rownames_to_column(var = "ID") %>% 
  select(ID, CCA1, CCA2) %>% 
  mutate(ID = sub(pattern = "site", replacement = "", x = ID),
         ID = sub(pattern = "change\\.", replacement = "", x = ID))
cca_arrows <- cca_arrows[-c(1:6),]
cca_arrows

## species
cca_species <- data.frame(figure06$fn.cca2g$CCA$v) %>% 
  rownames_to_column(var = "ID") %>% 
  select(ID, CCA1, CCA2) %>% 
  mutate(ID = sub(pattern = "_bio", replacement = "", x = ID),
         ID = sub(pattern = "\\.", replacement = " ", x = ID))
cca_species

## bromeliads
bromeliads <- data.frame(figure06$fn.cca2g$CCA$wa) %>% 
  select(CCA1, CCA2) %>% 
  mutate(site = rep(c("Argentina", "Cardoso", "Colombia", "French Guiana",
                      "Macae", "Puerto Rico", "Costa Rica"), each = 30))
bromeliads

## site centroids
site_centroid <- data.frame(figure06$fn.cca2g$CCA$centroids) %>% 
  rownames_to_column(var = "ID") %>% 
  select(ID, CCA1, CCA2) %>% 
  mutate(ID = sub(pattern = "site", replacement = "", x = ID),
         CCA2 = ifelse(ID == "argentina", CCA2 + 0.3, CCA2),
         ID = sub(pattern = "argentina", replacement = "Argentina", x = ID),
         ID = sub(pattern = "cardoso", replacement = "Cardoso", x = ID),
         ID = sub(pattern = "colombia", replacement = "Colombia", x = ID),
         ID = sub(pattern = "costarica", replacement = "Costa Rica", x = ID),
         ID = sub(pattern = "frenchguiana", replacement = "French Guiana", x = ID),
         ID = sub(pattern = "macae", replacement = "Macae", x = ID),
         ID = sub(pattern = "puertorico", replacement = "Puerto Rico", x = ID))
site_centroid

# colors ------------------------------------------------------------------------------------------------------------------------------

# c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")
# ARGENTINA, CARDOSO, COLOMBIA, COSTA RICA, FRENCH GUIANA, MACAE, PUERTO RICO

# panel A: bromeliads and sites -------------------------------------------------------------------------------------------------------

panelA <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70", size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70", size = 0.5, alpha = 0.8) +
  geom_point(data = bromeliads, mapping = aes(x = CCA1, y = CCA2, fill = site), colour = "grey80",
             stroke = 0.5, size = 3, shape = 21, alpha = 0.6) +
  geom_text(data = site_centroid, mapping = aes(x = CCA1, y = CCA2, label = ID),
            size = 4, fontface = 4) +
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "forestgreen",
                               "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = paste0("CCA 1 (", round(figure06$fn.cca2g$CCA$eig[1] * 100, digits = 1),"%)"),
                     limits = c(-2.5, 2.5), expand = c(0,0)) +
  scale_y_continuous(name = paste0("CCA 2 (", round(figure06$fn.cca2g$CCA$eig[2] * 100, digits = 1),"%)"),
                     limits = c(-3, 4.5), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "none",
        axis.title = element_text(colour = "black", face = 2, size = 9)) +
  annotate(geom = "text", x = 2.2, y = 4.2, label = "a", size = 3, fontface = 2)
panelA

# panel B: arrows and species ---------------------------------------------------------------------------------------------------------

## edit label position

cca_arrows <- cca_arrows %>% 
  mutate(CCA1_text = CCA1,
         CCA2_text = CCA2,
         ID = ifelse(ID == "k", yes = "argentina:k", no = ID),
         CCA1_text = ifelse(ID == "argentina:k", CCA1_text, CCA1_text),
         CCA2_text = ifelse(ID == "argentina:k", CCA2_text, CCA2_text),
         CCA1_text = ifelse(ID == "puertorico:k", CCA1_text + 0.1, CCA1_text),
         CCA2_text = ifelse(ID == "puertorico:k", CCA2_text + 0.05, CCA2_text),
         CCA2_text = ifelse(ID == "argentina:k", CCA2_text + 0.05, CCA2_text),
         CCA1_text = ifelse(ID == "macae:k", CCA1_text + 0.05, CCA1_text),
         CCA2_text = ifelse(ID == "macae:k", CCA2_text + 0.05, CCA2_text),
         CCA2_text = ifelse(ID == "costarica:k", CCA2_text + 0.05, CCA2_text),
         CCA1_text = ifelse(ID == "costarica:k", CCA1_text + 0.1, CCA1_text),
         CCA1_text = ifelse(ID == "frenchguiana:k", CCA1_text, CCA1_text),
         CCA2_text = ifelse(ID == "frenchguiana:k", CCA2_text - 0.06, CCA2_text),
         CCA2_text = ifelse(ID == "cardoso:k", CCA2_text - 0.05, CCA2_text),
         CCA1_text = ifelse(ID == "cardoso:k", CCA1_text + 0.1, CCA1_text),
         CCA2_text = ifelse(ID == "colombia:k", CCA2_text - 0.05, CCA2_text),
         CCA1_text = ifelse(ID == "colombia:k", CCA1_text + 0.02, CCA1_text),
         ID = sub(pattern = "puertorico", replacement = "Puerto Rico", x = ID),
         ID = sub(pattern = "cardoso", replacement = "Cardoso", x = ID),
         ID = sub(pattern = "macae", replacement = "Macae", x = ID),
         ID = sub(pattern = "frenchguiana", replacement = "French Guiana", x = ID),
         ID = sub(pattern = "costarica", replacement = "Costa Rica", x = ID),
         ID = sub(pattern = "argentina", replacement = "Argentina", x = ID),
         ID = sub(pattern = "colombia", replacement = "Colombia", x = ID))
cca_arrows

## editing species
cca_species <- cca_species %>% 
  mutate(CCA1_text = CCA1/figure06$fn.cca2g$CCA$rank,
         CCA2_text = CCA2/figure06$fn.cca2g$CCA$rank,
         CCA1_text = ifelse(ID == "gatherer", CCA1_text + 0.15, CCA1_text),
         CCA2_text = ifelse(ID == "piercer", CCA2_text - 0.05, CCA2_text),
         CCA1_text = ifelse(ID == "piercer", CCA1_text + 0.05, CCA1_text),
         CCA1_text = ifelse(ID == "engulfer", CCA1_text + 0.23, CCA1_text),
         CCA2_text = ifelse(ID == "engulfer", CCA2_text + 0.03, CCA2_text),
         CCA1_text = ifelse(ID == "filter feeder", CCA1_text - 0.2, CCA1_text),
         CCA2_text = ifelse(ID == "scraper", CCA2_text - 0.05, CCA2_text),
         CCA1_text = ifelse(ID == "shredder", CCA1_text + 0.04, CCA1_text),
         CCA2_text = ifelse(ID == "shredder", CCA2_text + 0.06, CCA2_text),
         ID = sub(pattern = "engulfer", replacement = "Engulfer", x = ID),
         ID = sub(pattern = "filter feeder", replacement = "Filter Feeder", x = ID),
         ID = sub(pattern = "gatherer", replacement = "Gatherer", x = ID),
         ID = sub(pattern = "piercer", replacement = "Piercer", x = ID),
         ID = sub(pattern = "scraper", replacement = "Scraper", x = ID),
         ID = sub(pattern = "shredder", replacement = "Shredder", x = ID))
cca_species

## to draw a cicle
circleFun <- function(center = c(0,0.1), diameter = 1.6, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

## circle to draw
circle <- circleFun()

panelB <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70", size = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70", size = 0.5, alpha = 0.8) +
  geom_path(data = circle, mapping = aes(x =  x, y = y), size = 0.5, colour = "grey70", alpha = 0.8) +
  annotate(geom = "segment", x = 0.05, xend = 0.17, y = -0.188, yend = -0.15924910, size = 0.4, colour = "grey70") +
  geom_segment(data = cca_arrows, mapping = aes(x = 0, xend = CCA1, y = 0, yend = CCA2),
               colour = "black", size = 0.5, arrow = arrow(length = unit(0.06, "inches"))) +
  geom_point(data = cca_species, mapping = aes(x = CCA1/figure06$fn.cca2g$CCA$rank,
                                               y = CCA2/figure06$fn.cca2g$CCA$rank), 
             size = 2, shape = 21, colour = "black", fill = "grey50") +
  geom_text(data = cca_species, mapping = aes(x = CCA1_text,
                                              y = CCA2_text, label = ID), size = 2.5, fontface = 2) +
  scale_x_continuous(name = paste0("CCA 1 (", round(figure06$fn.cca2g$CCA$eig[1] * 100, digits = 1),"%)"),
                      limits = c(-0.9, 0.9), expand = c(0,0)) +
  scale_y_continuous(name = paste0("CCA 2 (", round(figure06$fn.cca2g$CCA$eig[2] * 100, digits = 1),"%)"),
                     limits = c(-0.7, 1.05), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(colour = "black", face = 2, size = 9)) +
  annotate(geom = "text", x = 0.8, y = 0.97, label = "b", size = 3, fontface = 2) +
  annotate(geom = "text", x = 0.33, y = 0.03, label = paste0("Argentina:","Delta~k"), parse = TRUE, size = 2.5) +
  annotate(geom = "segment", x = 0.08, xend = 0.15, y = -0.02, yend = 0.02, size = 0.4, colour = "grey70") +
  annotate(geom = "text", x = 0.18, y = -0.49, label = paste0("Cardoso:","Delta~k"), parse = TRUE, size = 2.5) +
  annotate(geom = "text", x = 0.4, y = -0.34, label = paste0("Colombia:","Delta~k"), parse = TRUE, size = 2.5) +
  annotate(geom = "text", x = 0.54, y = 0.21139791, label = paste0("Costa~Rica:","Delta~k"), parse = TRUE, size = 2.5) +
  annotate(geom = "text", x = -0.4036756, y = -0.33, label = paste0("French~Guiana:","Delta~k"), parse = TRUE, size = 2.5) +
  annotate(geom = "text", x = 0.1119099, y = 0.32180750, label = paste0("Macae:","Delta~k"), parse = TRUE, size = 2.5) +
  annotate(geom = "text", x = -0.3494158, y = 0.50764944, label = paste0("Puerto~Rico:","Delta~k"), parse = TRUE, size = 2.5)
panelB

# saving panels -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_6.pdf", 
       plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB), ncol = 2)),
       width = 183, height = 90, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Figure_6.png", 
       plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB), ncol = 2)),
       width = 183, height = 90, units = "mm", dpi = 600)
