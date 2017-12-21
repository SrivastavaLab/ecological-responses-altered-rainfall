
# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)

# load data ---------------------------------------------------------------------------------------------------------------------------

concord.out.all <- read_csv( "Data/concord_out_all.csv")
concord.out.all <- filter(concord.out.all, Response != "sqrt.n15.bromeliad.final") #DS added sept 25th to remove redundant N variable

# splitting data frame ---------------------------------------------------------------------------------------------------------------

concord.rain <- filter(concord.out.all, Model == "rain")
concord.hydro <- filter(concord.out.all, Model == "hydro")

#colours-------------------------------------------------------------------------
#using the colourspace package for 9 qualitative colurs, get:
# "#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363"

# general theme -----------------------------------------------------------------------------------------------------------------------

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title = element_text(face = "bold", size = 9.4, colour = "black"),
                       axis.text = element_text(size = 8, colour = "black"), 
                       axis.line.x = element_line(colour = "black", size = 0.3),
                       axis.line.y = element_line(colour = "black", size = 0.3),
                       panel.grid = element_blank(),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 10, 0, 0)))

# PANEL A: rainfall -------------------------------------------------------------------------------------------------------------------

panelA <- ggplot(data = concord.rain, mapping = aes(x = Sites, y = Spearman, colour = Response)) +
  annotate(geom = "rect", xmin = 0.5, xmax = 5.5, ymin = -0.362, ymax = 0.362, alpha = 0.1, colour = "grey90", size = 0) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.3, size = 1, stroke = 0.5, show.legend = FALSE) +
  geom_smooth(span = 1.2, show.legend = FALSE, se = FALSE, size = 0.8) +
  scale_y_continuous(name = "Spearman rank correlation (Rainfall)", 
                     breaks = seq(-0.8, 0.6, by = 0.2),
                     labels = c("-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6"), limits = c(-0.8, 0.7)) +
  scale_x_continuous(name = "", expand = c(0,0), limits = c(0.5, 5.6)) +
  scale_colour_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  #scale_colour_manual(values = c("mediumpurple4", "darkorange2", "mediumblue", "lawngreen", "firebrick1","black", "grey", "black", "grey")) +
  scale_fill_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  #scale_fill_manual(values = c("mediumpurple4", "darkorange2", "mediumblue", "lawngreen", "firebrick1", "black", "grey","black", "grey")) +
  general_theme +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
  annotate(geom = "text", x = 5.4, y = 0.7, label = "a", size = 2.7, fontface = 2)
panelA

# PANEL B: hydrology ------------------------------------------------------------------------------------------------------------------
#note that piercers are not plotted this panel

panelB <- ggplot(data = concord.hydro, mapping = aes(x = Sites, y = Spearman, colour = Response)) +
  annotate(geom = "rect", xmin = 0.5, xmax = 5.5, ymin = -0.362, ymax = 0.362, alpha = 0.1, colour = "grey90", size = 0) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.3, size = 1, stroke = 0.5, show.legend = FALSE) +
  geom_smooth(span = 1.2, show.legend = FALSE, se = FALSE, size = 0.8) +
  scale_y_continuous(name = "Spearman rank correlation (Hydrology)", 
                     breaks = seq(-0.8, 0.7, by = 0.2),
                     labels = c("-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6"), limits = c(-0.8, 0.7)) +
  scale_x_continuous(name = "Number of Sites",
                     expand = c(0,0), limits = c(0.5, 5.6)) +
  scale_colour_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  scale_fill_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  general_theme +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
  annotate(geom = "text", x = 5.4, y = 0.7, label = "b", size = 2.7, fontface = 2)
panelB

# Nicholas’ tweak ---------------------------------------------------------------------------------------------------------------------

figleg <- ggplot() +
  coord_cartesian(ylim = c(0, 1.95), xlim = c(-5, 7)) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(t = 5.5, r = 0, b = 10, l = -20, unit = "pt")) +
  annotate(geom = "text", x = -0.8, y = 1.9, label = "Invertebrate functional groups", size = 3, fontface = 2) +
  annotate(geom = "segment", x = -2, xend = -1, y = 1.55, yend = 1.55, size = 1, colour = "#007EB1") +
  annotate(geom = "segment", x = -2, xend = -1, y = 1.3, yend = 1.3, size = 1, colour = "#00888C") +
  annotate(geom = "segment", x = -2, xend = -1, y = 1.05, yend = 1.05, size = 1, colour = "#008859") +
  annotate(geom = "segment", x = 1.2, xend = 2.2, y = 1.55, yend = 1.55, size = 1, colour = "#607700") +
  annotate(geom = "segment", x = 1.2, xend = 2.2, y = 1.3, yend = 1.3, size = 1, colour = "#A95700") +
  annotate(geom = "segment", x = 1.2, xend = 2.2, y = 1.05, yend = 1.05, size = 1, colour = "#B94363") +
  annotate(geom = "text", x = -3.1, y = 1.55, label = "Engulfer", size = 3) +
  annotate(geom = "text", x = -3.3, y = 1.3, label = "Filter feeder", size = 3) +
  annotate(geom = "text", x = -3.1, y = 1.05, label = "Gatherer", size = 3) +
  annotate(geom = "text", x = 0.2, y = 1.55, label = "Piercer", size = 3) +
  annotate(geom = "text", x = 0.2, y = 1.3, label = "Scraper", size = 3) +
  annotate(geom = "text", x = 0.2, y = 1.05, label = "Shredder", size = 3) +
  annotate(geom = "text", x = 5.2, y = 1.9, label = "Stocks and fluxes", size = 3, fontface = 2) +
  annotate(geom = "segment", x = 6.5, xend = 7.5, y = 1.55, yend = 1.55, size = 1, colour = "#4D67C3") +
  annotate(geom = "segment", x = 6.5, xend = 7.5, y = 1.3, yend = 1.3, size = 1, colour = "#008200") +
  annotate(geom = "segment", x = 6.5, xend = 7.5, y = 1.05, yend = 1.05, size = 1, colour = "#8C6900") +
  annotate(geom = "text", x = 4.5, y = 1.55, label = "Bacterial density", size = 3) +
  annotate(geom = "text", x = 4.5, y = 1.3, label = "CO[2]~concentration", size = 3, parse = TRUE) +
  annotate(geom = "text", x = 4.5, y = 1.05, label = "Nitrogen uptake", size = 3)
figleg

# saving plots ------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_2_Nicholas_tweak_opt3.png", 
       plot = grid.arrange(panelA, panelB, figleg, nrow = 3, heights = c(1, 1, 0.5)),
       width = 89, height = 183, units = "mm", dpi = 600)

#here, need to have "Number of sites under panel a and B only, and to have everything fit

# ggsave(filename = "figure/figures/Figure_2.pdf", 
#        plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB), ncol = 2), 
#                            bottom = textGrob("Number of Sites", 
#                                              gp = gpar(fontsize = 9, font = 2), vjust = 0)),
#        width = 183, height = 89, units = "mm", dpi = 300)
# 
# ggsave(filename = "figure/figures/Figure_2.png", 
#        plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB, fig_legend1, fig_legend2), ncol = 2), 
#                            bottom = textGrob("Number of Sites", 
#                                              gp = gpar(fontsize = 9, font = 2), vjust = 0)),
#        width = 183, height = 89, units = "mm", dpi = 600)
