
# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)

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
  theme(plot.margin = unit(c(10, 5.5, -10, 5.5), "pt")) +
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
  scale_x_continuous(name = "", expand = c(0,0), limits = c(0.5, 5.6)) +
  scale_colour_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  scale_fill_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  general_theme +
  theme(plot.margin = unit(c(10, 5.5, -10, 5.5), "pt")) +
  annotate(geom = "text", x = 5.4, y = 0.7, label = "b", size = 2.7, fontface = 2)
panelB


# panel C and D -----------------------------------------------------------------------------------------------------------------------------


#needs to have a label saying "Invertebrate functional groups" above it
# fig_legend1 <- data.frame(
#   sites = c("engulfer_bio", "filter.feeder_bio", "gatherer_bio", "piercer_bio",  "scraper_bio", "shredder_bio"),
#   sites_ordered = letters[1:6],
#   x = rep(0, 6),
#   y = rep(1, 6)) %>%
#   ggplot(mapping = aes(x = sites_ordered, y = x)) +
#   geom_errorbar(mapping = aes(ymin = x - y, ymax = x + y, colour = sites_ordered), width = 0, size = 0.8) +
#   scale_colour_manual(values = rev(c("#007EB1","#00888C","#008859","#607700","#A95700","#B94363"))) +
#   scale_x_discrete(labels = rev(c("Engulfer", "Filter feeder", "Gatherer","Piercer", "Scraper", "Shredder"))) +
#   scale_y_continuous(limits = c(-1.1, 1.1)) +
#   coord_flip() +
#   theme(panel.background = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
#         axis.title = element_blank(),
#         panel.grid = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = "none",
#         plot.margin = unit(x = c(40, 40, 0, 40), units = "pt"))
# fig_legend1

# Nicholas’ tweak ---------------------------------------------------------------------------------------------------------------------

#needs to have a label saying "Invertebrate functional groups" above it
fig_legend1a <- data.frame(
  sites = c("engulfer_bio", "filter.feeder_bio", "gatherer_bio"), 
  sites_ordered = letters[1:3],
  x = rep(0, 3),
  y = rep(1, 3)) %>% 
  ggplot(mapping = aes(x = sites_ordered, y = x)) +
  geom_errorbar(mapping = aes(ymin = x - y, ymax = x + y, colour = sites_ordered), width = 0, size = 0.5) +
  scale_colour_manual(values = rev(c("#007EB1","#00888C","#008859"))) +
  scale_x_discrete(labels = rev(c("Engulfer", "Filter feeder", "Gatherer"))) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(x = c(0, -10, 0, 40), units = "pt"))
fig_legend1a

#needs to have a label saying "Invertebrate functional groups" above it
fig_legend1b <- data.frame(
  sites = c("piercer_bio",  "scraper_bio", "shredder_bio"), 
  sites_ordered = letters[1:3],
  x = rep(0, 3),
  y = rep(1, 3)) %>% 
  ggplot(mapping = aes(x = sites_ordered, y = x)) +
  geom_errorbar(mapping = aes(ymin = x - y, ymax = x + y, colour = sites_ordered), width = 0, size = 0.5) +
  scale_colour_manual(values = rev(c("#607700","#A95700","#B94363"))) +
  scale_x_discrete(labels = rev(c("Piercer", "Scraper", "Shredder"))) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(x = c(0, 0, 0, 50), units = "pt"))
fig_legend1b

#needs to have a line saying "Stocks and fluxes" above it:
fig_legend2 <- data.frame(
  sites = c("bacteria.per.nl.final", "log.co2.final", "scaled.n15.bromeliad.final"), 
  sites_ordered = letters[1:3],
  x = rep(0, 3),
  y = rep(1, 3)) %>% 
  ggplot(mapping = aes(x = sites_ordered, y = x)) +
  geom_errorbar(mapping = aes(ymin = x - y, ymax = x + y, colour = sites_ordered), width = 0, size = 0.5) +
  scale_colour_manual(values = rev(c("#4D67C3","#008200","#8C6900"))) +
  scale_x_discrete(labels = rev(c("Bacterial density", expression(paste("CO"[2]," concentration")),  "Nitrogen uptake"))) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(x = c(0, 75, 0, 75), units = "pt"))
fig_legend2

# saving plots ------------------------------------------------------------------------------------------
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

## saving final plot with a big tweak
ggsave(filename = "figure/figures/Figure_2.pdf", 
       plot = grid.arrange(cowplot::plot_grid(arrangeGrob(fig_legend1a, fig_legend1b, ncol = 2, top = textGrob("Invertebrate functional group", gp = gpar(fontsize = 8, font = 2), hjust = 0.3)),
                                              arrangeGrob(fig_legend2, top = textGrob("Stocks and fluxes", gp = gpar(fontsize = 8, font = 2))),
                                              panelA, panelB, rel_heights = c(0.5, 2)), 
                           bottom = textGrob("Number of Sites", gp = gpar(fontsize = 9, font = 2), vjust = 0.5)),
       width = 183, height = 89, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Figure_2.png", 
       plot = grid.arrange(cowplot::plot_grid(arrangeGrob(fig_legend1a, fig_legend1b, ncol = 2, top = textGrob("Invertebrate functional group", gp = gpar(fontsize = 8, font = 2), hjust = 0.3)),
                                              arrangeGrob(fig_legend2, top = textGrob("Stocks and fluxes", gp = gpar(fontsize = 8, font = 2))),
                                              panelA, panelB, rel_heights = c(0.5, 2)), 
                           bottom = textGrob("Number of Sites", gp = gpar(fontsize = 9, font = 2), vjust = 0.5)),
       width = 183, height = 89, units = "mm", dpi = 600)
