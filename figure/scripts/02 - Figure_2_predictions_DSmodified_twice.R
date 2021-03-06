
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
  scale_x_continuous(name = "Number of Sites",
                     expand = c(0,0), limits = c(0.5, 5.6)) +
  scale_colour_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  scale_fill_manual(values = c("#4D67C3","#007EB1","#00888C","#008859","#008200","#607700","#8C6900","#A95700","#B94363")) +
  general_theme +
  theme(plot.margin = unit(c(-10, 5.5, 10, 5.5), "pt")) +
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

figleg <- ggplot() +
  coord_cartesian(ylim = c(1, 1.95), xlim = c(-5, 4.5)) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(t = 50, r = 20, b = 30, l = -10, unit = "pt")) +
  annotate(geom = "text", x = 0, y = 1.9, label = "Invertebrate functional groups", size = 3, fontface = 2) +
  annotate(geom = "segment", x = -1.3, xend = -0.3, y = 1.8, yend = 1.8, size = 1, colour = "#007EB1") +
  annotate(geom = "segment", x = -1.3, xend = -0.3, y = 1.7, yend = 1.7, size = 1, colour = "#00888C") +
  annotate(geom = "segment", x = -1.3, xend = -0.3, y = 1.6, yend = 1.6, size = 1, colour = "#008859") +
  annotate(geom = "segment", x = 3, xend = 4, y = 1.8, yend = 1.8, size = 1, colour = "#607700") +
  annotate(geom = "segment", x = 3, xend = 4, y = 1.7, yend = 1.7, size = 1, colour = "#A95700") +
  annotate(geom = "segment", x = 3, xend = 4, y = 1.6, yend = 1.6, size = 1, colour = "#B94363") +
  annotate(geom = "text", x = -2.8, y = 1.8, label = "Engulfer", size = 3) +
  annotate(geom = "text", x = -3, y = 1.7, label = "Filter feeder", size = 3) +
  annotate(geom = "text", x = -2.8, y = 1.6, label = "Gatherer", size = 3) +
  annotate(geom = "text", x = 1.5, y = 1.8, label = "Piercer", size = 3) +
  annotate(geom = "text", x = 1.5, y = 1.7, label = "Scraper", size = 3) +
  annotate(geom = "text", x = 1.5, y = 1.6, label = "Shredder", size = 3) +
  annotate(geom = "text", x = 0, y = 1.3, label = "Stocks and fluxes", size = 3, fontface = 2) +
  annotate(geom = "segment", x = -1.5, xend = -0.5, y = 1.2, yend = 1.2, size = 1, colour = "#4D67C3") +
  annotate(geom = "segment", x = 3.5, xend = 4.5, y = 1.2, yend = 1.2, size = 1, colour = "#008200") +
  annotate(geom = "segment", x = 0.8, xend = 1.8, y = 1.1, yend = 1.1, size = 1, colour = "#8C6900") +
  annotate(geom = "text", x = -3.5, y = 1.2, label = "Bacterial density", size = 3) +
  annotate(geom = "text", x = 1.5, y = 1.2, label = "CO[2]~concentration", size = 3, parse = TRUE) +
  annotate(geom = "text", x = -1.3, y = 1.1, label = "Nitrogen uptake", size = 3)
figleg

# saving plots ------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_2_Nicholas_tweak_opt1.png", 
       plot = grid.arrange(panelA, figleg, panelB, ncol = 2, widths = c(0.8, 0.7)),
       width = 183, height = 183, units = "mm", dpi = 600)

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
