## TO DO: FIX THE COLOR SCHEME THROUGHOUT PANELS, AS SOME SITES ARE MISSING AND ARE AUTOMATICALLY COLOURED BY THE VECTOR

# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)
library(scales)

# load data ---------------------------------------------------------------------------------------------------------------------------

figure_data <- readRDS("figure/data/figure04_summary_list.rds")

# general theme for ggplot2 -----------------------------------------------------------------------------------------------------------

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title = element_text(face = "bold", size = 10, colour = "black"),
                       axis.text = element_text(size = 8, colour = "black"), 
                       axis.ticks = element_line(colour = "black", size = 0.35),
                       axis.line.x = element_line(colour = "black", size = 0.35),
                       axis.line.y = element_line(colour = "black", size = 0.35),
                       panel.grid = element_blank(),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 10, 0, 0)))

# colors ------------------------------------------------------------------------------------------------------------------------------

# c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")
# ARGENTINA, CARDOSO, COLOMBIA, COSTA RICA, FRENCH GUIANA, MACAE, PUERTO RICO

## also
color_scheme <- data.frame(
  sites = c("Argentina", "Cardoso", "Colombia", "Costa Rica", "French Guiana", "Macae", "Puerto Rico"), 
  sites_ordered = letters[1:7],
  x = rep(0, 7),
  y = rep(1, 7))
color_scheme

# K SCALAR vs NITROGEN UPTAKE ---------------------------------------------------------------------------------------------------------

panelA <- ggplot() +
#  geom_point(data = figure_data$fig_N_rain$res, mapping = aes(x = k.scalar, y = (visregRes^(1/0.125)) - 4, fill = site),
#             colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_N_rain$fit, mapping = aes(x = k.scalar, y = (visregFit^(1/0.125)) - 4, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportional change in k", breaks = c(0.5, 1, 2), labels = c(0.5, 1, 2), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Detrital nitrogen uptake (",{delta}^15*N,")"))), 
                     trans = trans_new("eight", function(x) (x + 4)^0.125, function(x) (x^8) - 4, 
                                       breaks = extended_breaks(n = 5), 
                                       format = scientific_format(digits = 0)),
                     breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 4e7),
                     labels = c("1e2", "1e3", "1e4", "1e5", "1e6", "1e7", "4e7")) +
  scale_colour_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "a", x = 2, y = 4e7, fontface = 2, size = 4, col = "black")
panelA

# K SCALAR vs SHREDDER BIOMASS ----------------------------------------------------------------------------------------------------------

panelB <- ggplot() +
# geom_point(data = figure_data$fig_sh_rain$res, mapping = aes(x = k.scalar, y = exp(visregRes)/10, fill = site),
#            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_sh_rain$fit, mapping = aes(x = k.scalar, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportional change in k", breaks = c(0.5, 1, 2), labels = c(0.5, 1, 2), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Shredder biomass (mg plant"^-1,")"))),
                     breaks = c(0.1, 0.5, 1, 2, 4, 6, 8, 10),
                     trans = "log") +
  scale_colour_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "b", x = 2, y = 11, fontface = 2, size = 4, col = "black")
panelB

# K SCALAR vs GATHERER BIOMASS --------------------------------------------------------------------------------------------------------

panelC <- ggplot() +
  # geom_point(data = figure_data$fig_ga_rain$res, mapping = aes(x = k.scalar, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_ga_rain$fit, mapping = aes(x = k.scalar, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportional change in k", breaks = c(0.5, 1, 2), labels = c(0.5, 1, 2), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Gatherer biomass (mg plant"^-1,")"))), 
                     breaks = c(0.1, 0.3, 1, 2, 4, 8, 16, 32),
                     trans = "log") +
  scale_colour_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "c", x = 2, y = 32, fontface = 2, size = 4, col = "black")
panelC

# K SCALAR vs SCRAPER BIOMASS --------------------------------------------------------------------------------------------------------

panelD <- ggplot() +
  # geom_point(data = figure_data$fig_sc_rain$res, mapping = aes(x = k.scalar, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_sc_rain$fit, mapping = aes(x = k.scalar, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportional change in k", breaks = c(0.5, 1, 2), labels = c(0.5, 1, 2), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Scraper biomass (mg plant"^-1,")"))),
                     breaks = c(0.5, 1, 2, 4, 8, 16, 24),
                     limits = c(0.5, 24), trans = "log") +
  scale_colour_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "d", x = 2, y = 24, fontface = 2, size = 4, col = "black")
panelD

# MU SCALAR vs FILTER FEEDERS ---------------------------------------------------------------------------------------------------------

panelE <- ggplot() +
  # geom_point(data = figure_data$fig_ff_rain$res, mapping = aes(x = k.scalar, y = exp(visregRes)/100, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_ff_rain$fit, mapping = aes(x = mu.scalar, y = exp(visregFit)/100, colour = site), size = 0.8) +
  scale_x_continuous(name = expression(bold(paste("Proportional change in ",{mu}))),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3),
                     labels = c("", "0.1", "0.2", "0.4", "0.5", "0.6", "0.8", "1", "1.5", "2", "2.5", "3"), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Filter feeder biomass (mg plant"^-1,")"))), 
                     breaks = c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2), 
                     labels = c("1e-4", "1e-3", "1e-2", "0.1", "0.5", "1.0", "2.0"), trans = "log") +
  scale_colour_manual(values = c("plum4", "darkorange2", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "e", x = 3, y = 2, fontface = 2, size = 4, col = "black")
panelE

# MU SCALAR vs PIERCER BIOMASS --------------------------------------------------------------------------------------------------------

panelF <- ggplot() +
  # geom_point(data = figure_data$fig_ff_rain$res, mapping = aes(x = k.scalar, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_pi_rain$fit, mapping = aes(x = mu.scalar, y = exp(visregFit)/10, colour = site), size = 1) +
  scale_x_continuous(name = expression(bold(paste("Proportional change in ",{mu}))),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3),
                     labels = c("", "0.1", "0.2", "0.4", "0.5", "0.6", "0.8", "1", "1.5", "2", "2.5", "3"), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Piercer biomass (mg plant"^-1,")"))),
                     breaks = c(0.03, 0.1, 0.3, 0.9, 2, 4, 7, 10), limits = c(0.03, 10.1), trans = "log") +
  scale_colour_manual(values = c("plum4", "darkorange2", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "f", x = 3, y = 10, fontface = 2, size = 4, col = "black")
panelF

# MU SCALAR vs BACTERIAL DENSITY ------------------------------------------------------------------------------------------------------

panelG <- ggplot() +
  # geom_point(data = figure_data$fig_ba_rain$res, mapping = aes(x = k.scalar, y = exp(visregRes)/100, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_ba_rain$fit, mapping = aes(x = mu.scalar, y = exp(visregFit)/100, colour = site), size = 0.8) +
  scale_x_continuous(name = expression(bold(paste("Proportional change in ",{mu}))),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3),
                     labels = c("", "0.1", "0.2", "0.4", "0.5", "0.6", "0.8", "1", "1.5", "2", "2.5", "3"), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Bacterial density (cells ",{mu}*L^-1,")"))),
                     breaks = c(1, 2, 4, 8, 16, 32, 50), trans = "log") +
  scale_colour_manual(values = c("darkorange2","deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("darkorange2", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "g", x = 3, y = 50, fontface = 2, size = 4, col = "black")
panelG

# LEGEND ------------------------------------------------------------------------------------------------------------------------------

fig_legend <- ggplot(data = color_scheme, mapping = aes(x = sites_ordered, y = x)) +
  geom_errorbar(mapping = aes(ymin = x - y, ymax = x + y, colour = sites_ordered), width = 0, size = 1.5) +
  scale_colour_manual(values = rev(c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1"))) +
  scale_x_discrete(labels = rev(color_scheme$sites)) +
  scale_y_continuous(limits = c(-2, 2)) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
fig_legend

# SAVING PANELS -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Appendix_4-1_responses_by_rainfall.pdf",
       plot = grid.arrange(panelA, panelB, panelC, panelD, panelE, panelF, panelG, fig_legend, ncol = 2),
       width = 210, height = 297, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Appendix_4-1_responses_by_rainfall.png",
       plot = grid.arrange(panelA, panelB, panelC, panelD, panelE, panelF, panelG, fig_legend, ncol = 2),
       width = 210, height = 297, units = "mm", dpi = 600)
