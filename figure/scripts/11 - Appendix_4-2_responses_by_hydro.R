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

# CV WATER DEPTH vs DECOMPOSITION -----------------------------------------------------------------------------------------------------

panelA <- ggplot() +
  #  geom_point(data = figure_data$fig_decomp_hydro$res, mapping = aes(x = cv.depth, y = visregRes^2, fill = site),
  #             colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_decomp_hydro$fit, mapping = aes(x = cv.depth, y = visregFit^2, colour = site), size = 0.8) +
  scale_x_continuous(name = "CV water depth", breaks = seq(from = 0, to = 260, by = 30)) +
  scale_y_continuous(name = "Litter decomposition (% mass loss)", 
                     breaks = c(0.01, 0.1, 0.2, 0.4, 0.6, 0.8), 
                     labels = c("0.01", "0.1", "0.2", "0.4", "0.6", "0.8"), limits = c(0.01, 0.8), trans = "sqrt") +
  scale_colour_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "a", x = 220, y = 0.8, fontface = 2, size = 4, col = "black")
panelA

# CV WATER DEPTH vs FILTER FEEDERS ----------------------------------------------------------------------------------------------------

panelB <- ggplot() +
  # geom_point(data = figure_data$fig_ff_hydro$res, mapping = aes(x = cv.depth, y = exp(visregRes)/100, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_ff_hydro$fit, mapping = aes(x = cv.depth, y = exp(visregFit)/100, colour = site), size = 0.8) +
  scale_x_continuous(name = "CV water depth", breaks = seq(from = 0, to = 260, by = 30)) +
  scale_y_continuous(name = expression(bold(paste("Filter feeder biomass (mg plant"^-1,")"))),
                     breaks = c(1e-4, 1e-2, 0.3, 1, 3, 5), 
                     labels = c("1e-4", "1e-2", "0.3", "1.0", "3.0", "5.0"), limits = c(0.00009, 7), trans = "log") +
  scale_colour_manual(values = c("plum4", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "b", x = 220, y = 5, fontface = 2, size = 4, col = "black")
panelB

# LONGEST DROUGHT PERIOD vs NITROGEN UPTAKE -------------------------------------------------------------------------------------------

panelC <- ggplot() +
  # geom_point(data = figure_data$fig_N_hydro$res, mapping = aes(x = long_dry, y = (visregRes^(1/0.125)) - 4, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_N_hydro$fit, mapping = aes(x = long_dry, y = (visregFit^(1/0.125)) - 4, colour = site), size = 0.8) +
  scale_x_continuous(name = "Longest drought period", breaks = seq(0, 40, 10)) +
  scale_y_continuous(name = expression(bold(paste("Detrital nitrogen uptake (",{delta}^15*N,")"))), 
                     trans = trans_new("eight", function(x) (x + 4)^0.125, function(x) (x^8) - 4, 
                                       breaks = extended_breaks(n = 5), 
                                       format = scientific_format(digits = 0)),
                     breaks = c(-9e-16, 1e1, 1e2, 5e2, 1e3, 2e3),
                     labels = c("-9e-16", "1e1", "1e2", "5e2", "1e3", "2e3")) +
  scale_colour_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "c", x = 40, y = 2.1e3, fontface = 2, size = 4, col = "black")
panelC

# PROPORTION OF DRIEDOUT vs SCRAPER BIOMASS -------------------------------------------------------------------------------------------

panelD <- ggplot() +
  # geom_point(data = figure_data$fig_sc_hydro$res, mapping = aes(x = prop.driedout.days, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_sc_hydro$fit, mapping = aes(x = prop.driedout.days, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportion of drought days", breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(name = expression(bold(paste("Scraper biomass (mg plant"^-1,")"))), 
                     breaks = c(0.01, 0.1, 1, 2, 4, 8, 16, 21),
                     labels = c("0.01", "0.1", "1.0", "2.0", "4.0", "8.0", "16.0", "21.0"),
                     limits = c(0.01, 21), trans = "log") +
  scale_colour_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "d", x = 0.8, y = 21, fontface = 2, size = 4, col = "black")
panelD

# PROPORTION OF DRIEDOUT DAYS vs GATHERER BIOMASS -------------------------------------------------------------------------------------

panelE <- ggplot() +
  # geom_point(data = figure_data$fig_ga_hydro$res, mapping = aes(x = last_wet, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_ga_hydro$fit, mapping = aes(x = last_wet, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Days since last overflow", breaks = seq(from = 0, to = 60, by = 10)) +
  scale_y_continuous(name = expression(bold(paste("Gatherer biomass (mg plant"^-1,")"))),
                     breaks = c(0.1, 0.3, 1, 2, 4, 8, 16, 24), trans = "log") +
  scale_colour_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "e", x = 65, y = 24, fontface = 2, size = 4, col = "black")
panelE

# PROPORTION OF OVERFLOW DAYS vs SHREDDER BIOMASS -------------------------------------------------------------------------------------

panelF <- ggplot() +
  # geom_point(data = figure_data$fig_sh_hydro$res, mapping = aes(x = prop.overflow.days, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_sh_hydro$fit, mapping = aes(x = prop.overflow.days, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportion of overflow days", breaks = seq(from = 0, to = 0.8, by = 0.1)) +
  scale_y_continuous(name = expression(bold(paste("Shredder biomass (ln mg plant"^-1,")"))),
                     breaks = c(0.1, 0.3, 1, 2, 4, 8, 12), trans = "log") +
  scale_colour_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "f", x = 0.8, y = 12, fontface = 2, size = 4, col = "black")
panelF

# PROPORTION OF OVERFLOW DAYS vs PIERCER BIOMASS --------------------------------------------------------------------------------------

panelG <- ggplot() +
  # geom_point(data = figure_data$fig_pi_hydro$res, mapping = aes(x = prop.overflow.days, y = exp(visregRes)/10, fill = site),
  #            colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = figure_data$fig_pi_hydro$fit, mapping = aes(x = prop.overflow.days, y = exp(visregFit)/10, colour = site), size = 0.8) +
  scale_x_continuous(name = "Proportion of overflow days", breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(name = expression(bold(paste("Piercer biomass (mg plant"^-1,")"))),
                     breaks = c(1e-9, 5e-8, 2e-6, 2e-4, 2e-2, 1, 30), 
                     labels = c("1e-9", "5e-8", "2e-6", "2e-4", "2e-2", "1.0", "30.0"), trans = "log") +
  scale_colour_manual(values = c("plum4", "mediumblue", "lawngreen")) +
  scale_fill_manual(values = c("plum4", "mediumblue", "lawngreen")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "g", x = 0.8, y = 30, fontface = 2, size = 4, col = "black")
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

ggsave(filename = "figure/figures/Appendix_4-2_responses_by_hydro.pdf",
       plot = grid.arrange(panelA, panelB, panelC, panelD, panelE, panelF, panelG, fig_legend,ncol = 2),
       width = 210, height = 297, units = "mm", dpi = 300)


ggsave(filename = "figure/figures/Appendix_4-2_responses_by_rainfall.png",
       plot = grid.arrange(panelA, panelB, panelC, panelD, panelE, panelF, panelG, fig_legend, ncol = 2),
       width = 210, height = 297, units = "mm", dpi = 600)