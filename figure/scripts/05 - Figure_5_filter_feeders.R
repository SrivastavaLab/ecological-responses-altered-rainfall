# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)

# load data ---------------------------------------------------------------------------------------------------------------------------

figure5ab_data <- readRDS("figure/data/figure04_summary_list.rds") # for panels A and B
figure5c_data<-readRDS("figure/data/figure05_hydro_list.rds") # for panel C

# separating the data -----------------------------------------------------------------------------------------------------------------

fig5a_data <- figure5ab_data$fig_ff_rain #this is Fig 5a
fig5b_data <- figure5ab_data$fig_ff_hydro #this is Fig 5b
fig5c_data <- figure5c_data$fig_cvdepth_rain #this is Fig 5c

# general theme for ggplot2 -----------------------------------------------------------------------------------------------------------

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title = element_text(face = "bold", size = 12, colour = "black"),
                       axis.text = element_text(size = 10, colour = "black"), 
                       axis.line.x = element_line(colour = "black", size = 0.5),
                       axis.line.y = element_line(colour = "black", size = 0.5),
                       panel.grid = element_blank(),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 5, 0, 0)))

# colors ------------------------------------------------------------------------------------------------------------------------------

# c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")
# ARGENTINA, CARDOSO, COLOMBIA, COSTA RICA, FRENCH GUIANA, MACAE, PUERTO RICO

# panel A -----------------------------------------------------------------------------------------------------------------------------

panelA <- ggplot() +
  geom_point(data = fig5a_data$res, mapping = aes(x = mu.scalar, y = exp(visregRes)/100, fill = site),
             colour = "black", shape = 21, size = 1.5, alpha = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(data = fig5a_data$fit, mapping = aes(x = mu.scalar, y = exp(visregFit)/100, colour = site), size = 0.8) +
  scale_x_continuous(name = expression(bold(paste("Proportional change in ",{mu}))),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3),
                     labels = c("", "0.1", "0.2", "0.4", "", "0.6", "", "1", "1.5", "2", "", "3"), trans = "log") +
  scale_y_continuous(name = expression(bold(paste("Filter feeder biomass (mg plant"^-1,")"))), 
                     breaks = c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8), 
                     labels = c("1e-4", "1e-3", "1e-2", "0.1", "0.5", "1.0", "2.0", "4", "8"), trans = "log") +
  scale_colour_manual(values = c("mediumpurple4", "darkorange2", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "a", x = 3, y = 8, fontface = 2, size = 4, col = "black")
panelA


# panel B -----------------------------------------------------------------------------------------------------------------------------

panelB <- ggplot() +
  geom_point(data = fig5b_data$res, mapping = aes(x = cv.depth, y = exp(visregRes)/100, fill = site),
             colour = "black", shape = 21, size = 1.5, alpha = 0.3) +
  geom_line(data = fig5b_data$fit, mapping = aes(x = cv.depth, y = exp(visregFit)/100, colour = site), size = 0.8) +
  scale_x_continuous(name = "CV water depth", breaks = seq(from = 0, to = 260, by = 30)) +
  scale_y_continuous(name = expression(bold(paste("Filter feeder biomass (mg plant"^-1,")"))),
                     breaks = c(1e-4, 1e-2, 0.3, 1, 3, 5), 
                     labels = c("1e-4", "1e-2", "0.3", "1.0", "3.0", "5.0"), limits = c(0.00009, 5), trans = "log") +
  scale_colour_manual(values = c("mediumpurple4", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("mediumpurple4", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none") +
  annotate(geom = "text", label = "b", x = 220, y = 5, fontface = 2, size = 4, col = "black")
panelB

# panel C -----------------------------------------------------------------------------------------------------------------------------

panelC <- ggplot() +
  geom_point(data = fig5c_data$res, mapping = aes(x = mu.scalar, y = exp(visregRes), fill = site),
             colour = "black", shape = 21, size = 1.5, alpha = 0.3) +
  geom_line(data = fig5c_data$fit, mapping = aes(x = mu.scalar, y = exp(visregFit), colour = site), size = 0.8) +
  scale_x_continuous(name = expression(bold(paste("Proportional change in ",{mu}))),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3),
                     labels = c("", "0.1", "0.2", "0.4", "", "0.6", "", "1", "1.5", "2", "", "3"), trans = "log") +
  scale_y_continuous(name = "CV water depth", breaks = seq(from = 0, to = 300, by = 50)) +
  scale_colour_manual(values = c("mediumpurple4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("mediumpurple4", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(0, 16, 0, 0), vjust = -0.2)) +
  annotate(geom = "text", label = "c", x = 3, y = 300, fontface = 2, size = 4, col = "black")
panelC

# panel D -----------------------------------------------------------------------------------------------------------------------------


fig_legend <- data.frame(
  sites = c("Argentina", "Cardoso", "Colombia", "Costa Rica", "French Guiana", "Macae", "Puerto Rico"), 
  sites_ordered = letters[1:7],
  x = rep(0, 7),
  y = rep(1, 7)) %>% 
  ggplot(mapping = aes(x = sites_ordered, y = x)) +
  geom_errorbar(mapping = aes(ymin = x - y, ymax = x + y, colour = sites_ordered), width = 0, size = 0.8) +
  scale_colour_manual(values = rev(c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1"))) +
  scale_x_discrete(labels = rev(c("Argentina", "Cardoso", "Colombia", "Costa Rica", "French Guiana", "Macae", "Puerto Rico"))) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  coord_flip() +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(x = c(50, 50, 50, 50), units = "pt"))
fig_legend

# saving figures ----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_5.pdf", 
       plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB, panelC, fig_legend), ncol = 2)),
       width = 183, height = 183, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Figure_5.png", 
       plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB, panelC, fig_legend), ncol = 2)),
       width = 183, height = 183, units = "mm", dpi = 600)
