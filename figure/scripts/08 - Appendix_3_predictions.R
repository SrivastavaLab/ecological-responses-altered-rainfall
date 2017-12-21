
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

# panel A -----------------------------------------------------------------------------------------------

panelA <- ggplot(data = concord.rain, mapping = aes(x = Sites, y = CCC, colour = Response)) +
  annotate(geom = "rect", xmin = 0.5, xmax = 5.5, ymin = -0.29, ymax = 0.29, alpha = 0.1, colour = "grey90", size = 0) +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.3, size = 1, stroke = 0.5, show.legend = FALSE) +
  geom_smooth(span = 1.2, show.legend = FALSE, se = FALSE, size = 0.8) +
  scale_y_continuous(name = "Concordance correlation coefficient (rainfall)", 
                     breaks = seq(-0.8, 0.6, by = 0.2),
                     labels = c("-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6"), limits = c(-0.8, 0.6)) +
  scale_x_continuous(name = "", expand = c(0,0), limits = c(0.5, 5.6)) +
  general_theme +
  theme(plot.margin = unit(c(5.5, 5.5, -10, 5.5), "pt")) +
  annotate(geom = "text", x = 5.4, y = 0.6, label = "a", size = 2.7, fontface = 2)
panelA

# panel B -----------------------------------------------------------------------------------------------

panelB <- ggplot(data = concord.hydro, mapping = aes(x = Sites, y = CCC, colour = Response)) +
  annotate(geom = "rect", xmin = 0.5, xmax = 5.5, ymin = -0.29, ymax = 0.29, alpha = 0.1, colour = "grey90", size = 0) +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.3, size = 1, stroke = 0.5, show.legend = FALSE) +
  geom_smooth(span = 1.2, show.legend = FALSE, se = FALSE, size = 0.8) +
  scale_y_continuous(name = "Concordance correlation coefficient (hydrology)", 
                     breaks = seq(-0.8, 0.6, by = 0.2),
                     labels = c("-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6"), limits = c(-0.8, 0.6)) +
  scale_x_continuous(name = "", expand = c(0,0), limits = c(0.5, 5.6)) +
  general_theme +
  theme(plot.margin = unit(c(5.5, 5.5, -10, 5.5), "pt")) +
  annotate(geom = "text", x = 5.4, y = 0.6, label = "b", size = 2.7, fontface = 2)
panelB

# saving plot -------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Appendix3.pdf", 
       plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB), ncol = 2), 
                           bottom = textGrob("Number of Sites", 
                                             gp = gpar(fontsize = 9, font = 2), vjust = 0)),
       width = 183, height = 89, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Appendix3.png", 
       plot = grid.arrange(arrangeGrob(grobs = list(panelA, panelB), ncol = 2), 
                           bottom = textGrob("Number of Sites", 
                                             gp = gpar(fontsize = 9, font = 2), vjust = 0)),
       width = 183, height = 89, units = "mm", dpi = 600)
