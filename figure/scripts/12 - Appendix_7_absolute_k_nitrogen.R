# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)
library(scales)

# load data ---------------------------------------------------------------------------------------------------------------------------

Appendix7_data <- readRDS("figure/data/Appendix7_summary_list.rds")

# general theme for ggplot2 -----------------------------------------------------------------------------------------------------------

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title = element_text(face = "bold", size = 10, colour = "black"),
                       axis.text = element_text(size = 8, colour = "black"), 
                       axis.line.x = element_line(colour = "black", size = 0.2),
                       axis.line.y = element_line(colour = "black", size = 0.2),
                       panel.grid = element_blank(),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 3, 0, 0)))

# colors ------------------------------------------------------------------------------------------------------------------------------

# c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")
# ARGENTINA, CARDOSO, COLOMBIA, COSTA RICA, FRENCH GUIANA, MACAE, PUERTO RICO

# figure ------------------------------------------------------------------------------------------------------------------------------
#DS corrected code to have linear values of abs k placed on a log x axis

appendix7 <- ggplot() +
  geom_point(data = Appendix7_data$abs_k_Nfig$res, mapping = aes(x = abs(intended.k), y = (visregRes^(1/0.125)) - 4, fill = site),
             colour = "black", shape = 21, size = 1.5, alpha = 0.3) +
  geom_line(data = Appendix7_data$abs_k_Nfig$fit, mapping = aes(x = abs(intended.k), y = (visregFit^(1/0.125)) - 4, colour = site), size = 0.6) +
  scale_y_continuous(name = expression(bold(paste("Detrital nitrogen uptake (",{delta}^15*N,")"))), 
                     trans = trans_new("eight", function(x) (x + 4)^0.125, function(x) (x^8) - 4, 
                                       breaks = extended_breaks(n = 10), 
                                       format = scientific_format(digits = 0)),
                     breaks = c(1e-1, 2e1, 5e1, 1e2, 2e2),
                     labels = c("1e-1", "2e1", "5e1", "1e2", "2e2")) +
  scale_x_continuous(name = "Absolute rainfall dispersion (k)", breaks = c(0.03, 0.1, 0.3,1), labels = c(0.03, 0.1, 0.3, 1), trans = "log") +
  scale_colour_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none")
appendix7

# saving figure -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Appendix_7.pdf", plot = appendix7, width = 89, height = 89, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Appendix_7.png", plot = appendix7, width = 89, height = 89, units = "mm", dpi = 600)

