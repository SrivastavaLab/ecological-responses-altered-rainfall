# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)

# load data ---------------------------------------------------------------------------------------------------------------------------

family_sensitivity <- read_csv("Data/family_sensitivity.csv")

# general theme for ggplot2 -----------------------------------------------------------------------------------------------------------

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title = element_text(face = "bold", size = 9.4, colour = "black"),
                       axis.text = element_text(size = 8, colour = "black"), 
                       axis.line.x = element_line(colour = "black", size = 0.3),
                       axis.line.y = element_line(colour = "black", size = 0.3),
                       panel.grid = element_blank(),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 8, 0, 0)),
                       plot.margin = unit(x = c(5.5, 9.5, 5.5, 5.5), units = "pt"))

# colors ------------------------------------------------------------------------------------------------------------------------------

# c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")
# ARGENTINA, CARDOSO, COLOMBIA, COSTA RICA, FRENCH GUIANA, MACAE, PUERTO RICO

# fix data ----------------------------------------------------------------------------------------------------------------------------

family_sensitivity$site<-c("Argentina", "Cardoso", "Colombia", "Costa Rica", "French Guiana", "Macae", "Puerto Rico")
family_sensitivity$textx <- family_sensitivity$ambient.driedout
family_sensitivity$texty <- family_sensitivity$sens.index
family_sensitivity$textx[1] <- family_sensitivity$textx[1] + 0.03
family_sensitivity$textx[3] <- family_sensitivity$textx[3] + 0.03
family_sensitivity$textx[4] <- family_sensitivity$textx[4] - 0.015
family_sensitivity$textx[5] <- family_sensitivity$textx[5] + 0.042
family_sensitivity$textx[6] <- family_sensitivity$textx[6] + 0.01
family_sensitivity$textx[7] <- family_sensitivity$textx[7] + 0.036
family_sensitivity$texty[4] <- family_sensitivity$texty[4] + 0.007
family_sensitivity$texty[6] <- family_sensitivity$texty[6] - 0.005

# figure ------------------------------------------------------------------------------------------------------------------------------

fig4 <- ggplot() +
  geom_smooth(data = family_sensitivity, mapping = aes(x = ambient.driedout, y = sens.index),
              method = lm, se = FALSE, size = 0.5, colour = "black") +
  geom_point(data = family_sensitivity, mapping = aes(x = ambient.driedout, y = sens.index, fill = site),
             size = 2, shape = 21, colour = "black", stroke = 0.2, show.legend = FALSE) +
  geom_text(data = family_sensitivity, mapping = aes(x = textx, y = texty, label = site), size = 2.8) +
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "forestgreen",
                               "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = "Ambient proportion of days without water") +
  scale_y_continuous(name = "Species pool sensitivity to rain") +
  general_theme
fig4

# saving figure -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_4.pdf", plot = fig4, width = 89, height = 89, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Figure_4.png", plot = fig4, width = 89, height = 89, units = "mm", dpi = 600)

