# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)
library(scales)

# load data ---------------------------------------------------------------------------------------------------------------------------

#this to be a scatterplot of 7 site-level mean correlations, with a darker point of
table.corr.fg <- read.csv("Data/tablecorrfg.csv") %>% 
  rename(site = X, r = ".") %>% #mean of site means: 0.029, 95% CI is +/-0.038
  mutate(site = c("Argentina", "Cardoso", "Colombia", "Costa Rica", "French Guiana", "Macae", "Puerto Rico"),
         null = rep(1, 7))

overall <- data.frame(mean_calc = 0.029, ci_calc = 0.038, x_position = 1.05)

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

# creating plot -----------------------------------------------------------------------------------------------------------------------

app6 <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5, size = 0.5) +
  geom_point(data = table.corr.fg , mapping = aes(x = null, y = r, fill = site),
             size = 3, shape = 21, colour = "black", stroke = 0.4, show.legend = FALSE) +
  geom_text(data = table.corr.fg , mapping = aes(x = null - 0.04, y = r, label = site), size = 2.8) +
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = "Site", limits = c(0.9, 1.12), expand = c(0,0)) +
  scale_y_continuous(name = expression(paste("Mean correlation between functional groups (", italic("r"), ")")), limits=c(-0.06, 0.11), expand = c(0,0)) +
  general_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  geom_errorbar(data = overall, mapping = aes(x = x_position, ymin = mean_calc - ci_calc, ymax = mean_calc + ci_calc), size = 0.5) +
  geom_point(data = overall , mapping = aes(x = x_position, y = mean_calc),
             size = 5, shape = 21, colour = "black", fill = "grey60", stroke = 0.8, show.legend = FALSE)
app6

# data for triangular matrix ----------------------------------------------------------------------------------------------------------

fgpairs <- read.csv("Data/fgpairs.csv") %>% 
  rename(site = X, correl = V3) # this to be a triangular matrix of values (correl = mean correlations)
fgpairs

# saving figure -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Appendix_6_correlation_funct_groups.pdf",
       plot = app6, width = 89, height = 89, dpi = 300, units = "mm")

ggsave(filename = "figure/figures/Appendix_6_correlation_funct_groups.png",
       plot = app6, width = 89, height = 89, dpi = 600, units = "mm")