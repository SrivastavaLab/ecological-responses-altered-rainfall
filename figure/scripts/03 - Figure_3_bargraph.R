# clean up the environment ------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)

# load data ---------------------------------------------------------------

figure_data <- readRDS("figure/data/figure04_summary_list.rds")

# general theme for ggplot2 -----------------------------------------------

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title = element_text(face = "bold", size = 10, colour = "black"),
                       axis.text = element_text(size = 8, colour = "black"), 
                       axis.line.x = element_line(colour = "black", size = 0.25),
                       axis.line.y = element_line(colour = "black", size = 0.25),
                       panel.grid = element_blank(),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 10, 0, 0)))

# color scheme ------------------------------------------------------------

## fixing the order of the levels
                     
figure_data$aic.percent$response<-c("Decomposition", "Nitrogen uptake", "CO2 concentration", "Shredder", "Filter feeder", "Scraper", "Gatherer",           
                                    "Engulfer", "Piercer", "Bacterial density","Total Invertebrates") 
View(figure_data$aic.percent)
  
figure_data$aic.percent$response <- forcats::fct_relevel(figure_data$aic.percent$response,
                                                         "CO2 concentration", "Decomposition", "Nitrogen uptake", 
                                                         "Bacterial density", "Total Invertebrates", "Engulfer",
                                                         "Piercer", "Shredder", "Gatherer", "Scraper", "Filter feeder")

# BAR GRAPH ---------------------------------------------------------------

## for rainfall models
fig3_panelA <- ggplot(data = figure_data$aic.percent, mapping = aes(x = response, y = draintrue, fill = raintype)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_bar(colour = "black", stat = "identity", size = 0.3, width = 0.6) +
  scale_y_continuous(name = "Effect of rainfall",
                     limits = c(-0.03, 0.19), breaks = seq(0, 0.2, 0.05), expand = c(0,0)) +
  scale_x_discrete("") +
  scale_fill_manual(values = c("purple3", "gold1", "white")) +
  general_theme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,4,1,-1), "lines"),
        axis.text.y = element_text(face = 2)) +
  annotate(geom = "segment", x = 3.5, xend = 3.5, y = -0.03, yend = 0.16, colour = "black", linetype = 2) +
  annotate(geom = "segment", x = 5.5, xend = 5.5, y = -0.03, yend = 0.16, colour = "black", linetype = 2) +
  coord_flip() +
  annotate(geom = "text", x = 11, y = 0.17, label = "a", size = 3, fontface = 2)
fig3_panelA

## for hydrological models
fig3_panelB <- ggplot(data = figure_data$aic.percent, mapping = aes(x = response, y = dhydrotrue, fill = raintype)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_bar(colour = "black", stat = "identity", size = 0.3, width = 0.6) +
  scale_y_continuous(name = "Effect of hydrology",
                     limits = c(-0.03, 0.33), breaks = seq(0, 0.25, 0.05), expand = c(0,0)) +
  scale_x_discrete("") +
  scale_fill_manual(values = c("purple3", "gold1", "white")) +
  general_theme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,-4), "lines"),
        axis.text.y = element_text(face = 2)) +
  annotate(geom = "segment", x = 3.5, xend = 3.5, y = -0.03, yend = 0.3, colour = "black", linetype = 2) +
  annotate(geom = "segment", x = 5.5, xend = 5.5, y = -0.03, yend = 0.3, colour = "black", linetype = 2) +
  coord_flip() +
  annotate(geom = "text", x = 11, y = 0.30, label = "b", size = 3, fontface = 2) +
  annotation_custom(grob = textGrob("Fluxes", rot = 270, 
                                    gp = gpar(fontface = 2, fontsize = 10, col = "black")), 
                    xmin = 2, xmax = 2, ymin = 0.34, ymax = 0.34) +
  annotation_custom(grob = textGrob("Stocks", rot = 270, 
                                    gp = gpar(fontface = 2, fontsize = 10, col = "black")), 
                    xmin = 4.5, xmax = 4.5, ymin = 0.34, ymax = 0.34) + 
  annotation_custom(grob = textGrob("Functional groups", rot = 270, 
                                    gp = gpar(fontface = 2, fontsize = 10, col = "black")), 
                    xmin = 8.5, xmax = 8.5, ymin = 0.34, ymax = 0.34)
fig3_panelB

gt <- ggplot_gtable(ggplot_build(fig3_panelB))
gt$layout$clip[gt$layout$name=="panel"] <- "off"

# saving figure -----------------------------------------------------------

ggsave(filename = "figure/figures/Figure_3.pdf",
       plot = grid.arrange(arrangeGrob(grobs = list(fig3_panelA, gt), ncol = 2, nrow = 1)),
       width = 183, height = 120, units = "mm", dpi = 300)

ggsave(filename = "figure/figures/Figure_3.png",
       plot = grid.arrange(arrangeGrob(grobs = list(fig3_panelA, gt), ncol = 2, nrow = 1)),
       width = 183, height = 120, units = "mm", dpi = 600)
