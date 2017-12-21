# clean up the environment ------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# load packages -----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)
library(scales)

# loading data ------------------------------------------------------------------------------------------------------------------------

FlushedRetained <- read.csv("Data/FlushedRetained")
names(FlushedRetained) <- c("Source", "Chironomidae", "Tipulidae", "Pseudostigmatidae", "Psychodidae", "Scirtidae", "Ceratopogonidae", "Naididae", "Culicidae", "Syrphidae")

# tidying data ------------------------------------------------------------------------------------------------------------------------

longfr <- FlushedRetained %>% 
  gather(key = Family, value = Proportion, 2:10)
longfr

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

# creating figure ---------------------------------------------------------------------------------------------------------------------

app9 <- ggplot(longfr, aes(x = reorder(Family, -Proportion), y = Proportion, fill = Source)) + 
  geom_bar(position = "dodge", stat = "identity", colour = "black", size = 0.2) +
  scale_y_continuous(name = "Proportion of individuals", limits = c(0, 0.4), expand = c(0, 0)) +
  scale_x_discrete(name = "Family") +
  labs(fill = "") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top") +
  general_theme
app9

# saving figure -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Appendix_9_flushed_inverts.pdf", 
       plot = app9, width = 120, height = 89, dpi = 300, units = "mm")

ggsave(filename = "figure/figures/Appendix_9_flushed_inverts.png", 
       plot = app9, width = 120, height = 89, dpi = 600, units = "mm")
