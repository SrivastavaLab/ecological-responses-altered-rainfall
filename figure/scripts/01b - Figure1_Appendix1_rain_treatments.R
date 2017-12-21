# cleaning up the environment ---------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# loading packages --------------------------------------------------------------------------------------------------------------------

library(tidyverse)

# loading data ------------------------------------------------------------------------------------------------------------------------

fulldata <- read.csv( "Data/fulldata.csv")
site_treatment_data <- read.csv( "Data/site_treatment_data.csv")

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

# creating functions ------------------------------------------------------------------------------------------------------------------

scaleFUN <- function(x) sprintf("%.2f", x)

# tidying site names ------------------------------------------------------------------------------------------------------------------

site_treatment_data$site <- sub(pattern = "argentina", replacement = "Argentina", x = site_treatment_data$site)
site_treatment_data$site <- sub(pattern = "cardoso", replacement = "Cardoso", x = site_treatment_data$site)
site_treatment_data$site <- sub(pattern = "colombia", replacement = "Colombia", x = site_treatment_data$site)
site_treatment_data$site <- sub(pattern = "costarica", replacement = "Costa Rica", x = site_treatment_data$site)
site_treatment_data$site <- sub(pattern = "frenchguiana", replacement = "French Guiana", x = site_treatment_data$site)
site_treatment_data$site <- sub(pattern = "macae", replacement = "Macae", x = site_treatment_data$site)
site_treatment_data$site <- sub(pattern = "puertorico", replacement = "Puerto Rico", x = site_treatment_data$site)

sitemean <- filter(fulldata,mu.scalar == "1", k.scalar == "1") %>%
  select(site, mu.scalar, k.scalar, intended.mu, intended.k)

sitecorners <- filter(fulldata,mu.scalar %in% c("0.1","3"),
                      k.scalar %in% c("0.5","2")) %>% 
  select(site, mu.scalar, k.scalar, intended.mu, intended.k)

sitecorners$order<-rep(c(1,2,4,3), 7)

sitecorners<-arrange(sitecorners, order)


# panel A -----------------------------------------------------------------------------------------------------------------------------

ggplot() +
  geom_point(data = fulldata, mapping = aes(x = intended.mu, y = intended.k, fill = site),
             colour = "grey80", shape = 21, size = 2, alpha = 0.5, na.rm = TRUE) + 
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = expression(bold(paste("Mean daily rainfall, ",mu, " (mm)"))), trans = "log2") +
  scale_y_continuous(name = "Rainfall dispersion, k", label = scaleFUN, trans = "log2") +
  coord_equal() +
  general_theme +
  theme(legend.position = "none")
#  geom_polygon(data = sitecorners, aes(color = site), fill = NA) + ## missing data = sitecorners

# panel B -----------------------------------------------------------------------------------------------------------------------------

ggplot(data = site_treatment_data) +
  facet_wrap(~ site) + 
  geom_point(mapping = aes(x = intended.mu, y = intended.k, fill = site),
             colour = "grey80", shape = 21, size = 2, alpha = 0.5, na.rm = TRUE) + 
  # missing data = sitecorners
  # geom_polygon(data = sitecorners, mapping = aes(color = site), fill = NA) +
  scale_fill_manual(values = c("plum4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = expression(bold(paste("Mean daily rainfall, ",mu, " (mm)"))), trans = "log2") +
  scale_y_continuous(name = "Rainfall dispersion, k", label = scaleFUN, trans = "log2") +
  coord_equal() +
  # missing data = sitemean
  # geom_point(data = sitemean, mapping = aes(color = site), size = 4, alpha = 1, na.rm = TRUE)+
  general_theme +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold", colour = "black", size = 8),
        legend.position = "none")
  
## mising data = original_data_parameters
# plot2<-plot1b + 
#   geom_point(pch = 21, size = 2.2, fill = "white", position = position_jitter(0.2, 0.2),
#              data = original_data_parameters, alpha = 0.8) + theme(legend.position="none")
# plot2
