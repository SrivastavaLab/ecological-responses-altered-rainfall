# cleaning up the environment ---------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# loading packages --------------------------------------------------------------------------------------------------------------------

library(tidyverse)

# loading data ------------------------------------------------------------------------------------------------------------------------

fulldata <- read.csv("Data/fulldata.csv")
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

sitemean <- filter(fulldata, mu.scalar == "1", k.scalar == "1") %>%
  select(site, mu.scalar, k.scalar, intended.mu, intended.k) %>% 
  rename(mu_mean = intended.mu, k_mean = intended.k)

sitecorners <- filter(fulldata, mu.scalar %in% c("0.1","3"),
                      k.scalar %in% c("0.5","2")) %>% 
  select(site, mu.scalar, k.scalar, intended.mu, intended.k)

sitecorners$order <- rep(c(1,2,4,3), 7)

sitecorners <- arrange(sitecorners, order)

polygon_data <- sitecorners %>% 
  select(site, intended.mu, intended.k, order) %>% 
  group_by(site) %>% 
  summarise(xmin = min(intended.mu),
            xmax = max(intended.mu),
            ymin = min(intended.k),
            ymax = max(intended.k)) %>% 
  left_join(y = sitemean, by = "site") %>% 
  select(-mu.scalar, -k.scalar)
polygon_data

polygon_data$site <- sub(pattern = "argentina", replacement = "Argentina", x = polygon_data$site)
polygon_data$site <- sub(pattern = "cardoso", replacement = "Cardoso", x = polygon_data$site)
polygon_data$site <- sub(pattern = "colombia", replacement = "Colombia", x = polygon_data$site)
polygon_data$site <- sub(pattern = "costarica", replacement = "Costa Rica", x = polygon_data$site)
polygon_data$site <- sub(pattern = "frenchguiana", replacement = "French Guiana", x = polygon_data$site)
polygon_data$site <- sub(pattern = "macae", replacement = "Macae", x = polygon_data$site)
polygon_data$site <- sub(pattern = "puertorico", replacement = "Puerto Rico", x = polygon_data$site)

# panel A -----------------------------------------------------------------------------------------------------------------------------

panelA <- ggplot() +
  geom_rect(data = polygon_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = site), fill = NA) +
  scale_colour_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  geom_point(data = fulldata, mapping = aes(x = intended.mu, y = intended.k, fill = site),
             colour = "grey80", shape = 21, size = 2, alpha = 0.5, na.rm = TRUE) + 
  geom_point(data = sitemean, mapping = aes(x = mu_mean, y = k_mean, fill = site),
              colour = "black", shape = 21, size = 4, alpha = 1, na.rm = TRUE) +
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = expression(bold(paste("Mean daily rainfall, ",mu, " (mm ",~day^-1,")"))), trans = "log2") +
  scale_y_continuous(name = "Rainfall dispersion, k", label = scaleFUN, trans = "log2") +
  coord_equal() +
  general_theme +
  theme(legend.position = "none")
panelA

# panel B -----------------------------------------------------------------------------------------------------------------------------

panelB <- ggplot(data = site_treatment_data, mapping = aes(x = intended.mu, y = intended.k, fill = site)) +
  facet_wrap(~ site) + 
  geom_point(colour = "grey80", shape = 21, size = 2, alpha = 0.5, na.rm = TRUE) + 
  geom_rect(data = polygon_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = site), fill = NA, inherit.aes = FALSE) +
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_x_continuous(name = expression(bold(paste("Mean daily rainfall, ",mu, " (mm)"))), trans = "log2") +
  scale_y_continuous(name = "Rainfall dispersion, k", label = scaleFUN, trans = "log2") +
  coord_equal() +
  geom_point(data = polygon_data, mapping = aes(x = mu_mean, y = k_mean, color = site), size = 4, alpha = 1, na.rm = TRUE) +
  scale_color_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold", colour = "black", size = 8),
        legend.position = "none")
panelB  

# mising data = original_data_parameters --
# Note from Andrew: as of Nov 2017 this "missing data" is created in 03_overall_figure_treatments.R

original_data_parameters <- read_csv("Data/original_precipitation_mu_k_parameters.csv") %>% 
  mutate(site = site_capital)

plot2<-panelB +
  geom_point(pch = 21, size = 2.2, fill = "white", position = position_jitter(0.2, 0.2),
             data = original_data_parameters, alpha = 0.8) + theme(legend.position="none")
plot2

# saving panels -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Figure_1_parameter_space.pdf", plot = panelA, width = 120, height = 89, dpi = 300, units = "mm")
ggsave(filename = "figure/figures/Figure_1_parameter_space_multipanel.pdf", plot = panelB, width = 120, height = 160, dpi = 300, units = "mm")

ggsave(filename = "figure/figures/Figure_1_parameter_space.png", plot = panelA, width = 120, height = 89, dpi = 600, units = "mm")
