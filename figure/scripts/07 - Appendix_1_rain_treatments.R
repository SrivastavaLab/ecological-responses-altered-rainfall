# cleaning up the environment ---------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

# loading packages --------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(fitdistrplus)

# loading data ------------------------------------------------------------------------------------------------------------------------

fulldata <- read.csv("Data/fulldata.csv")
site_treatment_data <- read.csv("Data/site_treatment_data.csv")

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

rename_sites <- function(x){
  x$site <- sub(pattern = "argentina", replacement = "Argentina", x = x$site)
  x$site <- sub(pattern = "cardoso", replacement = "Cardoso", x = x$site)
  x$site <- sub(pattern = "colombia", replacement = "Colombia", x = x$site)
  x$site <- sub(pattern = "costarica", replacement = "Costa Rica", x = x$site)
  x$site <- sub(pattern = "frenchguiana", replacement = "French Guiana", x = x$site)
  x$site <- sub(pattern = "macae", replacement = "Macae", x = x$site)
  x$site <- sub(pattern = "puertorico", replacement = "Puerto Rico", x = x$site)
  return(x)
}

# tidying data for graphic creation ---------------------------------------------------------------------------------------------------

sitemean <- filter(fulldata, mu.scalar == "1", k.scalar == "1") %>%
  dplyr::select(site, intended.mu, intended.k) %>% 
  rename(mu_mean = intended.mu, k_mean = intended.k)

polygon_data <- filter(fulldata, mu.scalar %in% c("0.1","3"),
                      k.scalar %in% c("0.5","2")) %>% 
  dplyr::select(site, intended.mu, intended.k) %>% 
  group_by(site) %>% 
  summarise(xmin = min(intended.mu),
            xmax = max(intended.mu),
            ymin = min(intended.k),
            ymax = max(intended.k)) %>% 
  left_join(y = sitemean, by = "site")
polygon_data

polygon_data <- rename_sites(polygon_data)
fulldata <- rename_sites(fulldata)

# original data used to create rainfall schedulles ------------------------------------------------------------------------------------

originals_data <- dir("Data/original_rainfall_data/", full.names = TRUE) %>% 
  set_names(., basename(.)) %>% 
  map(read_csv) %>% 
  map_df(~.x %>% gather(key = sampletime, value = ppt), .id = "sitename") %>% 
  mutate(ppt = round(ppt)) %>% 
  group_by(sitename, sampletime) %>% 
  nest %>% 
  mutate(data = map(data, unlist)) %>% 
  mutate(distr = map(data, safely(fitdist), distr = "nbinom", start = list(mu=9, size=0.7)))
originals_data

nbin_fits <- originals_data %>% 
  mutate(distr_outpus = map(distr, ~ .x[["result"]] %>% coef %>% enframe)) %>% 
  unnest(distr_outpus)
nbin_fits

correct_sitenames <- tribble(
  ~sitename,                     ~site,
  "Ppt_Argentina.csv",            "argentina",
  "Ppt_Pitilla_Oct9Dec7.csv",     "costarica",
  "Ppt_PuertoRico.csv",           "puertorico",
  "Ppt.Colombiaprecipitation.csv","colombia",
  "Ppt.French.Guiana.csv",        "frenchguiana",
  "PptCardoso.csv",               "cardoso",
  "PptMacae.csv",                 "macae"
)

original_data_parameters <- nbin_fits %>% 
  left_join(y = correct_sitenames, by = "sitename") %>% 
  dplyr::select(site, name, value, sampletime) %>% 
  spread(name, value) %>% 
  rename(intended.mu = mu, intended.k = size)
original_data_parameters

original_data_parameters <- rename_sites(original_data_parameters)
original_data_parameters

# creating figure ---------------------------------------------------------------------------------------------------------------------

app1 <- ggplot(data = site_treatment_data) +
  geom_point(data = site_treatment_data, aes(x = intended.mu, y = intended.k, fill = site), shape = 21, size = 2, alpha = 0.5) + 
  geom_rect(data = polygon_data, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = site), fill = NA) +
  scale_x_continuous(name = expression(bold(paste("Mean daily rainfall, ",mu, " (mm)"))), trans = "log2") +
  scale_y_continuous(name = "Rainfall dispersion, k", label = scaleFUN, trans = "log2") +
  geom_point(data = polygon_data, size = 4, alpha = 1, na.rm = TRUE, aes(x = mu_mean, y = k_mean, fill = site), shape = 21, size = 4) +
  coord_equal() +
  facet_wrap(~ site) + 
  geom_point(data = original_data_parameters, mapping = aes(x = intended.mu, y = intended.k), shape = 21, size = 3, 
             fill = "grey70", colour = "black", stroke = 0.6, alpha = 0.6) + 
  scale_fill_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  scale_colour_manual(values = c("mediumpurple4", "darkorange2", "forestgreen", "deepskyblue1", "mediumblue", "lawngreen", "firebrick1")) +
  general_theme +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 8, color = "black"),
        strip.background = element_blank())
app1

# saving panels -----------------------------------------------------------------------------------------------------------------------

ggsave(filename = "figure/figures/Appendix_1_parameter_space.pdf", plot = app1, width = 160, height = 160, dpi = 300, units = "mm")

ggsave(filename = "figure/figures/Appendix_1_parameter_space.png", plot = app1, width = 160, height = 160, dpi = 600, units = "mm")



