source("01_datasets_for_paper1.R")
source("02_functions.R")

##====an overall figure of treatments

sitemean <- filter(fulldata,mu.scalar == "1", k.scalar == "1") %>%
  select(site, mu.scalar, k.scalar, intended.mu, intended.k)

sitecorners <- filter(fulldata,mu.scalar %in% c("0.1","3"),
                      k.scalar %in% c("0.5","2")) %>% 
  select(site, mu.scalar, k.scalar, intended.mu, intended.k)

sitecorners$order<-rep(c(1,2,4,3), 7)

sitecorners<-arrange(sitecorners, order)

# Read in the data that was originally used to estimate the distributions

library(purrr)

originals_data <- dir("Data/original_rainfall_data/", full.names = TRUE) %>% 
  set_names(., basename(.)) %>% 
  map(read_csv) %>% 
  map_df(~.x %>% gather(sampletime, ppt), .id = "sitename")



fit_distrs <- originals_data %>% 
  mutate(ppt = round(ppt)) %>% 
  group_by(sitename, sampletime) %>% nest %>% 
  mutate(data = map(data, unlist)) %>% 
  mutate(distr = map(data, safely(fitdist), distr = "nbinom", start = list(mu=9, size=0.7)))

# How many datasets are in each? 
fit_distrs %>% 
  dplyr::select(sitename, sampletime) %>% distinct %>% 
  group_by(sitename) %>% 
  dplyr::summarize(n = n())

# 5 to 26


nbin_fits <- fit_distrs %>% 
  mutate(distr_outpus = map(distr, 
                            ~ .x[["result"]] %>% coef %>% enframe)) %>% 
  unnest(distr_outpus)

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

site_names <- c(
  argentina="Argentina",
  cardoso="Cardoso",
  colombia="Colombia",
  costarica="Costa Rica",
  frenchguiana="French Guiana",
  macae="Macae",
  puertorico="Puerto Rico"
)

# add these capitalized names to the original_data_parameters

original_data_params_name <- original_data_parameters %>% 
  mutate(site_capital = site_names[site])

site_treatment_data <- fulldata %>% 
  dplyr::select(intended.mu, intended.k, site) %>% 
  mutate(site = ordered(site, levels = c("argentina", "macae", "frenchguiana",
                  "cardoso", "colombia", "costarica",  
                  "puertorico")))



# writing out some of this data!  -----------------------------------------

# writing out original data parameters -- mu and k as it was calculated before the schedule:

original_data_params_name %>% 
  write_csv("Data/original_precipitation_mu_k_parameters.csv")

write.csv(site_treatment_data, "Data/site_treatment_data.csv")
write.csv(fulldata, "Data/fulldata.csv")

scaleFUN <- function(x) sprintf("%.2f", x)

plot1<-ggplot(data=fulldata, aes(x=intended.mu, y=intended.k, fill=site))+
    geom_point(size=2, alpha=1/4, na.rm=TRUE, aes(color=site))+ 
    scale_x_continuous(trans='log2') +
    scale_y_continuous(label=scaleFUN,trans='log2')+
    geom_point(data=sitemean,size=4, alpha=1, na.rm=TRUE, aes(color=site))+ 
    geom_polygon(aes(color = site), data = sitecorners, fill = NA) +
    coord_equal()+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
    labs(x = expression(paste("Mean daily rainfall, ",mu, " (mm)")), y = "Rainfall dispersion, k")


png(file = "figure/ExperimentalSetup.png", units = "px", width = 7780, height=7380, res=1440)

plot1 

dev.off()


# Printing out.. with different colours

png(file = "figure/ExperimentalSetup_brewer.png", units = "px", width = 7780, height=7380, res=1440)

plot1 + 
  scale_color_brewer(type ="qual", palette = 2)

dev.off()



#Appendix plot

png(file = "figure/HistoricalRain.png", units = "px", width = 7780, height=7380, res=1440)

plot1b<-ggplot()+
  geom_point(data=site_treatment_data, aes(x=intended.mu, y=intended.k, fill=site), size=2, alpha=1/4, na.rm=TRUE, aes(color=site))+ 
  geom_polygon(mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = site), data = polygon_data, fill = NA) +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(label=scaleFUN,trans='log2')+
  geom_point(data=sitemean,size=4, alpha=1, na.rm=TRUE, aes(color=site))+
  coord_equal()+
  facet_wrap(~site, labeller=labeller(site=site_names)) + 
  #theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x = expression(paste("Mean daily rainfall, ",mu, " (mm)")), y = "Rainfall dispersion, k") + 
  theme_dark()

plot2<-plot1b + 
  geom_point(pch = 21, size = 2.2, fill = "white", position = position_jitter(0.2, 0.2),
             data = original_data_parameters, alpha = 0.8) + theme(legend.position="none")
plot2

dev.off()




