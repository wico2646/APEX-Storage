library(ggplot2)
library(viridis)
library(dplyr)

attach(Alpha_flux_combined)

ggplot(data = Alpha_flux_combined, aes(x = Carex.RAbundance, y = CH4.mgC.m2.d, color = WT.Treatment)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() +
  ggtitle("Effect of water table and Carex abundance on methane flux") +
  xlab("Carex relative abundance") +
  ylab(expression("Methane flux (mgC • m"^-2 * " • day"^-1*")"))

ggplot(data = Alpha_flux_combined, aes(x = Veg.Treatment, y = NEE.mgC.m2.d, color = WT.category)) +
  geom_boxplot() +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() +
  ggtitle("Effect of water table and plant funtional group removal on NEE") +
  xlab("Vegetation removal treatment") +
  ylab(expression("NEE (mgC • m"^-2 * " • day"^-1*")"))

data_2019 <- Alpha_flux_combined %>% filter(Year == 2019)
data_2021 <- Alpha_flux_combined %>% filter(Year == 2021)
data_subset <- Alpha_flux_combined %>% filter(Year %in% c(2016, 2019, 2021))

ggplot(data = data_subset, aes(x = Carex.RAbundance, y = CH4.mgC.m2.d, color = WT.Treatment)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_color_viridis(discrete=TRUE) +
  theme_bw(base_size = 20) +
  ggtitle("Effect of water table and Carex abundance on methane flux") +
  xlab("Carex relative abundance") +
  ylab(expression("Methane flux (mgC • m"^-2 * " • day"^-1*")"))

ggplot(data = data_subset, aes(x = Veg.Treatment, y = NEE.mgC.m2.d, color = WT.Treatment)) +
  geom_boxplot() +
  scale_color_viridis(discrete=TRUE) +
  theme_bw(base_size = 20) +
  ggtitle("Effect of water table and plant funtional group removal on NEE") +
  xlab("Vegetation removal treatment") +
  ylab(expression("NEE (mgC • m"^-2 * " • day"^-1*")"))


ggplot(data = Alpha_flux_combined, aes(x = Julian.Day, y = Water.level, color = factor(Year))) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() 
  