# calculate SSTA based on Nieblas 2009 (ssta = daily - monthly mean)

library(tidyverse)
library(lubridate)
library(wesanderson)
library(mapdata)
library(marmap)
options(dplyr.width = Inf)

load('./extracted enviro data/satellite data list (-36).RData')
sst <- ll$sst

sst <- sst %>% mutate(year = year(date), month = month(date, label = T, abbr = T)) 

# calculate monthly mean
sst <- sst %>%
  group_by(year, month, x, y) %>% 
  mutate(monthly_mean = mean(v1, na.rm = T))
  
# calculate anomaly
sst <- sst %>% ungroup() %>% 
  mutate(ssta = v1 - monthly_mean)

# visualise
pal <- wes_palette('Zissou1')
map <- map_data("worldHires")
sst %>% ungroup() %>% 
  filter(year == 2016, month == 'Mar') %>% 
  ggplot(aes(x, y, fill = ssta)) + 
  geom_raster() + 
  facet_wrap(~date) + 
  geom_map(map_id = "Australia", map = map, fill = 'grey40') +
  scale_fill_gradient2(high = last(pal), low = first(pal), name = 'ssta')

# climatology anomaly
sst.cli <- sst %>% 
  group_by(month, x, y) %>% 
  summarise(sst = mean(v1, na.rm = T)) %>% 
  group_by(x,y) %>% 
  mutate(global_mean = mean(sst, na.rm = T), ssta = sst - global_mean)

# stf
stf <- sst.cli %>% filter(sst >= 12 & sst < 12+1)

# get bathy data
bathy <- getNOAA.bathy(136,141,-43, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000)
sb <- b2000 %>% group_by(x) %>% summarise(y = min(y))

## My theme
mytheme <- theme(text = element_text(size = 8.5),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 # plot.title = element_text(hjust = 0.5, margin = margin(0,0,1,0,'pt')) # center title,
                 plot.title = element_text(margin = margin(0,0,1,0,'pt')))

p1 <- sst.cli %>% ungroup() %>% 
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = ssta)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.5, fill = 'black') +
  geom_contour(aes(z = ssta), size = 0.2, colour = 'grey10') + 
  facet_wrap(~month) + 
  geom_map(map_id = "Australia", map = map) +
  scale_fill_gradient2(high = last(pal), low = first(pal), name = 'SSTA') + 
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(c)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5, size = 1) +
  mytheme

tiff(filename = './plots/SSTA climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p1)
dev.off()