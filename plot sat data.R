# satallite data climatology plots
library(dplyr)
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(marmap)
library(lubridate)
library(ggpubr)

rm(list = ls())
# load list of dataframes
load('satellite data list (-36).Rdata')
# l <- ll[c(1:3,6:8)]
# # Calculate climatology ---------------------------------------------------
# clm <- lapply(l, FUN = function(x) {
#   x %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>% 
#     group_by(x, y, month) %>% 
#     summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
# })
# 
# # log chl values
# clm$chl <- clm$chl %>% mutate(mean = log(mean), sd = log(sd))
# 
# # merge scv and scu df
# clm$scu$vm <- clm$scv$mean
# colnames(clm$scu)[4] <- 'um'
# # clm$scu <- clm$scu %>%  mutate(mag = sqrt(um^2 + vm^2))
# clm <- clm[-6]
# save(clm, file = 'climatology ssha sst chl qscurl sc dataframe.Rdata')
load('climatology ssha sst chl qscurl sc dataframe.Rdata')

## Calculate chl and sst and surface current anomalies
## CHL anomaly
# calculate climatology mean
# cm <- clm$chl %>% group_by(x, y) %>% summarise(mean = mean(mean, na.rm = T))
cm <- ll$chl %>% group_by(x, y) %>% summarise(mean = mean(v1, na.rm = T)) %>% mutate(mean = log(mean))
# calculmate seasonal chl anomaly
a <- clm$chl %>% group_by(x, y) %>% summarise(n = n())
n <- unique(a$n)
n
clm$achl <- Map(function(u, df){
  df <- df %>% select(-sd)
  df$cm <- NA
  r <- rep(1:nrow(cm), each = n)
  df$cm <- u$mean[r]
  df$mean <- df$mean - df$cm
  return(df)
}, u = list(cm), df = list(clm$chl))[[1]]

# is climatology mean the same as the mean of all observations for each cell? -
# slightly different may be rounding error
# tmp <- ll$chl %>% group_by(x,y) %>% mutate(v1 = log(v1)) %>% summarise(mean = mean(v1, na.rm = T))
# tmp2 <- Map(function(u, df){
#   df <- df %>% select(-sd)
#   df$cm <- NA
#   r <- rep(1:nrow(cm), each = n)
#   df$cm <- u$mean[r]
#   df$mean <- df$mean - df$cm
#   return(df)
# }, u = list(tmp), df = list(clm$chl))[[1]]

## SST anomaly
cm <- clm$sst %>% group_by(x, y) %>% summarise(mean = mean(mean, na.rm = T))
# calculmate seasonal chl anomaly
a <- clm$sst %>% group_by(x, y) %>% summarise(n = n())
n <- unique(a$n)
n
clm$ssta <- Map(function(u, df){
  df <- df %>% select(-sd)
  df$cm <- NA
  r <- rep(1:nrow(cm), each = n)
  df$cm <- u$mean[r]
  df$mean <- df$mean - df$cm
  return(df)
}, u = list(cm), df = list(clm$sst))[[1]]

## surface current anomaly
cm <- clm$scu %>% group_by(x, y) %>% summarise(umean = mean(um, na.rm = T), vmean = mean(vm, na.rm = T))
ll$scu$v2 <- ll$scv$v1
cm <- ll$scu %>% group_by(x,y) %>% summarise(umean = mean(v1, na.rm = T), vmean = mean(v2, na.rm = T))
a <- clm$scu %>% group_by(x, y) %>% summarise(n = n())
n <- unique(a$n)
n
clm$sca <- Map(function(u, df){
  df$umean <- NA
  df$vmean <- NA
  r <- rep(1:nrow(u), each = n)
  df$umean <- u$umean[r]
  df$vmean <- u$vmean[r]
  df$uma <- df$um - df$umean
  df$vma <- df$vm - df$vmean
  df <- df %>% select(-c(um,sd,vm,umean,vmean))
  return(df)
}, u = list(cm), df = list(clm$scu))[[1]]

## Upwelling
load('upwelling (2009 - 2017).Rdata')

uwclm <- lapply(list(uw), FUN = function(x) {
  x %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>%
    
    group_by(x, y, month) %>%
    filter(year(date) < 2017) %>% 
    summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
})[[1]]


#### make plots ####
# get map 
map <- map_data("worldHires")

# stf
stf <- clm$sst %>% filter(mean >= 12 & mean < 12+1)

# get bathy data
bathy <- getNOAA.bathy(136,141,-43, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000)
sb <- b2000 %>% group_by(x) %>% summarise(y = min(y))

## My theme
mytheme <- theme(text = element_text(size = 8),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 plot.title = element_text(hjust = 0.5, margin = margin(0,0,1,0,'pt')))

## Plot SSHA - looks good!
p1 <- ggplot(data = clm$ssha, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.08) +
  facet_wrap(~ month) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'SSHA')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(b)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme

tiff(filename = 'SSHA climatology2.tiff',  width=7, height=7, units= "in", res = 300)
print(p1)
dev.off()

## plot ACHL x surface currents - looks good!
scaler = 4
butext <- data.frame(x = 140.3, y = -37.4, month = as.factor("Feb"))

p2 <- ggplot(data = clm$achl, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  geom_segment(data = clm$sca, 
               arrow = arrow(length = unit(0.1, 'cm')), 
               aes(x = x, y = y, xend = x + uma * scaler, yend = y + vma * scaler), 
               lwd = 0.2, colour = 'black', show.legend = FALSE) +
  facet_wrap(~ month) +
  scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'ACHL')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(a)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme

p2b <- p2 + geom_text(data = butext,  label = 'BC', colour = 'white', size =3)

tiff(filename = 'CHLA x Surface Currents climatology3.tiff',  width=7, height=7, units= "in", res = 300)
print(p2b)
dev.off()

## plot SSTA - meh
# ggplot(data = clm$ssta, aes(x, y)) +
#   geom_raster(aes(fill = ssta)) +
#   facet_wrap(~ month) +
#   scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'SSTA')  +
#   geom_map(map_id = "Australia", map = map, colour = "grey") +
#   ylim(-43, -36) +
#   labs(y = 'Lat', x = 'Lon') + 
#   geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
#   mytheme

## plot SST
p3 <- ggplot(data = clm$sst, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.08) +
  facet_wrap(~ month) +
  # scale_fill_distiller(name = 'SST', type = 'div', palette = 'RdBu') +
  scale_fill_gradientn(colours = c('midnightblue', 'white','firebrick2'), name = 'SST')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(c)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme
tiff(filename = 'SST climatology2.tiff',  width=7, height=7, units= "in", res = 300)
print(p3)
dev.off()

## plot upwelling
uwclm2 <- uwclm
uwclm2$mean[uwclm2$mean < 0] <- NA

p4 <- ggplot(data = uwclm, aes(x, y)) +
  geom_raster(aes(fill = mean > 0)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.08) +
  facet_wrap(~ month) +
  # scale_fill_viridis() +
  scale_fill_brewer(type = 'qual', palette = 2, name = 'Upwelling') +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(d)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme

tiff(filename = 'Upwelling climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p4)
dev.off()

## combine SSHA, ACHL, SST
# ggarrange(p2,p1,p3, ncol = 1, nrow = 3, align = 'v')


# plot speciific year ------------------------------------------------
l <- ll[c(1:3,6:8)]
yr <- 2007
clm16 <- lapply(l, FUN = function(x) {
  x %>% filter(year(date) == yr) %>% 
    mutate(month = month(date, abbr = TRUE, label = TRUE)) %>%
    group_by(x, y, month) %>%
    summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
})
# log chl values
clm16$chl <- clm16$chl %>% mutate(mean = log(mean), sd = log(sd))
# merge scv and scu df
clm16$scu$vm <- clm16$scv$mean
colnames(clm16$scu)[4] <- 'um'
clm16 <- clm16[-6]
# calculmate seasonal chl anomaly
cm <- ll$chl %>% filter(year(date) == yr) %>% group_by(x, y) %>% summarise(mean = mean(v1, na.rm = T)) %>% mutate(mean = log(mean))
a <- clm16$chl %>% group_by(x, y) %>% summarise(n = n())
n <- unique(a$n)
n
clm16$achl <- Map(function(u, df){
  df <- df %>% select(-sd)
  df$cm <- NA
  r <- rep(1:nrow(cm), each = n)
  df$cm <- u$mean[r]
  df$mean <- df$mean - df$cm
  return(df)
}, u = list(cm), df = list(clm16$chl))[[1]]
## surface current anomaly
ll$scu$v2 <- ll$scv$v1
cm <- ll$scu %>% filter(year(date) == yr) %>% group_by(x,y) %>% summarise(umean = mean(v1, na.rm = T), vmean = mean(v2, na.rm = T))
a <- clm16$scu %>% group_by(x, y) %>% summarise(n = n())
n <- unique(a$n)
n
clm16$sca <- Map(function(u, df){
  df$umean <- NA
  df$vmean <- NA
  r <- rep(1:nrow(u), each = n)
  df$umean <- u$umean[r]
  df$vmean <- u$vmean[r]
  df$uma <- df$um - df$umean
  df$vma <- df$vm - df$vmean
  df <- df %>% select(-c(um,sd,vm,umean,vmean))
  return(df)
}, u = list(cm), df = list(clm16$scu))[[1]]
stf <- clm16$sst %>% filter(mean >= 12 & mean < 12+1)
scaler = 4
tiff(filename = 'achlxsc2016monthly.tiff',  width=7, height=7, units= "in", res = 300)
ggplot(data = clm16$achl, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  geom_segment(data = clm16$sca,
               arrow = arrow(length = unit(0.1, 'cm')),
               aes(x = x, y = y, xend = x + uma * scaler, yend = y + vma * scaler),
               lwd = 0.2, colour = 'black', show.legend = FALSE) +
  # geom_segment(data = clm16$scu, 
  #              arrow = arrow(length = unit(0.1, 'cm')), 
  #              aes(x = x, y = y, xend = x + um * scaler, yend = y + vm * scaler), 
  #              lwd = 0.2, colour = 'black', show.legend = FALSE) +
  facet_wrap(~ month) +
  scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'ACHL')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = yr) + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme
dev.off()

tiff(filename = 'ssha2016monthly.tiff',  width=7, height=7, units= "in", res = 300)
ggplot(data = clm16$ssha, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.08) +
  facet_wrap(~ month) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'SSHA')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = yr) + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme
dev.off()

# plot sst ssha chl with for loop ----------------------------------------------
## plot for sst, ssha, chl
# pn <- c('p1', 'p2', 'p3')
# vn <- c('ssha', 'sst', 'chl')
# for(i in 1:length(pn)){
# p <- ggplot(data = clm[[i]], aes(x, y)) +
#   geom_raster(aes(fill = mean)) +
#   geom_contour(aes(z = sd), color = 'white', size = 0.1) +
#   facet_wrap(~ month) +
#   scale_fill_viridis() +
#   geom_map(map_id = "Australia", map = map, colour = "grey") +
#   ylim(-43, -36) +
#   labs(title = vn[i]) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'black', alpha = 0.5, weight = 2)
# tiff(filename = paste(vn[i], '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
# print(p)
# dev.off()
# assign(pn[i],p)
# }

# #### plot surface current ####
# scaler = 5
# p <- ggplot(data = clm$scu, 
#             aes(x = x, y = y)) + 
#   geom_raster(aes(fill = mag)) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.3, weight = 2)+
#   geom_segment(arrow = arrow(length = unit(0.1, 'cm')), aes(xend = x + um * scaler, yend = y + vm * scaler), lwd = 0.5) +
#   facet_wrap(~ month) + 
#   scale_fill_viridis() + 
#   geom_map(map_id = "Australia", map = map, colour = "grey") + 
#   labs(title = 'surface current')
#   
# tiff(filename = paste('surface current', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
# p
# dev.off()
# 
# 

# plot GLOB ---------------------------------------------------------------

tmp <- ll$glob
tmp <- tmp  %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>% 
  group_by(x, y, month) %>% 
  summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))

tmp <- tmp %>% mutate(mean_old = mean) %>% mutate(mean = log(mean), sd = log(sd))

(p <- ggplot(data = tmp, aes(x, y)) + 
    geom_raster(aes(fill = mean)) +
    geom_contour(aes(z = sd), color = 'white', size = 0.1) +
    facet_wrap(~ month) +
    scale_fill_viridis() +
    geom_map(map_id = "Australia", map = map, colour = "grey") +
    ylim(-43, -36) +
    labs(title = 'chl (glob oc)') +
    geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'black', alpha = 0.5, weight = 2))

tiff(filename = paste('chl glob', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
print(p)
dev.off()


# #### plot ascat curl ####
tmp <- ll$ascurl
tmp <- tmp  %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>% 
  group_by(x, y, month) %>% 
  summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))

p2 <- ggplot(tmp, aes(x = x, y = y)) +
  geom_tile(aes(fill = mean )) +
  geom_contour(aes(z = sd), color = 'white', size = 0.2) +
  geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.5, weight = 2) +
  facet_wrap(~ month) +
  ylim(-43, -36) +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  labs(title = 'ASCAT curl (2010 - 2016') 

mi <- (tmp %>%  filter(mean < 0))$mean
mibreaks <- c( median(mi),min(mi))
ma <- (tmp %>%  filter(mean > 0))$mean
mabreaks <- c( max(ma),median(ma))
my_colpal <- diverging_pal(11)[c(4:6,10:11)]
brlab <- as.character(signif(c(mabreaks, 0, mibreaks)*10^5, 3))
scname <- 'curl 10^-5'
(p2 <- p2 + scale_fill_gradientn(colours = my_colpal, 
                         values = rescale(c(mabreaks, 0, mibreaks)),
                         guide = "legend",
                         breaks = c(mabreaks, 0, mibreaks),
                         labels = brlab,
                         name = scname) )
  
tiff(filename = paste('ascat curl (2010 - 2016)', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
p2
dev.off()



# plot quikscat curl ------------------------------------------------------

p4 <- ggplot(data = clm$qscurl, aes(x, y)) +
  geom_raster(aes(fill = mean < 0)) +
  facet_wrap(~ month) +
  scale_fill_discrete(name = 'Curl < 0') + 
  # scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', name = 'Curl')  +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = 'QUIKSCAT Wind Stress Curl (1999 - 2009)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme

tiff(filename = paste('quikscatcurl (1999 - 2009) climatology2', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
p4
dev.off()


# tmp <- ll$qscurl
# tmp <- tmp  %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>% 
#   group_by(x, y, month) %>% 
#   summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
# 
# p2 <- ggplot(tmp, aes(x = x, y = y)) +
#   geom_tile(aes(fill = mean )) +
#   geom_contour(aes(z = sd), color = 'white', size = 0.2) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.5, weight = 2) +
#   facet_wrap(~ month) +
#   ylim(-43, -36) +
#   geom_map(map_id = "Australia", map = map, colour = "grey") +
#   labs(title = 'QUIKSCAT curl (1999 - 2009)') 
# 
# mi <- (tmp %>%  filter(mean < 0))$mean
# mibreaks <- c( median(mi),min(mi))
# ma <- (tmp %>%  filter(mean > 0))$mean
# mabreaks <- c( max(ma),median(ma))
# my_colpal <- diverging_pal(11)[c(4:6,10:11)]
# brlab <- as.character(signif(c(mabreaks, 0, mibreaks)*10^5, 3))
# scname <- 'curl 10^-5'
# (p2 <- p2 + scale_fill_gradientn(colours = my_colpal, 
#                                  values = rescale(c(mabreaks, 0, mibreaks)),
#                                  guide = "legend",
#                                  breaks = c(mabreaks, 0, mibreaks),
#                                  labels = brlab,
#                                  name = scname) )
# 
# tiff(filename = paste('quikcat curl (1999 - 2009)', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
# p2
# dev.off()

# plot my wind stress curl ----------------------------------------------------------
tmp <- ll$wcurl
tmp <- tmp  %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>% 
  group_by(x, y, month) %>% 
  summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))


p <- ggplot(data = tmp, aes(x, y)) + 
  geom_tile(aes(fill = mean > 0)) +
  geom_contour(aes(z = sd), color = 'white', size = 0.1) +
  facet_wrap(~ month) + 
  geom_map(map_id = "Australia", map = map, colour = "grey") + 
  ylim(-43, -36) + 
  labs(title = 'wind stress curl') + 
  geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'black', alpha = 0.5, weight = 2)
# p + scale_fill_brewer(palette = 'Blues', direction = -1)

tiff(filename = paste('mywindstresscurl', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
p
dev.off()


# plot wind vector field --------------------------------------------------
wind <- ll$uwind
colnames(wind)[4] <- 'u'
wind$v <- ll$vwind$v1
rho_air = 1.2
cd = 1.2e-3

wclm <- wind  %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>% 
  mutate(wmag = sqrt(u^2+v^2), u = rho_air*cd*wmag*u, v = rho_air*cd*wmag*v) %>% 
  group_by(x, y, month) %>% 
  summarise(um = mean(u, na.rm = TRUE), vm = mean(v, na.rm = TRUE)) 


scaler = 8
p <- ggplot(data = wclm, aes(x = x, y = y)) +
  geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.3, weight = 2)+
  geom_segment(arrow = arrow(length = unit(0.1, 'cm')), aes(xend = x + um * scaler, yend = y + vm * scaler), lwd = 0.5) +
  facet_wrap(~ month) +
  scale_fill_viridis() +
  xlim(136, 141) + 
  ylim(-43, -36) +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  labs(title = 'wind field')


# #### combine wind stress curl and surface current plots ####
# scaler = 6
# p <- ggplot(data = clm$curl) + 
#   geom_tile(data = clm$curl, aes(x = x, y = y, fill = mean < 0)) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.7) +
#   geom_segment(data = clm$scu, 
#                arrow = arrow(length = unit(0.1, 'cm')), 
#                aes(x = x, y = y, alpha = mag, xend = x + um * scaler, yend = y + vm * scaler), 
#                lwd = 0.5) +
#   facet_wrap(~ month) + 
#   ylim(-43, -36) + 
#   geom_map(map_id = "Australia", map = map, colour = "grey") + 
#   labs(title = 'wind stress curl and surface current') 
# 
# tiff(filename = paste('wind stress curl and surface current', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
# p
# dev.off()
# 
#### chl x surface current plots ####
scaler = 5
p <- ggplot(data = clm$chl) +
  geom_tile(data = clm$chl, aes(x = x, y = y, fill = mean)) +
  geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.5) +
  geom_segment(data = clm$scu, 
               arrow = arrow(length = unit(0.1, 'cm')), 
               aes(x = x, y = y, xend = x + um * scaler, yend = y + vm * scaler), 
               lwd = 0.2, colour = 'black', size = 1) +
  facet_wrap(~ month) +
  ylim(-43, -36) +
  xlim(136, 141) +
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  labs(title = 'chl x surface current') +
  scale_fill_viridis()

tiff(filename = paste('chl and surface current', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
p
dev.off()
# 
# #### combine chl and wind stress curl ####
# p <- ggplot(data = clm$chl) + 
#   geom_point(data = clm$curl, aes(x = x, y = y, alpha = mean < 0)) +
#   geom_tile(data = clm$chl, aes(x = x, y = y, fill = mean), alpha = 0.9) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.7) +
#   facet_wrap(~ month) + 
#   ylim(-43, -36) + 
#   geom_map(map_id = "Australia", map = map, colour = "grey") + 
#   labs(title = 'wind stress curl and surface current') + 
#   scale_fill_viridis()
# 
# 
# # sst sd vs chl -----------------------------------------------
# 
# p <- ggplot(data = clm$chl) + 
#   geom_tile(data = clm$chl, aes(x = x, y = y, fill = mean), alpha = 0.9) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.7) +
#   geom_contour(data = clm$sst, aes(x = x, y =y, z = sd), color = 'white', size = 0.1) +
#   facet_wrap(~ month) + 
#   ylim(-43, -36) + 
#   geom_map(map_id = "Australia", map = map, colour = "grey") + 
#   labs(title = 'wind stress curl and surface current') + 
#   scale_fill_viridis()
# 
# 
# # sst sd vs surface current -----------------------------------------------
# 
# p <- ggplot(data = clm$sst) + 
#   geom_tile(data = clm$sst, aes(x = x, y = y, fill = sd))+
#   # geom_contour(data = clm$sst, aes(x = x, y = y, z = sd)) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.7) +
#   geom_segment(data = clm$scu, 
#                arrow = arrow(length = unit(0.1, 'cm')), 
#                aes(x = x, y = y, xend = x + um * scaler, yend = y + vm * scaler), 
#                lwd = 0.5) +
#   facet_wrap(~ month) + 
#   ylim(-43, -36) + 
#   geom_map(map_id = "Australia", map = map, colour = "grey") + 
#   labs(title = 'sst and surface current') + 
#   scale_fill_viridis()
# 
# tiff(filename = paste('sstsd and surface current', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
# p
# dev.off()
# 
# # plot hovmollerish surface current plots
# tmp <- clm$scu %>% group_by(x) %>% mutate(region = if(x < 137) {
#   region = 1} else if(x >= 137 & x < 138){
#     region = 2} else if(x >= 138 & x < 139){
#       region = 3} else {region = 4})
# tmp <- tmp %>% group_by(month, y, region) %>% mutate(um = mean(um, na.rm = TRUE),
#                                          vm = mean(vm, na.rm = TRUE))
# tmp <- tmp %>%  mutate(mag = sqrt(um^2 + vm^2))
# tmp2 <- tmp %>% group_by(month, y, region) %>% summarise(um = first(um),
#                                                          vm = first(vm), 
#                                                          mag = first(mag),
#                                                          x = first(x))
# 
# scaler = 5
# p <- ggplot(data = tmp, 
#             aes(x = x, y = y)) + 
#   geom_raster(aes(fill = mag)) +
#   geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.3, weight = 2)+
#   geom_segment(arrow = arrow(length = unit(0.1, 'cm')), aes(xend = x + um * scaler, yend = y + vm * scaler), lwd = 0.5) +
#   facet_wrap(~ month) + 
#   scale_fill_viridis() + 
#   geom_map(map_id = "Australia", map = map, colour = "grey") + 
#   labs(title = 'surface current')+
#   geom_vline(aes(xintercept = 137), linetype = 'dashed')+
#   geom_vline(aes(xintercept = 138), linetype = 'dashed')+
#   geom_vline(aes(xintercept = 139), linetype = 'dashed')
# tiff(filename = paste('surface current subregion split', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
# p
# dev.off()
