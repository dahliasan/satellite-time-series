# satallite data climatology plots
library(tidyverse)
library(viridis)
# library(maps)
library(mapdata)
library(marmap)
library(lubridate)
library(ggpubr)
library(wesanderson)

rm(list = ls())
# load list of dataframes
load('./extracted enviro data/satellite data list (-36).RData')
l <- ll[c("ssha", "sst", "chl", "scu", "scv", "ssta")]

# calculate climatology (already saved output, don't need to run again) ---------------------------------------------------
clm <- lapply(l, FUN = function(x) {
  x %>% mutate(month = month(date, abbr = TRUE, label = TRUE)) %>%
    group_by(x, y, month) %>%
    summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
})

# log chl values before calculating anomaly
clm$chl <- clm$chl %>% mutate(mean = log(mean), sd = log(sd))

## Calculate climatology anomalies for some variables
calcClimAnom <- function(climdat, overalldat, log = FALSE){
  # climdat = monthly climatology dataframe colname = x, y, month, mean
  # overalldat = original time series dataframe which you calculated climatology from colnames = x, y, date, v1
  
  # overall time series mean 
  cm <- overalldat %>% 
    group_by(x, y) %>%
    summarise(mean = mean(v1, na.rm = T))
  
  if(log) cm <- mutate(cm, mean = log(mean))
    
  # calculate climatology anomalies (= monthly values - overall time series mean)
  left_join(climdat, cm, by = c('x', 'y')) %>% 
    mutate(mean = mean.x - mean.y) %>% 
    dplyr::select(x, y, month, mean)
}

## monthly chl anomaly climatology
clm$chlA <- calcClimAnom(clm$chl, ll$chl, log = T)

## monthly SST anomaly climatology
clm$ssta <- calcClimAnom(clm$sst, ll$sst)

## monthly surface current u and v anomaly climatology
clm$scua <- calcClimAnom(clm$scu, ll$scu)
clm$scva <- calcClimAnom(clm$scv, ll$scv)

# merge surface current u and v into a single data frame 
clm$sc <- left_join(clm$scu, clm$scv, by = c('x', 'y', 'month')) %>% 
  rename(um = mean.x, vm = mean.y, usd = sd.x, vsd = sd.y) %>% 
  mutate(mag = sqrt(um^2 + vm^2))
clm$sca <- left_join(clm$scua, clm$scva, by = c('x', 'y', 'month')) %>% 
  rename(um = mean.x, vm = mean.y) %>% 
  mutate(mag = sqrt(um^2 + vm^2))

# remove unwanted variables from list and save for easy loading in the future
clm <- clm[c("ssha", "sst", "chl", "chlA", "ssta", "sc", "sca")]
save(clm, file = 'climatology_ssha_sst_chl_chlA_ssta_sc_sca_dataframe.RData')

# calc upwelling wind index climatology
load("./extracted enviro data/bonneycoast_10m_upwelling_wind_index_1997-2017.RData" )
uw_clm <- uw %>% 
  mutate(month = as.Date(format(date, '1997-%m-01'))) %>% 
  group_by(month) %>% 
  summarise(upwl_wnd = mean(upwell_wind, na.rm = T), upwl_wnd_sd = sd(upwell_wind, na.rm = T))

load("~/satellite time series/extracted enviro data/climatology_ssha_sst_chl_chlA_ssta_sc_sca_dataframe.RData")
# Main plots (for publication) -------------------------------------------------------------------
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
mytheme <- theme(text = element_text(size = 8.5),
                 panel.background = element_rect(fill = "white", colour = 'black'),
                 # plot.title = element_text(hjust = 0.5, margin = margin(0,0,1,0,'pt')) # center title,
                 plot.title = element_text(margin = margin(0,0,1,0,'pt')))

## my colour palette
pal <- wes_palette('Zissou1')

## SSHA - looks good!
p1 <- ggplot(data = clm$ssha, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.2) +
  facet_wrap(~ month) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.5, fill = 'black') +
  scale_fill_gradient2(high = last(pal), low = first(pal), name = 'SSHA')  +
  geom_map(map_id = "Australia", map = map) +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(b)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5, size = 1) +
  mytheme

tiff(filename = './plots/SSHA climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p1)
dev.off()

## ACHL x currents - looks good!
scaler = 4
butext <- data.frame(x = 140.3, y = -37.4, month = as.factor("Feb"))

p2 <- ggplot(data = clm$chlA, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.5, fill = 'black') +
  geom_segment(data = clm$sca,
               arrow = arrow(length = unit(0.1, 'cm')),
               aes(x = x, y = y, xend = x + um * scaler, yend = y + vm * scaler),
               lwd = 0.2, colour = 'black', show.legend = FALSE) +
  facet_wrap(~ month) +
  scale_fill_gradient2(high = first(pal), low = last(pal), name = 'ChlA')  +
  geom_map(map_id = "Australia", map = map) +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(a)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5, size = 1) +
  mytheme

p2b <- p2 + geom_text(data = butext,  label = 'BC', colour = 'white', size =3)

tiff(filename = './plots/ChlA x Surface Currents climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p2b)
dev.off()


# SSTA --------------------------------------------------------------------

# l$ssta %>% 
#   filter(year(date) == 2016) %>% 
#   mutate(month = month(date)) %>% 
#   group_by(x,y, month) %>% 
#   summarise(mean = mean(v1, na.rm = T)) %>% 
clm$ssta %>% 
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  facet_wrap(~ month) +
  scale_fill_gradient2(high = last(pal), low = first(pal), name = 'SSTA')  +
  # scale_fill_viridis() + 
  geom_map(map_id = "Australia", map = map, colour = "grey") +
  lims(y = c(-43, -36), x = c(136, 141)) +
  labs(y = 'Lat', x = 'Lon') +
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme


# SST ---------------------------------------------------------------------
p3 <- ggplot(data = clm$sst, aes(x, y)) +
  geom_raster(aes(fill = mean)) +
  geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
  geom_contour(aes(z = sd), color = 'grey10', size = 0.08) +
  facet_wrap(~ month) +
  scale_fill_gradientn(colours = c(first(pal), 'white',last(pal)), name = 'SST')  +
  geom_map(map_id = "Australia", map = map) +
  ylim(-43, -36) +
  xlim(136, 141) +
  labs(y = 'Lat', x = 'Lon', title = '(c)') + 
  geom_path(data = sb, aes(x = x, y = y), color = 'black', alpha = 0.5) +
  mytheme
tiff(filename = 'SST climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p3)
dev.off()


# upwelling wind stress ---------------------------------------------------
# line and error bar climatology
uw_clm %>% 
  ggplot(aes(month, upwl_wnd, ymin = upwl_wnd - upwl_wnd_sd, ymax = upwl_wnd + upwl_wnd_sd)) + 
  geom_line() +
  geom_point() + 
  geom_errorbar(width = 3) + 
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  mytheme

# boxplot climatology
df <- uw %>% mutate(month = as.Date(format(date, '1997-%m-01'))) 
by_month <- df %>% group_by(month) %>% nest()
ylim <- purrr::map(by_month$data,  ~boxplot.stats(.$upwell_wind)$stats[c(1,5)]) %>% unlist()
ylim1 <- c(min(ylim), max(ylim))
p0 <- df %>% 
  ggplot(aes(month, upwell_wind, group = month))+
  # geom_boxplot(outlier.shape=NA) +
  # coord_cartesian(ylim = ylim1*1.05) + # removes outliers 
  geom_boxplot() +
  scale_x_date(date_breaks = '1 month', date_labels = '%b') +
  labs(y = expression(paste("Wind stress (", Nm^-2, ")", sep = ' ')), x = 'Month') + 
  mytheme 

tiff(filename = './plots/bonneyCoast_upwellingWindStress_climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p0)
dev.off()

# monthly mean sum of windstress 
u <- uw %>% mutate(month = month(date, abbr = T, label = T), year = year(date)) %>% 
  group_by(month, year) %>% 
  summarise(upwell_wind = sum(upwell_wind, na.rm = T))

p0 <- u %>% 
  ggplot(aes(month, upwell_wind, group = month))+
  geom_boxplot() +
  labs(y = expression(paste("Average cumulative wind stress (", Nm^-2, ")", sep = ' ')), x = 'Month') + 
  mytheme +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', size = 0.5)


tiff(filename = './plots/bonneyCoast_cumulativeMonthlyUpwellingWindStress_climatology.tiff',  width=7, height=7, units= "in", res = 300)
print(p0)
dev.off()

## combine SSHA, ACHL, SST
# ggarrange(p2,p1,p3, ncol = 1, nrow = 3, align = 'v')


# speciific year ------------------------------------------------
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
# calculate seasonal chl anomaly
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

# GLOB ---------------------------------------------------------------

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



# ascat curl --------------------------------------------------------------
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



# quikscat curl ------------------------------------------------------

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


# quikscat + ascat wind stress curl ----------------------------------------------------------
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


# wind direction --------------------------------------------------
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
 
# chl x currents ----------------------------------------------------------
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