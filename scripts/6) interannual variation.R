## Calculate interannual variability ##
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(marmap)
library(lubridate)
library(ggpubr)

## Load list of raw data
load('./extracted enviro data/satellite data list (-36).RData')

## Extract desired variables
ll <- ll[c('ssha', 'sst', 'chl')]

## load season classification function
getSeason <- function(DATES) {
  wi <- as.Date("2012-06-01", format = "%Y-%m-%d") # Winter
  sp <- as.Date("2012-09-01",  format = "%Y-%m-%d") # Spring
  su <- as.Date("2012-12-01",  format = "%Y-%m-%d") # Summer
  au <- as.Date("2012-03-01",  format = "%Y-%m-%d") # Fall

  ## use dummy year since we're only interested at seasonal level
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))

  ifelse (d >= wi & d < sp, "Winter",
          ifelse (d >= sp & d < su, "Spring",
                  ifelse (d >= au & d < wi, "Autumn", "Summer")))
}

## average time series values by season
s <- lapply(ll, function(x){
  x$date[which(month(x$date) == 12)] <- x$date[which(month(x$date) == 12)] + 365 # add 365 days to december month as it's part on the next year's summer period. 
  x %>% mutate(year = year(date), season = getSeason(date)) %>%
    filter(year > 1997 & year < 2017) %>%
    group_by(x, y, season, year) %>%
    summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
  })

## save(s, file = 'satellite data seasonaly (-36).Rdata')

## Calculate yearly values
yearly <- lapply(ll, FUN = function(x) {
  x %>% mutate(year = year(date)) %>%
    group_by(x, y, year) %>%
    filter(year > 1997 & year < 2017) %>%
    summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
})

# save(yearly, file = 'satellite data yearly (-36).Rdata')

load('./extracted enviro data/satellite data yearly (-36).Rdata')
load('./extracted enviro data/satellite data seasonaly (-36).Rdata')

# log chl values
yearly$chl <- yearly$chl %>% mutate(mean = log(mean), sd = log(sd))
s$chl <- s$chl %>% mutate(mean = log(mean), sd = log(sd))


# Calculate mean for entire time series for each cell  --------------------
## Yearly
tsm <- lapply(yearly, function(df) df %>% group_by(x, y) %>% summarise(mean = mean(mean, na.rm = T)))

yy <- Map(function(u, df){
  df$tsmean <- NA
  r <- rep(1:nrow(u), each = n)
  df$tsmean <- u$mean[r]
  df$iav <- df$mean - df$tsmean
  return(df)
}, u = tsm[1:3], df = yearly[1:3])

## Seasonal
tsm <- lapply(s, function(df) df %>% group_by(x, y, season) %>% summarise(mean = mean(mean, na.rm = T)))

# Check that each cell-season combination has the same number of observations (n) 
# lapply(s, function(df){
#   a <- df %>% group_by(x, y, season) %>% summarise(n = n())})
n <- unique(unlist(lapply(s, function(df){
  a <- df %>% group_by(x, y, season) %>% summarise(n = n())
  unique(a$n)
  })))
n

# Run below if above check is correct i.e. same number of observations for each group
## Calculate annual anomalies
# u <- tsm$scu
# df <- s$scu
ss <- Map(function(u, df){
  df$tsmean <- NA
  # r <- unlist(Map(function(x, y, z) {
  #   which(u$x == x & u$y == y & u$season == z)}, df$x, df$y, df$season))
  r <- rep(1:nrow(u), each = n)
  df$tsmean <- u$mean[r]
  df$iav <- df$mean - df$tsmean
  df$season <- factor(df$season, levels = c('Summer', 'Autumn', 'Winter', 'Spring'))
  return(df)
}, u = tsm[1:3], df = s[1:3])

# Make plots --------------------------------------------------------------
# get map 
map <- map_data("worldHires")
# Get shelf break boundary
bathy <- getNOAA.bathy(136,141,-43, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000)
sb <- b2000 %>% group_by(x) %>% summarise(y = min(y))

## ggplot chl x surface currents
# to use this function below need to run above to get objects map, sb, stf, smin, smax
chlxScplot <- function(pd, sc, stf, scaler = 3, p.title = NULL, fn = 'plot1'){
  
  p <- ggplot(data = pd, aes(x, y)) +
    geom_raster(aes(fill = iav)) +
    facet_wrap(~ year) +
    geom_map(map_id = "Australia", map = map) +
    geom_segment(data = sc, 
                 arrow = arrow(length = unit(0.08, 'cm')), 
                 aes(x = x, y = y, xend = x + uiav * scaler, yend = y + viav * scaler), 
                 lwd = 0.1, colour = 'black', size = 1, show.legend = FALSE) +
    ylim(-43, -36) +
    xlim(136, 141) +
    geom_path(data = sb, aes(x = x, y = y), color = 'black') +
    geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
    theme(text = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, margin = margin(0,0,1,0,'pt'))) + 
    labs(title = p.title, x = 'Lon', y = 'Lat', fill = 'CHLAA') +
    scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue')  
    # scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', limits = c(smin, smax))  
  
  tiff(filename = paste0(fn, '.tiff'),  width=7, height=7, units= "in", res = 300)
  print(p)
  dev.off()
  
  return(p)
  
}
IAVplot <- function(pd, stf, p.title, fn, legend.name){
  p <- ggplot(data = pd, aes(x, y)) +
    geom_raster(aes(fill = iav)) +
    facet_wrap(~ year) +
    geom_map(map_id = "Australia", map = map, colour = "grey") +
    ylim(-43, -36) +
    xlim(136, 141) +
    geom_path(data = sb, aes(x = x, y = y), color = 'black') +
    geom_raster(data = stf, aes(x = x, y = y), alpha = 0.1, fill = 'black') +
    labs(title = p.title) +
    theme(text = element_text(size = 8.5),
          plot.title = element_text(margin = margin(0,0,1,0,'pt'))) + 
    labs(title = p.title, x = 'Lon', y = 'Lat', fill = legend.name) +
    scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue')
    # scale_fill_gradient2(high = 'firebrick2', low = 'midnightblue', limits = c(smin, smax))
  
  tiff(filename = paste0(fn, '.tiff'),  width=7, height=7, units= "in", res = 300)
  print(p)
  dev.off()
  
  return(p)
}

# Seasonal plots ----------------------------------------------------------
# get stf data
stf <- s$sst %>% filter(mean >= 12 & mean < 12+1)
stf <- split(stf, stf$season)
stf <- list(stf$Summer, stf$Autumn, stf$Winter, stf$Spring)
# split chl and sc by season
chlss <- split(ss$chl, f = ss$chl$season)
scss <- split(ss$sc, f = ss$sc$season)
# set scale limits
smin <- min(ss$chl$iav, na.rm = TRUE)
smax <- max(ss$chl$iav, na.rm = TRUE)
pt <- paste(names(chlss), 'CHL x Surface Currents')
pt <- list(pt[1], pt[2], pt[3], pt[4])
chlsc_plots <- Map(chlxScplot, pd = chlss, stf = stf, sc = scss, p.title = pt, fn = pt)

## ggplot - ssha, sst
# split chl and sst by season
sshass <- split(ss$ssha, f = ss$chl$season)
sstss <- split(ss$sst, f = ss$sc$season)
# plot ssha
# set scale limits
smin <- min(ss$ssha$iav, na.rm = TRUE)
smax <- max(ss$ssha$iav, na.rm = TRUE)
pt <- paste(names(sshass), 'SSHA')
pt <- list(pt[1], pt[2], pt[3], pt[4])
ssha_plots <- Map(IAVplot, pd = sshass, stf = stf, p.title = pt, fn = pt, legend.name = 'SSHAA')
# plot sst
smin <- min(ss$sst$iav, na.rm = TRUE)
smax <- max(ss$sst$iav, na.rm = TRUE)
pt <- paste(names(sshass), 'SST')
pt <- list(pt[1], pt[2], pt[3], pt[4])
sst_plots <- Map(IAVplot, pd = sstss, stf = stf, p.title = pt, fn = pt, legend.name = 'SSTAA')
  
## Calculate SD of annual anomaly (used for final manuscript)
ggplotSDIAV <- function(df, title, varname = '', pal = NULL, rm.out = FALSE){
  pal <- wesanderson::wes_palette('Zissou1')[-c(2,4)]
  df <- df %>% group_by(x, y, season) %>% summarise(meanIAV = mean(iav), sdIAV = sd(iav))
  if(rm.out == TRUE){
    hs <- hist(df$sdIAV, plot = FALSE)
    df <- df %>% filter(sdIAV <= hs$breaks[last(which(hs$counts == 1))])
  }
  mytheme <- theme(text = element_text(size = 8.5),
                   panel.background = element_rect(fill = "white", colour = 'black'),
                   plot.title = element_text( margin = margin(0,0,1,0,'pt')),
                   legend.position = "bottom",
                   legend.margin = margin(-1,0,0,0,'mm'))
  p <- ggplot(data = df, aes(x, y)) +
    geom_raster(aes(fill = sdIAV)) +
    facet_wrap(~ season) +
    geom_map(map_id = "Australia", map = map) +
    ylim(-43, -36) +
    xlim(136, 141) +
    geom_path(data = sb, aes(x = x, y = y), color = 'black') +
    labs(title = title, x = 'Lon', y = 'Lat', fill = paste0(varname)) +
    scale_fill_gradientn(colours = pal, name = 'iav') +
    mytheme
  return(p)
}
p1 <- ggplotSDIAV(ss$chl, '(a) Chl-a', pal = pal)
p2 <- ggplotSDIAV(ss$ssha, '(b) SSHA', pal = pal)
p3 <- ggplotSDIAV(ss$sst, '(c) SST', pal = pal)
gg <- ggarrange(p1,p2,p3, ncol = 3, nrow = 1, align = 'h')
tiff(filename = './plots/chl ssha sst sd iav.tiff',  width=6, height=3.5, units= "in", res = 300)
print(gg)
dev.off()


## Get range of annual anomalies
range <- lapply(ss[1:3], function(df) {
  a <- df %>% group_by(season) %>% summarise(min_AA = min(iav, na.rm = T), max_AA = max(iav, na.rm = T))
  a <- as.data.frame(a)
  return(a)
})
range$um <- lapply(ss[4], function(df) {
  a <- df %>% group_by(season) %>% summarise(min_AA = min(uiav, na.rm = T), max_AA = max(uiav, na.rm = T))
  a <- as.data.frame(a)
  return(a)
})[[1]]
range$vm <- lapply(ss[4], function(df) {
  a <- df %>% group_by(season) %>% summarise(min_AA = min(viav, na.rm = T), max_AA = max(viav, na.rm = T))
  a <- as.data.frame(a)
  return(a)
})[[1]]

range <- do.call(rbind, range)
range$min_AA <- signif(range$min_AA, 3)
range$max_AA <- signif(range$max_AA, 3)
range$range <- range$max_AA - range$min_AA
write.csv(range, file = 'annual anomalies range.csv')



# Annual plots ------------------------------------------------------------
# get stf data
stf <- yearly$sst %>% filter(mean >= 12 & mean < 12+1)
chlxScplot(pd = yy$chl, sc = yy$sc, stf = stf, scaler = 3, p.title = '(a)', fn = 'CHL x Surface Currents Annual Anomalies.tiff')
IAVplot(pd = yy$ssha, stf = stf, p.title = '(b)', fn = 'SSHA Annual Anomalies.tiff', legend.name = 'SSHAA' )
IAVplot(pd = yy$sst, stf = stf, p.title = '(c)', fn = 'SST Annual Anomalies.tiff', legend.name = 'SSTAA' )       

ggplotSDIAV <- function(df, title, varname = '', pal = 'Spectral'){
  df <- df %>% group_by(x, y) %>% summarise(meanIAV = mean(iav), sdIAV = sd(iav))
  mytheme <- theme(text = element_text(size = 6),
                   panel.background = element_rect(fill = "white", colour = 'black'),
                   plot.title = element_text(hjust = 0.5, margin = margin(0,0,1,0,'pt')),
                   legend.position = "bottom",
                   legend.margin = margin(-1,0,0,0,'mm'))
  p <- ggplot(data = df, aes(x, y)) +
    geom_raster(aes(fill = sdIAV)) +
    geom_map(map_id = "Australia", map = map, colour = "grey") +
    ylim(-43, -36) +
    xlim(136, 141) +
    geom_path(data = sb, aes(x = x, y = y), color = 'black') +
    labs(title = title, x = 'Lon', y = 'Lat', fill = paste0(varname)) +
    scale_fill_distiller(type = 'div', palette = pal) +
    mytheme
  return(p)
}
p1 <- ggplotSDIAV(yy$chl, '(a) CHL iav')
p2 <- ggplotSDIAV(yy$ssha, '(b) SSHA iav', pal = 'RdYlBu')
p3 <- ggplotSDIAV(yy$sst, '(c) SST iav', pal = 'RdYlBu')
gg <- ggarrange(p1,p2,p3, ncol = 3, nrow = 1, align = 'h')
tiff(filename = 'annual chl ssha sst sd iav.tiff',  width=6, height=3.5, units= "in", res = 300)
print(gg)
dev.off()

range <- lapply(yy[1:3], function(df) {
  a <- df %>% ungroup() %>% summarise(min_AA = min(iav, na.rm = T), max_AA = max(iav, na.rm = T))
  a <- as.data.frame(a)
  return(a)
})
range$um <- lapply(yy[4], function(df) {
  a <- df %>% ungroup() %>% summarise(min_AA = min(uiav, na.rm = T), max_AA = max(uiav, na.rm = T))
  a <- as.data.frame(a)
  return(a)
})[[1]]
range$vm <- lapply(yy[4], function(df) {
  a <- df %>% ungroup () %>% summarise(min_AA = min(viav, na.rm = T), max_AA = max(viav, na.rm = T))
  a <- as.data.frame(a)
  return(a)
})[[1]]

range <- do.call(rbind, range)
range$min_AA <- signif(range$min_AA, 3)
range$max_AA <- signif(range$max_AA, 3)
range$range <- range$max_AA - range$min_AA
write.csv(range, file = 'annual anomalies range2.csv')



# inter-annual variation (annual annomaly) for wind stress ----------------------------------
# special case for a single location time series

# load upwelling wind stress separately 
load( "./extracted enviro data/bonneycoast_10m_upwelling_wind_index_1997-2017.RData")

# seasonal time series
x <- uw %>% rename(v1 = upwell_wind)
x$date[which(month(x$date) == 12)] <- x$date[which(month(x$date) == 12)] + 365 # add 365 days to december month as it's part on the next year's summer period. 
uwsy <- x %>% mutate(year = year(date), season = getSeason(date)) %>%
  filter(year > 1997 & year < 2017) %>%
  group_by(season, year) %>%
  summarise(mean = mean(v1, na.rm = TRUE))

# overall mean of time series
tsmean <- mean(uw$upwell_wind)

# seasonal inter-annual anomaly
uwsy <- uwsy %>% group_by(season, year) %>% mutate(anom = mean - tsmean)
uw_iav <- uwsy %>% 
  group_by(season) %>%
  summarise(sd = sd(anom, na.rm = T), min = min(anom), max = max(anom))
  
