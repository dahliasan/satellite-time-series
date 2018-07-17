### bfast (trend breaks in time series) analysis zonal shelf vs oceanic waters ###
library(tidyverse)
library(ggplot2)
library(bfast)
library(lubridate)
library(marmap)


#load data 
rm(list = ls())
load('./extracted enviro data/satellite data list (-36).RData')
load('aao.Rdata')
load('./raw enviro data/soi_bom.RData') # monthly time series
# load('oni.Rdata')

# remove tbl_df attribute from dataframes 
aao <- as.data.frame(aao)
soi <- as.data.frame(soi)
ll <- lapply(ll, as.data.frame)


# remove leap days for sst and ssha (daily data)
ll$ssha <- ll$ssha[-with(ll$ssha, which(month(date) == 2 & day(date) == 29)),]
ll$sst <- ll$sst[-with(ll$sst, which(month(date) == 2 & day(date) == 29)),]
aao <- aao[-with(aao, which(month == 2 & day == 29)),]
soi$date <- as.Date(soi$date)
soi$month <- month(soi$date)

# log chl data
ll$chl <- ll$chl %>% mutate(v1 = log(v1))

# calculate chl anomaly 
ll$chlA <- ll$chl %>% 
  group_by(x, y) %>% 
  mutate(overall.mean = mean(v1, na.rm = T)) %>% 
  group_by(x, y, date) %>% 
  mutate(anom = v1 - overall.mean) %>% 
  dplyr::select(-overall.mean, -v1) %>% 
  rename(v1 = anom)


# select desired data type, Only interested in SSHA, SST, CHL
ll2 <- list(ll$ssha, ll$sst, ll$chl)
names(ll2) <- c('ssha', 'sst', 'chl')

# get bathy data (boundary for shelf)
bathy <- getNOAA.bathy(136,141,-43, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000)
sb <- b2000 %>% group_by(x) %>% summarise(y = min(y))

# Extract only shelf values
z <- lag(sb$x)
z[1] <- 136
shf <- lapply(ll2, function(df){
  tmp <- do.call(rbind, Map(function(x, y, z){
    df[df$x <= x & df$x > z & df$y >= y,]
  }, sb$x, sb$y, z))
  return(tmp)
})
# save(shf, file = 'shelf_ssha_sst_chl.Rdata')
load('shelf_ssha_sst_chl.Rdata')
shf <- lapply(shf, as.data.frame)

# choose x (lon) boundaries
shf2 <- lapply(shf, function(x, minlon = 138, maxlon = 140){
  x %>% filter(x >= minlon & x <= maxlon)
})
# calculate shelf regional mean 
shfmean <- lapply(shf2, function(x){
  x %>% group_by(date) %>% summarise(v1 = mean(v1, na.rm =TRUE))
})

# calculate oceanic region mean
oc <- lapply(ll2, function(x){
  x %>% filter(y <= -39 & y >= -42)
})
oc2 <- lapply(oc, function(x, minlon = 138, maxlon = 140){
  x %>% filter(x >= minlon & x <= maxlon)
})
ocmean <- lapply(oc2, function(x){
  x %>% group_by(date) %>% summarise(v1 = mean(v1, na.rm =TRUE))
})

# check how many obervations per year in each data set. Complete temporal periods are 1998 - 2016.
lapply(shfmean, function(x){
  x %>% group_by(year(date)) %>% summarise(n = n())})
lapply(ocmean, function(x){
  x %>% group_by(year(date)) %>% summarise(n = n())})

# start timeseries for sst, ssha and chl from 1998-01-01
shf98 <- lapply(shfmean, function(x){
  x %>% filter(year(date) > 1997 & year(date) < 2017)})
oc98 <- lapply(ocmean, function(x){
  x %>% filter(year(date) > 1997 & year(date) < 2017)})
aao98 <- aao %>% filter(year(date) > 1997 & year(date) < 2017)
soi98 <- soi %>% filter(year(date) > 1997 & year(date) < 2017)
# oni98 <- oni %>% filter(year(date) > 1997 & year(date) <2017)

# convert into 8 day for non-chl data
vars <- c('ssha', 'sst')
convert8day <- function(x) {
  x %>% group_by(year(date)) %>% 
    mutate(date = cut(date, 46)) %>% group_by(date) %>% 
    #variable to be summarise should be named 'v1'
    summarise(v1 = mean(v1, na.rm = TRUE)) %>% mutate(date = as.Date(date)) 
} 
shf98[vars] <- lapply(shf98[vars], ~convert8day(.))
oc98[vars] <- lapply(oc98[vars], ~convert8day(.))
aao98 <- convert8day(aao98)

# load upwelling alongshore wind stress data
load( "./extracted enviro data/bonneycoast_10m_upwelling_wind_index_1997-2017.RData")

# only get data from 1998 onwards, and convert to 8 day time series
uw98 <- uw %>% 
  filter(year(date) > 1997) %>% 
  rename(v1 = upwell_wind) %>% 
  convert8day()

#--- Run all code above before continuing ---#
# BFAST analysis starts here ----------------------------------------------
# make ts objects
shf98ts <- lapply(shf98, function(x){
  ts(x$v1, start = c(1998,1), frequency = 46)
})
oc98ts <- lapply(oc98, function(x){
  ts(x$v1, start = c(1998,1), frequency = 46)
})
aao98ts <- with(aao98, ts(v1, start = c(1998,1), frequency = 46))
soi98ts <- with(soi98, ts(v1, start = c(1998,1), frequency = 12))
uw98ts <- with(uw98, ts(v1, start = c(1998,1), frequency = 46))

# make bfast objects
# 9 May 2018:  the default h = 0.15, which I've changed to h = 0.21 which is = ~ 4 years minimum for a segment, also tried max breaks = 2
shf98bf <- lapply(shf98ts, function(x){
  bfast(x, season="harmonic", max.iter=1, h = 0.15)
})
oc98bf <- lapply(oc98ts, function(x){
  bfast(x,  season="harmonic", max.iter=1, h = 0.15)
})
aao98bf <- list(bfast(aao98ts, season="harmonic", max.iter=1, h = 0.15))
names(aao98bf) <- 'sam'
soi98bf <- bfast(soi98ts, season="harmonic", max.iter=1)
uw98bf <- bfast(uw98ts, season="harmonic", max.iter=1, h = 0.15)

# save all bfast objects
bfl <- list(shf98bf, oc98bf, aao98bf)
save(bfl, file = 'bfast zonal objects 6Sep2017 h0.26.Rdata')
save(soi98bf, file = './output/bfast/bfast soi 1998-2016.Rdata')
save(uw98bf, file = 'bfast_upwellWindStress_1998-2017.Rdata')
