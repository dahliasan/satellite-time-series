# make dataframes from raster

library(dplyr)
library(tidyr)
library(raster)
library(stringr)

# load raster files
ssha <- stack("daily_ssha_136_141_-36_-43.grd")
chl <- stack("8day_chl_136_141_-36_-43.grd")
sst <- stack("daily_sst_136_141_-36_-43.grd")
glob <- stack('globoc_8day_136_141_-36_-43.grd')
ascat <- stack('ascatcurl_mday_136_141_-36_-43.grd')
qs <- stack('quikscatcurl_mday_136_141_-36_-43.grd')
scu <- stack("oscar_um_5day_136_141_-36_-43.grd")
scv <- stack("oscar_vm_5day_136_141_-36_-43.grd")
uwnd <- stack('uwnd_5-8day136_141_-36_-43.grd')
vwnd <- stack('vwnd_5-8day136_141_-36_-43.grd')
ssta <- brick("./extracted enviro data/ssta1day_129.5_145.5_-36_-46_1997-01-01_2017-12-31_.grd") 

# convert to raster to dataframe
vn <- c('ssha', 'sst', 'chl', 'glob','ascurl', 'qscurl', 'scu', 'scv', 'uwind', 'vwind', 'ssta')
l <- list(ssha, sst, chl, glob, ascat, qs, scu, scv, uwnd, vwnd, ssta)

# create list of data frames
ll <- lapply(l, FUN = function(x) {
  df <- tbl_df(data.frame(rasterToPoints(x)))
  t <- names(x)
  t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
  t <- as.character(gsub('\\.', '-', t))
  colnames(df)[3:length(df)] <- t
  df <- df %>% gather("date", 'v1', 3:length(df))
  df$date <- as.Date(df$date)
  return(df)
})

# assign correct variable name 
names(ll) <- vn

save(ll, file = './extracted enviro data/satellite data list (-36).Rdata')

## ekman upwelling data (2009 - present)
uw <- stack('upwelling_1day_2009-present.grd')
uw2 <- lapply(list(uw), FUN = function(x) {
    df <- tbl_df(data.frame(rasterToPoints(x)))
    t <- names(x)
    t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
    t <- as.character(gsub('\\.', '-', t))
    colnames(df)[3:length(df)] <- t
    df <- df %>% gather("date", 'v1', 3:length(df))
    df$date <- as.Date(df$date)
    return(df)
  })[[1]]
  uw <- uw2
  save(uw, file = 'upwelling (2009 - 2017).Rdata')

## ekman upwelling data (12 Sep 2017)
# dir(pattern = 'upwell')
# uw <- stack('upwelling_1day_136_141_-36_-43.grd')
# uw2 <- lapply(list(uw), FUN = function(x) {
#   df <- tbl_df(data.frame(rasterToPoints(x)))
#   t <- names(x)
#   t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
#   t <- as.character(gsub('\\.', '-', t))
#   colnames(df)[3:length(df)] <- t
#   df <- df %>% gather("date", 'v1', 3:length(df))
#   df$date <- as.Date(df$date)
#   return(df)
# })
# 
# uw <- uw2[[1]]
# save(uw, file = 'upwelling (1999 - 2017).Rdata')
