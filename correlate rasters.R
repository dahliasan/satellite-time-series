### Correlate SSHA and CHL rasters ###
library(raster)
library(lubridate)
library(stringr)
library(spatialEco)
rm(list = ls())

ssha <- stack("daily_ssha_136_141_-36_-43.grd")
chl <- stack("8day_chl_136_141_-36_-43.grd")
myextent <- extent(136, 141, -43, -36)
s <- raster(myextent, nrows = nrow(ssha), ncols = ncol(ssha), crs = ssha@crs)
# chl2 <- aggregate(chl, fact = 6)
chl2 <- resample(chl, s, method="bilinear")
chl2 <- log(chl2)


## convert to monthly averages
# r <- ssha

ry <- lapply(list(ssha, chl2), function(r){
  t <- names(r)
  t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
  t <- as.character(gsub('\\.', '-', t))
  t <- as.Date(t)
  r <- r[[-which(year(t) < 1998 | year(t) > 2016)]]
  t <- names(r)
  t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
  t <- as.character(gsub('\\.', '-', t))
  t <- as.Date(t)
  xy <- stackApply(r, month(t), mean)
  names(xy) <- unique(month(t))
  xy <- setZ(xy, unique(month(t)))
  return(xy)
})

r.cor <- rasterCorrelation(ry[[1]][[1]], ry[[2]][[1]], s = 3, type = "pearson")
plot(r.cor)

rcorall <- ry[[1]]
for(i in 1:12){
  print(i)
r.cor <- corLocal(ry[[1]][[i]], ry[[2]][[i]], test=TRUE, method = 'spearman')
# plot(r.cor) # r and p-value map will be shown
# only correlation cells where the p-value < 0.05 will be ommited
rcor.mask <- mask(r.cor[[1]], r.cor[[2]] < 0.05, maskvalue=FALSE)
names(rcor.mask) <- i
# plot(rcor.mask)
rcorall[[i]] <- rcor.mask
}


## aggregate total time series
ry2 <- lapply(list(ssha, chl2), function(r){
  t <- names(r)
  t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
  t <- as.character(gsub('\\.', '-', t))
  t <- as.Date(t)
  r <- r[[-which(year(t) < 1998 | year(t) > 2016)]]
  t <- names(r)
  t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
  t <- as.character(gsub('\\.', '-', t))
  t <- as.Date(t)
  xy <- mean(r, na.rm = TRUE)
  return(xy)
})

r.cor <- corLocal(ry2[[1]], ry2[[2]], test=TRUE, method = 'pearson')
plot(r.cor) # r and p-value map will be shown
# only correlation cells where the p-value < 0.05 will be ommited
rcor.mask <- mask(r.cor[[1]], r.cor[[2]] < 0.05, maskvalue=FALSE)
plot(rcor.mask)
