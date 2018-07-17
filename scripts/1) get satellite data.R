# Get satellite data 

# load libraries 
library(raster)
library(raadtools)
library(zoo)
library(rasterVis)
# library(ncdf4)
# library(xtractomatic)
library(viridis)
# library(rerddap)
library(anytime)
library(grid)
library(gridExtra)
library(animation)
library(remote) 
library(tidyverse)
library(lubridate)

rm(list = ls())
# set the spatial boundaries 
# myextent <- extent(136, 141, -43, -36)
myextent <- extent(129.5, 145.5, -46, -36) # for GLS kud analysis
start <- "1997-01-01"
end <- "2018-01-01"
extentfn <- paste(myextent@xmin, myextent@xmax, myextent@ymax, myextent@ymin, sep = '_')


# ssha --------------------------------------------------------------------
sshaf <- sshfiles(ssha = TRUE)
mysshafiles <- sshaf %>% filter(date >= as.POSIXct(start), date <= as.POSIXct(end)) %>% as.data.frame()

# crop to your extent
ssha <- crop(stack(mysshafiles$fullname, varname = "sla"), myextent)

# set names and time scale
tt <- as.Date(mysshafiles$date)
ssha <- setZ(ssha, tt)
names(ssha) <- as.character(tt)
writeRaster(ssha, filename = paste("daily_ssha", extentfn, start, end, '.grd', sep = '_'), overwrite=TRUE)

# sst ---------------------------------------------------------------------

sstf <- sstfiles()
mysstfiles <- sstf %>% filter(date >= as.POSIXct(start), date <= as.POSIXct(end))
mysstfiles <- as.data.frame(mysstfiles)

# range(mysstfiles$date)
# [1] "1997-01-01 11:00:00 AEDT" "2016-12-30 11:00:00 AEDT"

# get daily sst data for each month
sst <- crop(stack(mysstfiles$fullname, varname = "sst"), myextent)

# set names and time scale
tt <- as.Date(mysstfiles$date)
sst <- setZ(sst, tt)
names(sst) <- as.character(tt)
writeRaster(sst, filename = paste("daily_sst", extentfn, start, end, '.grd', sep = '_'), overwrite=TRUE)

# chl (one souce) ---------------------------------------------------------
fn <- paste('modischl8day', extentfn, '2003-01-05', '2018-01-01', '.tsv', sep = '_')
# download.file("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla8day.nc?chlorophyll[(2003-01-05):1:(2018-01-01T00:00:00Z)][(-36):1:(-46)][(129.5):1:(145.5)]",
              # destfile = fn)
download.file('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla8day.tsv?chlorophyll[(2003-01-05T00:00:00Z):1:(2018-01-01)][(80.02083):1:(80.02083)][(140.0208):1:(140.0208)]',
              destfile = fn)

mdr <- brick("modischl8day_129.5_145.5_-36_-46_2003-01-05_2018-01-01_.nc")

# Convert unix time to dates 
t <- data.frame(read_tsv("modischl8day_129.5_145.5_-36_-46_2003-01-05_2018-01-01_.tsv"))[-1,] # round a bout method to get time because .nc file doesn't have timestamps for some reason
t$time <- as.Date(strptime(t$time, format = '%Y-%m-%dT%H:%M:%SZ'))
t <- t$time
names(mdr) <- t
mdr <- setZ(mdr, t)

t <- seq(as.Date('2003-01-05'), by = '8 day', length.out = nlayers(mdr))

# resample raster to match sst and ssha (0.25)
r1 <- resample(mdr, sst, method="bilinear")

# join seawif and modis rasters
chl <- r1

fn <- paste('modischl8day', extentfn, '2003-01-05', '2018-01-01', '.grd', sep = '_')
writeRaster(chl, filename = fn, overwrite=TRUE)


# chl (two sources) ---------------------------------------------------------------------
# Combine data from two different sources e.g. SeaWIFs + MODIS 
## Xtractomatic method
# xpos <- c(136, 141)
# ypos <- c(-36.5, -43)
# tpos <- c('2003-01-01', '2016-12-31') 
## search for data 
# list1 <- list('varname', 'chl')
# list2 <- list('datasetname', '1day')
# mylist <- list(list1, list2)
# searchData(mylist)
# chlsw <- xtracto_3D(xpos, ypos, c('1997-09-04', '2002-12-31'), 'erdSW1chla1day', verbose = TRUE)
# chlmd <- xtracto_3D(xpos, ypos, tpos, 'erdMH1chla1day', verbose = TRUE)

# Direct download method (preferred)
# download.file("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdSW1chla8day.nc?chlorophyll[(1997-09-02):1:(2003-01-01T00:00:00Z)][(-36):1:(-43)][(136):1:(141)]",
#               destfile = 'seawifschl8day.nc')
download.file("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla8day.nc?chlorophyll[(2003-01-05):1:(2018-01-01T00:00:00Z)][(-36):1:(-46)][(129.5):1:(145.5)]",
              destfile = paste('modischl8day', extentfn, '2003-01-05', '2018-01-01', '.nc', sep = '_'))

swr <- stack('seawifschl8day (-36).nc')
mdr <- stack('modischl8day (-36).nc')

# get time names 
options("scipen" = 100, "digits" = 4)
t <- as.numeric(sapply(names(swr), function(x) str_split(x, 'X')[[1]][2]))
t <- anydate(t)
names(swr) <- t
swr <- setZ(swr, t)

# get time names
t <- as.numeric(sapply(names(mdr), function(x) str_split(x, 'X')[[1]][2]))
t <- anydate(t)
names(mdr) <- t
mdr <- setZ(mdr, t)
# mdr <- mdr[[-1]] # because sw and md date overlap 

# resample raster because sw and md not same extent
s <- raster(myextent, nrows = nrow(swr), ncols = ncol(swr), crs = swr@crs)
r1 <- resample(swr, s, method="bilinear")
r2 <- resample(mdr, s, method="bilinear")

# join seawif and modis rasters
chl <- stack(r1, r2)

writeRaster(chl, filename = paste("8day_chl_", extentfn, '.grd', sep = ''), overwrite=TRUE)

# make spatial resolution smaller
# chl_lowres <- aggregate(chlr, fact=3)
# chl_lowres <- setZ(chl_lowres, tt)
# writeRaster(chl_lowres, filename = "weekly_chl_137-139_-36.25--43 (low res).grd", overwrite=TRUE)


# glob colour -------------------------------------------------------------
library(RCurl)
url <- "ftp://ftp_hermes:hermes%@ftp.hermes.acri.fr/943062448/"
filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
gcftp <- read.csv('gcftp.csv')
colnames(gcftp) <- 'file'
gcftp$file <- as.character(gcftp$file)
filenames<- str_split(filenames[1], "\\n")[[1]]
filenames <- filenames[-length(filenames)]
# for(i in 1:length(filenames)){
#   download.file(paste(url, filenames[i], sep = ''), destfile = filenames[i])
# }

flist <- dir(pattern = 'L3m')
t <- lapply(str_split(flist, '_'), function(x) str_split(x[2], '-')[[1]][1])
t <- lapply(t, function(x) as.Date(x, format = '%Y%m%d'))
t <- do.call('c', t)
# save(flist, file = 'glob_flist.Rdata')
glob <- stack(flist)
names(glob)<- t
writeRaster(glob, filename=paste("globoc_8day_", extentfn, '.grd', sep = ''), overwrite=TRUE)

# oscar sea current -------------------------------------------------------

windu <- brick('/user/xhdfoo/.cache/R/rerddap/4f98d16b382f89c924553710c0738b15.nc',
               varname = 'um' )
t <- as.numeric(sapply(names(windu), function(x) str_split(x, 'X')[[1]][2]))
t <- anydate(t)
names(windu) <- t
windu <- setZ(windu, t)
writeRaster(windu, filename=paste("oscar_um_5day_", extentfn, '.grd', sep = ''), overwrite=TRUE)

windv <- brick('/user/xhdfoo/.cache/R/rerddap/4f98d16b382f89c924553710c0738b15.nc',
               varname = 'vm' )
t <- as.numeric(sapply(names(windv), function(x) str_split(x, 'X')[[1]][2]))
t <- anydate(t)
names(windv) <- t
windv <- setZ(windv, t)
writeRaster(windv, filename=paste("oscar_vm_5day_", extentfn, '.grd', sep = ''), overwrite=TRUE)

windsc <- sqrt(windu^2 + windv^2)

# get time names 
tt <- names(windsc)
tt <- sapply(tt, function(x) str_split(x, 'X')[[1]][2])
tt <- gsub('\\.', '-', tt)
tt <- as.Date(tt)


# wind stress curl --------------------------------------------------------
options("scipen"=100, "digits"=4)

# QUIKSCAT 1997 - 2009
download.file("coastwatch.pfeg.noaa.gov/erddap/griddap/erdQSstress1day_LonPM180.nc?curl[(1999-07-21):1:(2009-11-19T00:00:00Z)][(0.0):1:(0.0)][(-46):1:(-34)][(129):1:(146)]",
              destfile = 'quikscatcurl_1day.nc')
# download.file('http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQSstressmday_LonPM180.nc?curl[(1999-08-16T12:00:00Z):1:(2009-10-16T12:00:00Z)][(0.0):1:(0.0)][(-43):1:(-36)][(136):1:(143)]',
#               destfile = 'quikscatcurl_mday.nc')

curl1 <- brick('quikscatcurl_mday.nc')
tt <- names(curl1)
tt2 <- as.numeric(sapply(tt, function(x) str_split(x, 'X')[[1]][2]))
tt2 <- anydate(tt2)
names(curl1) <- tt2
writeRaster(curl1, filename=paste("quikscatcurl_mday_", extentfn, '.grd', sep = ''), overwrite=TRUE)

# ASCAT curl 2010 - 2016
download.file("coastwatch.pfeg.noaa.gov/erddap/griddap/erdQMstress8day_LonPM180.nc?curl[(2013-08-30):1:(2017-12-31T00:00:00Z)][(0.0):1:(0.0)][(-46):1:(-34)][(129):1:(146)]",
              destfile = '~/satellite time series/raw enviro data/ascatcurl_8day_129_146_-34_-46_20130830-20171231.nc')
# download.file('http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAstressmday.nc?curl[(2009-10-16T12:00:00Z):1:(2016-12-31T12:00:00Z)][(0.0):1:(0.0)][(-43):1:(-36)][(136):1:(143)]',
              # destfile = 'asscatcurl_mday.nc')
curl2 <- brick('ascatcurl_mday.nc')
tt <- names(curl2)
tt2 <- as.numeric(sapply(tt, function(x) str_split(x, 'X')[[1]][2]))
tt2 <- anydate(tt2)
names(curl2) <- tt2
writeRaster(curl2, filename=paste("ascatcurl_mday_", extentfn, '.grd', sep = ''), overwrite=TRUE)

# # resample raster because sw and md not same extent
# s <- raster(myextent, nrows = nrow(curl2), ncols = ncol(curl2), crs = curl2@crs)
# r1 <- resample(curl1lowres, s, method="bilinear")
# r2 <- resample(curl2, s, method="bilinear")
# 
# # join rasters
# curl <- stack(r1, r2)
# writeRaster(curl, filename=paste("windstresscurl_1day_", extentfn, '.grd', sep = ''), overwrite=TRUE)

# calculate climatology
# tt <- names(curl)
# tt <- sapply(tt, function(x) str_split(x, 'X')[[1]][2])
# tt <- gsub('\\.', '-', tt)
# tt <- as.Date(tt)
# curlclim <- stackApply(curl, month(tt), mean)
# tt2 <- unique(month(tt, abbr = TRUE, label = TRUE))
# names(curlclim) <- tt2
# writeRaster(curlclim, filename=paste("climatology_windstresscurl_1day_", extentfn, '.grd', sep = ''), overwrite=TRUE)


# upwelling ------------------------------------------------------------------

options("scipen"=100, "digits"=4)

# download 2009 - present ekman upwelling
download.file('coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAstress1day_LonPM180.nc?upwelling[(2009-10-07):1:(2017-09-06)][(0.0):1:(0.0)][(-43):1:(-36.5)][(136):1:(141)]',
              destfile = 'upwelling_8day_2009-present.nc')
uw <- stack('upwelling_1day_2009-present.nc')
tt <- names(uw)
tt2 <- as.numeric(sapply(tt, function(x) str_split(x, 'X')[[1]][2]))
tt2 <- anydate(tt2)
names(uw) <- tt2
writeRaster(uw, filename="upwelling_1day_2009-present.grd", overwrite=TRUE)

# download.file("coastwatch.pfeg.noaa.gov/erddap/griddap/erdQSstress1day_LonPM180.nc?upwelling[(1999-07-21):1:(2009-11-19T00:00:00Z)][(0.0):1:(0.0)][(-43):1:(-36.5)][(136):1:(141)]",
              # destfile = 'upwelling_1day(1).nc')
upwell1 <- stack('upwelling_1day(1).nc')
tt <- names(upwell1)
tt2 <- as.numeric(sapply(tt, function(x) str_split(x, 'X')[[1]][2]))
tt2 <- anydate(tt2)
names(upwell1) <- tt2
upwell1lowres <- aggregate(upwell1, fact=2)


# download.file("coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAstress1day_LonPM180.nc?upwelling[(2009-11-20):1:(2017-05-10T00:00:00Z)][(0.0):1:(0.0)][(-43):1:(-36.5)][(136):1:(141)]",
              # destfile = 'upwelling_1day(2).nc')
upwell2 <- stack('upwelling_1day(2).nc')
tt <- names(upwell2)
tt2 <- as.numeric(sapply(tt, function(x) str_split(x, 'X')[[1]][2]))
tt2 <- anydate(tt2)
names(upwell2) <- tt2

# resample raster because sw and md not same extent
s <- raster(myextent, nrows = nrow(upwell2), ncols = ncol(upwell2), crs = upwell2@crs)
r1 <- resample(upwell1lowres, s, method="bilinear")
r2 <- resample(upwell2, s, method="bilinear")

# join rasters
upwell <- stack(r1, r2)
writeRaster(upwell, filename=paste("upwelling_1day_", extentfn, '.grd', sep = ''), overwrite=TRUE)


# CCMP + ASCAT wind speeds ------------------------------------------------
options("scipen"=100, "digits"=4)

download.file("coastwatch.pfeg.noaa.gov/erddap/griddap/jplCcmp35aWindPentad_LonPM180.nc?uwnd[(1997-01-01T00:00:00Z):1:(2011-12-27T00:00:00Z)][(-43):1:(-36)][(136):1:(141)]",
              destfile = 'ccmp_uwnd.nc')
download.file("coastwatch.pfeg.noaa.gov/erddap/griddap/jplCcmp35aWindPentad_LonPM180.nc?vwnd[(1997-01-01T00:00:00Z):1:(2011-12-27T00:00:00Z)][(-43):1:(-36)][(136):1:(141)]",
              destfile = 'ccmp_vwnd.nc')
download.file("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAwind8day_LonPM180.nc?x_wind[(2011-01-01T00:00:00Z):1:(2016-12-26T00:00:00Z)][(10.0):1:(10.0)][(-43):1:(-36)][(136):1:(141)]",
              destfile = 'ascat_uwnd.nc')
download.file("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAwind8day_LonPM180.nc?y_wind[(2011-01-01T00:00:00Z):1:(2016-12-26T00:00:00Z)][(10.0):1:(10.0)][(-43):1:(-36)][(136):1:(141)]",
              destfile = 'ascat_vwnd.nc')

u1 <- stack('ccmp_uwnd.nc')
u2 <- stack('ascat_uwnd.nc')
v1 <- stack('ccmp_vwnd.nc')
v2 <- stack('ascat_vwnd.nc')

# join u wind rasters 
t <- names(u1)
tt <- as.numeric(sapply(t, function(x) str_split(x, 'X')[[1]][2]))
tt <- anydate(tt)
names(u1) <- tt

t <- names(u2)
tt <- as.numeric(sapply(t, function(x) str_split(x, 'X')[[1]][2]))
tt <- anydate(tt)
names(u2) <- tt

# resample rasters to join because not same extent
s <- raster(myextent, nrows = nrow(u1), ncols = ncol(u1), crs = u1@crs)
r1 <- resample(u1, s, method="bilinear")
r2 <- resample(u2, s, method="bilinear")

# make into single raster
uwnd <- stack(r1, r2)
writeRaster(uwnd, filename=paste("uwind_5-8day", extentfn, '.grd', sep = ''), overwrite=TRUE)


# join v wind rasters 
t <- names(v1)
tt <- as.numeric(sapply(t, function(x) str_split(x, 'X')[[1]][2]))
tt <- anydate(tt)
names(v1) <- tt

t <- names(v2)
tt <- as.numeric(sapply(t, function(x) str_split(x, 'X')[[1]][2]))
tt <- anydate(tt)
names(v2) <- tt

# resample rasters to join because not same extent
s <- raster(myextent, nrows = nrow(v1), ncols = ncol(v1), crs = v1@crs)
r1 <- resample(v1, s, method="bilinear")
r2 <- resample(v2, s, method="bilinear")

# make into single raster
vwnd <- stack(r1, r2)
writeRaster(vwnd, filename=paste("vwind_5-8day", extentfn, '.grd', sep = ''), overwrite=TRUE)

# calculate curl
rho_air = 1.2
cd = 1.2e-3
wmag = (uwnd^2+vwnd^2)^.5
f <- raster(myextent, nrows = nrow(taux), ncols = ncol(taux), crs = taux@crs)
yval <- yFromCell(taux, 1:ncell(taux))
pi = 3.142
omega = 7.292e-5
values(f) <- 2*omega*sin(yval*pi/180) # coriolis parameter
pw <- 1024 # density of water

#wind stress
taux = rho_air*cd*wmag*uwnd
tauy = rho_air*cd*wmag*vwnd

#wind curl
multiFocal <- function(x, w = matrix(1, nr=3, nc=3), ...) { 
  
  if(is.character(x)) { 
    x <- brick(x) 
  } 
  # The function to be applied to each individual layer 
  fun <- function(ind, x, w, ...){ 
    focal(x[[ind]], w=w, ...) 
  } 
  
  n <- seq(nlayers(x)) 
  list <- lapply(X=n, FUN=fun, x=x, w=w, ...) 
  
  out <- stack(list) 
  return(out) 
} 
dy <- matrix(c(1,1,1), nrow = 3, ncol = 1)
dx <- matrix(c(1,1,1), nrow = 1, ncol = 3)
dTxdy <- multiFocal(taux, w = dy, function(x) diff(rev(x), lag = 2)/(2*0.2414))
dTydx <- multiFocal(tauy, w = dx, function(x) diff(x, lag = 2)/(2*0.2381))
curl_tau = dTydx - dTxdy
names(curl_tau) <- names(taux)
wcurl <- curl_tau/-pw
names(wcurl) <- names(taux)
writeRaster(curl_tau, filename=paste("wcurl_5-8day", extentfn, '.grd', sep = ''), overwrite=TRUE)

# calculate curl method 2
uwnd <- stack('ccmp_uwnd.nc')
vwnd <- stack('ccmp_vwnd.nc')
# uwnd <- stack('ascat_uwnd.nc')
# vwnd <- stack('ascat_vwnd.nc')

# rename rasters
t <- names(uwnd)
tt <- as.numeric(sapply(t, function(x) str_split(x, 'X')[[1]][2]))
tt <- anydate(tt)
names(uwnd) <- tt
t <- names(vwnd)
tt <- as.numeric(sapply(t, function(x) str_split(x, 'X')[[1]][2]))
tt <- anydate(tt)
names(vwnd) <- tt

# calculate wind stress
rho_air = 1.2
cd = 1.2e-3
wmag = sqrt(uwnd^2+vwnd^2)
taux = rho_air*cd*wmag*uwnd
tauy = rho_air*cd*wmag*vwnd
dx <- res(taux)[1]
dy <- res(taux)[2]

# calculate curl
tauxx <- as.array(taux)
tauyy <- as.array(tauy)

#dTydx
dTydx <- array(dim = dim(tauyy))
for(i in 1:dim(tauyy)[3]){
  out <- apply(tauyy[,,i], 1, function(y){
    (lead(y) - lag(y)) / (2 * dx)
  }) # BY ROWS
  out <- t(out)
  dTydx[,,i] <- out
  
  # front edges
  front <- apply(tauyy[,1:2,i], 1, function(y){
    (lead(y) - y) / (2 * dx)
  })
  front <- t(front)
  dTydx[,1,i] <- front[,1]
  
  # back edges
  nc <- ncol(tauyy)
  back <- apply(tauyy[,(nc-1):nc,i], 1, function(y){
    (y - lag(y)) / (2 * dx)
  })
  back <- t(back)
  dTydx[,nc,i] <- back[,2]
}

#dTxdy
dTxdy <- array(dim = dim(tauxx))
for(i in 1:dim(tauxx)[3]){
  out <- apply(tauxx[,,i], 2, function(y){
    (lag(y) - lead(y)) / (2 * dy)
  }) # BY COLUMNS
  dTxdy[,,i] <- out
  
  # top edges
  front <- apply(tauxx[1:2,,i], 2, function(y){
    (lead(y) - y) / (2 * dy)
  })
  dTxdy[1,,i] <- front[1,]
  
  # back edges
  nr <- nrow(tauxx)
  back <- apply(tauxx[(nr-1):nr,,i], 2, function(y){
    (y - lag(y)) / (2 * dy)
  })
  dTxdy[nr,,i] <- back[2,]
}

# calculate curl
curlval <- dTydx - dTxdy
curl <- tauy
values(curl) <- curlval

# writeRaster(curlr, filename=paste("wcurl_5-8day", extentfn, '.grd', sep = ''), overwrite=TRUE)


# for(i in 1:dim(tauyy)[3]){
# out <- apply(tauyy[,,i], 1, function(y){
#   ((4/3) * (lead(y) - lag(y)) - (1/3) * (lead(y, 2) - lag(y, 2))/2) / (2 * 0.25)
# })
# out <- t(out)
# dTydy[,,i] <- out
# }
# 
# dTxdx <- array(dim = dim(tauxx))
# for(i in 1:dim(tauxx)[3]){
#   out <- apply(tauxx[,,i], 2, function(x){
#     ((4/3) * (lag(x) - lead(x)) - (1/3) * (lag(x, 2) - lead(x, 2))/2) / (2 * 0.25)
#   })
#   dTxdx[,,i] <- out
# }


# climatologies -----------------------------------------------------------
# load data
dir(pattern = '.grd')
glob <- stack("globoc_8day_136_141_-36_-43.grd")
qscurl <- stack("quikscatcurl_mday_136_141_-36_-43.grd" )
ascurl <- stack("ascatcurl_mday_136_141_-36_-43.grd")

# calculate climatology 
clim <- function(r, vn){
  # r = raster
  # vn = variable name string
  rasdates <- function(x) {
    t <- names(x)
    t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
    t <- as.character(gsub('\\.', '-', t))
    return(as.Date(t))
  }
  rclm <- stackApply(r, lubridate::month(rasdates(r)), mean)
  names(rclm) <- unique(lubridate::month(rasdates(r), label = T, abbr = T))
  writeRaster(rclm, filename = sprintf('%s_%s_%s.%s', 'climatology', vn, extentfn, 'grd'), overwrite=TRUE)
  return(rclm)
}
clim(glob, 'chlglob')
clim(qscurl, 'qscurl')
clim(ascurl, 'ascurl')


# antarctic oscillation AKA southern annular mode  --------------------------------------------------
library(stringr)
download.file('ftp://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.aao.index.b790101.current.ascii',
              destfile = 'aao.ascii')

aao <- tbl_df(read.csv('aao.ascii', col.names = 'v1'))
aao$v1 <- as.character(aao$v1)
aao <- apply(as.matrix(aao), 1, function(x){
  s <- str_split(x, pattern ="\\s+")
  out <- data.frame(year = s[[1]][1], month = s[[1]][2], day = s[[1]][3], v1 = s[[1]][4])
  return(out)
})

aao <- tbl_df(do.call('rbind', aao))
test <- aao %>% mutate_all(funs(as.character(.)))
test <- test %>% mutate_all(funs(as.numeric(.)))
aao <- test
aao <- aao %>% mutate(date = as.Date(sprintf('%s-%s-%s', year, month, day)))
save(aao, file = 'aao.Rdata')


# SOI  --------------------------------------------------------------------
library(stringr)
library(tidyr)
# download.file('https://www.esrl.noaa.gov/psd/data/correlation/soi.data',
#               destfile = 'soi.data')
download.file('ftp://ftp.bom.gov.au/anon/home/ncc/www/sco/soi/soiplaintext.html',
              destfile = './raw enviro data/soi_bom.data')
soi <- tbl_df(read.csv('./raw enviro data/soi_bom.data', col.names = 'v1', skip = 11))
soi <- soi[-c(143:146),]
soi$v1 <- as.character(soi$v1)
# x= soi[1,]
s <- split(soi$v1, 1:nrow(soi)) %>% purrr::map(~str_split(., "\\s+")[[1]][-1])
soi <- s %>% purrr::map(function(x){
  data.frame(year = x[1], 
             Jan = x[2],
             Feb = x[3],
             Mar = x[4],
             Apr = x[5],
             May = x[6],
             Jun = x[7],
             Jul = x[8],
             Aug = x[9],
             Sep = x[10],
             Oct = x[11],
             Nov = x[12],
             Dec = x[13])
}) %>% 
  reduce(bind_rows)
soi <- soi %>% mutate_all(funs(as.numeric(.)))
soi <- tbl_df(soi)
test <- soi %>% gather(month, v1, 2:13) %>% 
  mutate(day = 1, date = as.Date(sprintf('%s-%s-%s', year, month, day), format = '%Y-%b-%d'))
soi <- test %>%  arrange(date) %>% filter(year < 2017)
save(soi, file = './raw enviro data/soi_bom.RData')


# ONI ---------------------------------------------------------------------
library(stringr)
library(tidyr)
download.file('https://www.esrl.noaa.gov/psd/data/correlation/oni.data',
              destfile = 'oni.data')
oni <- tbl_df(read.csv('oni.data', col.names = 'v1'))
oni <- oni[-c(68:78),]
oni$v1 <- as.character(oni$v1)
x= oni[1,]
oni <- apply(as.matrix(oni), 1, function(x){
  s <- list(str_split(x, pattern ="\\s+")[[1]][-1])
  out <- data.frame(year = s[[1]][1], 
                    Jan = s[[1]][2],
                    Feb = s[[1]][3],
                    Mar = s[[1]][4],
                    Apr = s[[1]][5],
                    May = s[[1]][6],
                    Jun = s[[1]][7],
                    Jul = s[[1]][8],
                    Aug = s[[1]][9],
                    Sep = s[[1]][10],
                    Oct = s[[1]][11],
                    Nov = s[[1]][12],
                    Dec = s[[1]][13])
  return(out)
})

oni <- tbl_df(do.call('rbind', oni))
test <- oni %>% gather(month, v1, 2:13) %>% 
  mutate(day = 1, date = as.Date(sprintf('%s-%s-%s', year, month, day), format = '%Y-%b-%d'))
oni <- test %>%  arrange(date)
oni$v1 <- as.numeric(oni$v1)
save(oni, file = 'oni.Rdata')


# Wind stress upwelling ---------------------------------------------------
## Raadtools method
# dates <- seq(as.POSIXct('1997-01-01 00:00:00', tz = 'GMT'), as.POSIXct('2017-12-31 18:00:00', 'GMT'), by = '6 hours')
# w <- tibble(date = dates)
# w$lon <- 136.875
# w$lat <- -36.2
# w$u <- raadtools::extract(readwind, w[, c("lon", "lat", "date")], uonly=TRUE)   
# w$v <- raadtools::extract(readwind, w[, c("lon", "lat", "date")], vonly=TRUE)   
# save(w, file = './extracted enviro data/wind_6H_neptuneIs_1997-2017.RData') # unsure if raadtools wind data is wind at 10 m data
# load('./extracted enviro data/wind_6H_neptuneIs_1997-2017.RData')

## Download from NCEP source website method
dest <- './raw enviro data/ncep 6h 10m wind/'
sf_uwnd <- paste('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/uwnd.10m.gauss.',
                 seq(1997,2017, by = 1), '.nc', sep ='')
sf_vwnd <- paste('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/vwnd.10m.gauss.',
                 seq(1997,2017, by = 1), '.nc', sep ='')
wd <- "./ncep 6-h 10-m wind/"
map(sf_uwnd, ~download.file(., destfile = paste(dest, str_split(., '/')[[1]][7], sep = '')))
map(sf_vwnd, ~download.file(., destfile = paste(dest, str_split(., '/')[[1]][7], sep = '')))

# create raster stacks
uwnd <- stack(dir('./raw enviro data/ncep 6h 10m wind/', pattern = 'uwnd', full.names = T))
vwnd <- stack(dir('./raw enviro data/ncep 6h 10m wind/', pattern = 'vwnd', full.names = T))

# crop to focus location
# extent <- extent(136.12, 136.2, -35.5, -35.3) # neptune island 
extent <- extent(140, 141, -40, -39) # bonney coast
uwnd1 <- crop(uwnd, extent)
vwnd1 <- crop(vwnd, extent)

# convert to dataframe
dt <- seq(as.POSIXct('1997-01-01 00:00:00', tz = 'GMT'), as.POSIXct('2017-12-31 18:00:00', tz = 'GMT'), by = '6 hours')
w0 <- purrr::map(list(uwnd1, vwnd1), function(x){
  df <- data.frame(rasterToPoints(x))
  df <- df %>% tidyr::gather("date", 'v1', 3:ncol(df))
  df$date <- dt
  tbl_df(df)
})

w1 <- left_join(w0[[1]], w0[[2]], by = c('x', 'y', 'date'))
colnames(w1) <- c('lon', 'lat', 'date', 'u', 'v')
w <- w1
# save(w, file = './extracted enviro data/neptuneIs_ncep_6h_10m_wind_1997-2017.RData')
# save(w, file = './extracted enviro data/bonneycoast_ncep_6h_10m_wind_1997-2017.RData')


# calculate mean daily wind stress
w1 <- w %>% 
  mutate(day = as.Date(date)) %>%
  group_by(day) %>% 
  summarise(u = mean(u, na.rm = T), v = mean(v, na.rm = T)) %>% 
  rename(date = day)

# calculate wind speed and direction from u and v components
source('./scripts/wind_dir_speed.R')
w1 <- w1 %>% group_by(date) %>%  mutate(w = wind.dir.speed(u, v)[1], wdir = wind.dir.speed(u, v)[2])

# calculate upwelling wind index
cd <- 1.2 * 10^-3
pho <- 1.22
coast_angle = 315
# beta = (360 - coast_angle) * pi/180
# w1 <- w1 %>% mutate(w = sqrt(u^2 + v^2), alpha = atan(v/u), wind_index = pho * cd * w^2 *sin(alpha - beta)) # possibly wrong calculation
w1 <- w1 %>% mutate(upwell_wind = pho * cd * w^2 * cos((wdir - coast_angle) * pi/180))
uw <- w1
save(uw, file = "./extracted enviro data/bonneycoast_10m_upwelling_wind_index_1997-2017.RData")

w1 %>% mutate(year = year(date)) %>% 
  group_by(year) %>% 
  summarise(mean = mean(windstress)) %>% 
  filter(year < 2007) %>% 
  ggplot(aes(year, mean)) + geom_point() + geom_line()

w1 %>% 
  mutate(year = year(date), month = month(date), period = ifelse(month < 11 & month > 4, 'nonupwelling', 'upwelling')) %>% 
  filter(year < 2007) %>%
  group_by(period) %>% 
  summarise(mean = mean(upwell_wind)) %>% 
  ggplot(aes(period, mean)) + 
  geom_col()

# SST Anomaly (errdap) ----------------------------------------------------
fn <- paste('ssta1day', extentfn, '1997-01-01', '1998-01-01', '.nc', sep = '_')
url <- map(1997:2017, function(t1){
  paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst2Agg.nc?anom[(', t1, '-01-01T00:00:00Z):1:(', t1, '-12-31T00:00:00Z)][(0.0):1:(0.0)][(-46):1:(-36)][(125.5):1:(145.5)]', sep = '')
}) %>% unlist()
fn <- paste('ssta', 1997:2017, '.nc', sep = '')

map2(url, fn, ~download.file(.x, destfile = .y))

# combine all .nc files into raster stack
options("scipen" = 10000000, "digits" = 10)
r <- stack(fn)

# Convert unix time to dates 
t <- as.numeric(sapply(names(r), function(x) str_split(x, 'X')[[1]][2]))
t <- anydate(t)

# resample raster to match sst and ssha (0.25)
r1 <- resample(r, sst, method="bilinear")

# join seawif and modis rasters
names(r1) <- t
r <- setZ(r1, t)
ssta <- r1
fn <- paste('ssta1day', extentfn, '1997-01-01', '2017-12-31', '(NOAA)','.grd', sep = '_')
writeRaster(ssta, filename = fn, overwrite=TRUE)
              


# SST Anomaly (manual)--------------------------------------------------------------
sst <- brick("daily_sst_129.5_145.5_-36_-46_1997-01-01_2018-01-01_.grd" )

# create dates of layers
t <- seq(as.Date('1997-01-01'), length = nlayers(sst), by='day')

# remove leap days
sst <- sst[[-c(which(month(t) == 2 & day(t) == 29))]]
t <- t[-c(which(month(t) == 2 & day(t) == 29))]
sst <- setZ(sst, t)

# calculate anomaly
sst_global_mean <- mean(sst)
ssta_daily <- sst - sst_global_mean
names(ssta_daily) <- t
ssta_mc <- stackApply(ssta_daily, month(t), mean)
ssta_daily <- setZ(ssta_daily, t)
fn <- paste('ssta1day', extentfn, '1997-01-01', '2017-12-31', '.grd', sep = '_')

ssta <- ssta_daily

writeRaster(ssta, filename = fn, overwrite = TRUE)


# bathymetry --------------------------------------------------------------
bathy <- readbathy(xylim = myextent)
r1 <- resample(bathy, sst, method="bilinear")

fn <- paste('bathy', extentfn, '.grd', sep = '_')
writeRaster(r1 , filename = fn, overwrite = TRUE)








