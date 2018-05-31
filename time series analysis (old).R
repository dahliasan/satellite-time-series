### extract time series from satellite data ###

library(raadtools)
library(dplyr)
library(tidyr)
library(rgdal)
library(ggplot2)
library(viridis)
library(xts)
library(zoo)
library(lubridate)
library(strucchange)
library(bfast)
# library(WaveletComp)
library(stringr)
library(anytime)
library(forecast)
library(gridExtra)
library(grid)


# load raster files
ssharas <- stack("daily_ssha_136_141_-36.5_-43.grd")
chlras <- stack("weekly_chl_136_141_-36.5_-43.grd")
sstras <- stack("daily_sst_136_141_-36.5_-43.grd")
curlras <- stack("windstresscurl_1day_136_141_-36.5_-43.grd")
scuras <- stack("oscar_um_5day_136_141_-36.5_-43.grd")
scvras <- stack("oscar_vm_5day_136_141_-36.5_-43.grd")

# set focus locations 
x <- 138.5
y <- -36.75
e1 <- extent(x, x +0.01, y, y + 0.01)

x <- 138.5
y <- -39
e2 <- extent(x, x+0.01, y, y + 0.01)

x <- 138.5
y <- -41
e3 <- extent(x, x+0.01, y, y + 0.01)

extfn <- paste(e1@xmin, e1@ymin, e2@ymin, e3@ymin, sep = '_')

# extract raster for specific location
ssharasp1 <- crop(stack(ssharas), e1)
chlrasp1 <- crop(stack(chlras), e1)
sstrasp1 <- crop(stack(sstras), e1)
curlrasp1 <- crop(stack(curlras), e1)
scurasp1 <- crop(stack(scuras), e1)
scvrasp1 <- crop(stack(scvras), e1)

ssharasp2 <- crop(stack(ssharas), e2)
chlrasp2 <- crop(stack(chlras), e2)
sstrasp2 <- crop(stack(sstras), e2)
curlrasp2 <- crop(stack(curlras), e2)
scurasp2 <- crop(stack(scuras), e2)
scvrasp2 <- crop(stack(scvras), e2)

ssharasp3 <- crop(stack(ssharas), e3)
chlrasp3 <- crop(stack(chlras), e3)
sstrasp3 <- crop(stack(sstras), e3)
curlrasp3 <- crop(stack(curlras), e3)
scurasp3 <- crop(stack(scuras), e3)
scvrasp3 <- crop(stack(scvras), e3)

# get dfs ----------------------------------------------------

sstl <- list(sstrasp1, sstrasp2, sstrasp3)
sshal <- list(ssharasp1, ssharasp2, ssharasp3)
chll <- list(chlrasp1, chlrasp2, chlrasp3)
curll <- list(curlrasp1, curlrasp2, curlrasp3)
scul <- list(scurasp1, scurasp2, scurasp3)
scvl <- list(scvrasp1, scvrasp2, scvrasp3)

# convert raster to df 
ras2df <- function(x) {
  df <- as.data.frame(x)
  df <- tbl_df(t(df))
  # t <- seq(as.Date(startdate), length = nlayers(x), by = timescale)
  t <- names(x)
  t <- sapply(t, function(x) str_split(x, 'X')[[1]][2])
  t <- gsub('\\.', '-', t)
  t <- as.Date(t)
  df$date <- t
  return(df)
}
dfsstl <- lapply(sstl, ras2df)
dfsshal <- lapply(sshal, ras2df)
dfchll <- lapply(chll, ras2df)
dfcurll <- lapply(curll, ras2df)
dfscul <- lapply(scul, ras2df)
dfscvl <- lapply(scvl, ras2df)

# chl: pad and interpolate ####
varnames <- 'chl'
locfac <- c('a', 'b', 'c')

for(j in 1:length(locfac)) {
  # pad time series from 1997 - 2002
  x <- dfchll[[j]]
  tm <- seq(as.Date('1997-12-29'), length = 6, by = '1 year')
  tmd <- as.data.frame(list(date = tm))
  x <- tbl_df(merge(tmd, x, all=T))
  x %>% group_by(year(date)) %>% summarise(n = n())
  x <- x[order(x$date),]

  # assign location factor
  x$loc <- locfac[j]  
  
  dfchll[[j]] <- x
  
}

# interpolate missing values 
dfchll <- lapply(dfchll, function(x) 
{x$interp <- approx(x = x$date, y = x$V1, xout = x$date)$y
# if(length(which(is.na(x$interp))) > 0){
#   x <- x[-which(is.na(x$interp)),]} # remove NAs that still exist after interpolation
return(x)})

# correct colnames
dfchll <- lapply(dfchll, function(x) {
  colnames(x)[2] <- varnames
  return(x)})




# curl: pad and interpolate ------------------------------------------------------
varnames <- 'curl'
locfac <- c('a', 'b', 'c')

# make complete data frames 
for(j in 1:length(locfac)) {
  x <- dfcurll[[j]]
  # pad time series
  padts <- function(x) {
    tr <- as.Date(range(x$date))
    td <- x$date[2] - x$date[1]
    ta <- seq(tr[1], tr[2], by = td)
    tadf <- data.frame(list(date=ta))
    xa <- tbl_df(merge(tadf, x, all=T))
    return(xa)
  }
  x <- padts(x)
  x <- x[order(x$date),]
 
  # assign location factor
  x$loc <- locfac[j]  
  
  dfcurll[[j]] <- x
}

# interpolate missing values 
dfcurll <- lapply(dfcurll, function(x) 
{x$interp <- approx(x = x$date, y = x$V1, xout = x$date)$y
# if(length(which(is.na(x$interp))) > 0){
#   x <- x[-which(is.na(x$interp)),]} # remove NAs that still exist after interpolation
return(x)})
  
dfcurll <- lapply(dfcurll, function(x) {
  colnames(x)[2] <- varnames
  return(x)})

dfcurll[[j]] %>% group_by(year(date)) %>% summarise(n = n())


# the rest: interpolate ---------------------------------------------------

varnames <- c('sst', 'ssha', 'um', 'vm')
locfac <- c('a', 'b', 'c')
dflist <- list(dfsstl, dfsshal, dfscul, dfscvl)

# make complete data frames 
for(i in 1:length(varnames)) {
  for(j in 1:length(locfac)) {
    # pad time series 
    # padts <- function(x) {
    #   tr <- as.Date(range(x$date))
    #   td <- x$date[2] - x$date[1]
    #   ta <- seq(tr[1], tr[2], by = td)
    #   tadf <- data.frame(list(date=ta))
    #   xa <- tbl_df(merge(tadf, x, all=T))
    #   return(xa)
    # }
    # dflist[[i]][[j]] <- padts(dflist[[i]][[j]])
    # 
    # assign location factor
    dflist[[i]][[j]]$loc <- locfac[j]  
  }
  
  # interpolate missing values 
  dflist[[i]] <- lapply(dflist[[i]], function(x) 
  {x$interp <- approx(x = x$date, y = x$V1, xout = x$date)$y
  # if(length(which(is.na(x$interp))) > 0){
  #   x <- x[-which(is.na(x$interp)),]} # remove NAs that still exist after interpolation
  return(x)})
}

for(i in 1:length(varnames)) {
  dflist[[i]] <- lapply(dflist[[i]], function(x) {
    colnames(x)[1] <- varnames[i]
    return(x)})
}


# create dataframes -------------------------------------------------------
sst <- do.call(rbind, dflist[[1]])
ssha <- do.call(rbind, dflist[[2]])
scu <- do.call(rbind, dflist[[3]])
scv <- do.call(rbind, dflist[[4]])
curl <- do.call(rbind, dfcurll)
chl <- do.call(rbind, dfchll)
sc <- scu
sc$vm <- scv$vm
sc$interp <- with(sc, sqrt(um^2 + vm^2))

chl$loc <- as.factor(chl$loc)
sst$loc <- as.factor(sst$loc)
ssha$loc <- as.factor(ssha$loc)
curl$loc <- as.factor(curl$loc)
scu$loc <- as.factor(scu$loc)
scv$loc <- as.factor(scv$loc)
sc$loc <- as.factor(sc$loc)

# inspect dataframe 
sc %>% group_by(loc, year(date)) %>% summarise(n = n())
chl %>% group_by(loc, year(date)) %>% summarise(n = n())
sst %>% group_by(loc, year(date)) %>% summarise(n = n())
ssha %>% group_by(loc, year(date)) %>% summarise(n = n())
curl %>% group_by(loc, year(date)) %>% summarise(n = n())

# remove leap days for daily datasets
sst <- sst[-with(sst, which(month(date) == 2 & day(date) == 29)),]
ssha <- ssha[-with(ssha, which(month(date) == 2 & day(date) == 29)),]
curl <- curl[-with(curl, which(month(date) == 2 & day(date) == 29)),]

# climatology boxplots ----------------------------------------------------
tiff(file = paste("ssha_boxplot_", extfn, ".tiff", sep = ''), width=7.5, height=7.5, units="in", res=300)

ggplot(data = ssha, aes(x = as.factor(month(date)), y = interp)) + 
  geom_boxplot(aes(fill = loc), outlier.size = 0.5) + 
  labs(x = 'month', y = 'ssha') + 
  theme(legend.position = "top")
dev.off()

tiff(file = paste("sst_boxplot_", extfn, ".tiff", sep = ''), width=7.5, height=7.5, units="in", res=300)
ggplot(data = sst, aes(x = as.factor(month(date)), y = interp)) + 
  geom_boxplot(aes(fill = loc), outlier.size = 0.5) + 
  labs(x = 'month', y = 'sst') + 
  theme(legend.position = "top")
dev.off()

tiff(file = paste("chl_boxplot_", extfn, ".tiff", sep = ''), width=7.5, height=7.5, units="in", res=300)
ggplot(data = chl, aes(x = as.factor(month(date)), y = interp)) + 
  geom_boxplot(aes(fill = loc), outlier.size = 0.5) + 
  labs(x = 'month', y = 'chl') + 
  theme(legend.position = "top") +
  geom_hline(yintercept = 0.6, linetype = 'dashed')
dev.off()

tiff(file = paste("lgchl_boxplot_", extfn, ".tiff", sep = ''), width=7.5, height=7.5, units="in", res=300)
ggplot(data = chl, aes(x = as.factor(month(date)), y = log(interp))) + 
  geom_boxplot(aes(fill = loc), outlier.size = 0.5) + 
  labs(x = 'month', y = 'log(chl)') + 
  theme(legend.position = "top") +
  geom_hline(yintercept = log(0.6), linetype = 'dashed')
dev.off()

curls <- curl %>% 
  mutate(month = month(date, abbr = TRUE, label = TRUE), 
         year = year(date), 
         negcurl = as.integer(interp < 0)) %>% 
  group_by(loc, year, month) %>% 
  summarise(negcurl = sum(negcurl)) %>% group_by(loc, month) %>% 
  summarise(nnc = mean(negcurl), sd = sd(negcurl))

tiff(file = paste("curlupwellingevents_barplot_", extfn, ".tiff", sep = ''), width=7.5, height=7.5, units="in", res=300)
dodge <- position_dodge(width=0.9)
ggplot(data = curls, aes( x = month, y = nnc, fill = loc)) + 
  geom_bar(stat = 'identity', position = dodge) + 
  geom_errorbar(aes(ymax = nnc + sd, ymin = nnc -sd), alpha= 0.5, width = 0.4, position = dodge) + 
  labs(x = 'month', y = 'negative curl events (upwelling)', title = extfn) + 
  theme(legend.position = "top") + 
  theme_minimal() + 
  scale_fill_grey()
dev.off()


scs <- sc %>% mutate(month = month(date, abbr = TRUE, label = TRUE), 
                     year = year(date)) %>% 
  group_by(loc, month) %>% 
  summarise(um = mean(um), vm = mean(vm)) %>% 
  mutate(y = as.integer(loc), x = as.integer(month))

tiff(file = paste("winddirection_arrowplot_", extfn, ".tiff", sep = ''), width=7.5, height=7.5, units="in", res=300)
scaler <- 5
ggplot(data = scs, 
       aes(x = month, y = loc, xend = x + um * scaler, yend = y + vm * scaler)) + 
  geom_segment(arrow = arrow(length = unit(0.1, 'cm'))) + 
  labs(x = 'month', y = 'location', title = 'a = -36.75, b = -39, c = -41')
dev.off()


# create ts objects -------------------------------------------------------
# check if ts are in sequence 
which(!(chl %>% group_by(loc) %>% mutate(ddiff = date - lead(date)))$ddiff == -8)
which(!(sst %>% group_by(loc) %>% mutate(ddiff = date - lead(date)))$ddiff == -1)
which(!(ssha %>% group_by(loc) %>% mutate(ddiff = date - lead(date)))$ddiff == -1)
which(!(curl %>% group_by(loc) %>% mutate(ddiff = date - lead(date)))$ddiff == -1)
which(!(sc %>% group_by(loc) %>% mutate(ddiff = date - lead(date)))$ddiff == -5)

# create ts for ssta and ssha

sstp1ts <- ts(dflist[[1]][[1]]$sst, start = c(1997,1), frequency = 365)
sstp2ts <- ts(dflist[[1]][[2]]$sst, start = c(1997,1), frequency = 365)
sstp3ts <- ts(dflist[[1]][[3]]$sst, start = c(1997,1), frequency = 365)
sshap1ts <- ts(dflist[[2]][[1]]$ssha, start = c(1997,1), frequency = 365)
sshap2ts <- ts(dflist[[2]][[2]]$ssha, start = c(1997,1), frequency = 365)
sshap3ts <- ts(dflist[[2]][[3]]$ssha, start = c(1997,1), frequency = 365)
chlp1ts <- ts(dflist[[3]][[1]]$interp, start = c(1997,33), frequency = 46)
chlp2ts <- ts(dflist[[3]][[2]]$interp, start = c(1997,33), frequency = 46)
chlp3ts <- ts(dflist[[3]][[3]]$interp, start = c(1997,33), frequency = 46)
curlp1ts <- ts(dfcurll[[1]]$interp, start = c(1999,202), frequency = 365)
curlp2ts <- ts(dfcurll[[2]]$interp, start = c(1999,202), frequency = 365)
curlp3ts <- ts(dfcurll[[3]]$interp, start = c(1999,202), frequency = 365)



# create 8 day ts with same date range ------------------------------------
mindate <- "2000-01-01"
maxdate <- min(c(range(ssha$date)[2], range(sst$date)[2], range(chl$date)[2],
                 range(curl$date)[2]))
chl2 <- chl %>%  group_by(loc) %>% filter(date >= mindate & date <= maxdate)
chl2$chl <- chl2$interp
chl2 <- chl2[,-4]
sst2 <- sst %>%  group_by(loc) %>% filter(date >= mindate & date <= maxdate ) %>%
  mutate(yr = year(date)) %>% group_by(loc, yr) %>% 
  mutate(date = cut(date, 46)) %>% group_by(loc, date) %>% 
  summarise(sst = mean(interp))
ssha2 <- ssha %>%  group_by(loc) %>% filter(date >= mindate & date <= maxdate ) %>%
  mutate(yr = year(date)) %>% group_by(loc, yr) %>% 
  mutate(date = cut(date, 46)) %>% group_by(loc, date) %>% 
  summarise(ssha = mean(interp))
curl2 <- curl %>%  group_by(loc) %>% filter(date >= mindate & date <= maxdate ) %>%
  mutate(yr = year(date)) %>% group_by(loc, yr) %>% 
  mutate(date = cut(date, 46)) %>% group_by(loc, date) %>% 
  summarise(curl = mean(interp))
sc2 <- sc %>%  group_by(loc) %>% filter(date >= mindate & date <= maxdate ) %>%
  mutate(yr = year(date)) %>% group_by(loc, yr) %>% 
  mutate(date = cut(date, 46)) %>% group_by(loc, date) %>% 
  summarise(sc = mean(interp))

chl2 %>% group_by(loc, year(date)) %>% summarise(n = n())
sst2 %>% group_by(loc, year(date)) %>% summarise(n = n())
ssha2 %>% group_by(loc, year(date)) %>% summarise(n = n())
curl2 %>% group_by(loc, year(date)) %>% summarise(n = n())
sc2 %>% group_by(loc, year(date)) %>% summarise(n = n())


sstp1ts <- ts(sst2$sst[sst2$loc == 'a'], start = c(2000,1), frequency = 46)
sstp2ts <- ts(sst2$sst[sst2$loc == 'b'], start = c(2000,1), frequency = 46)
sstp3ts <- ts(sst2$sst[sst2$loc == 'c'], start = c(2000,1), frequency = 46)
sshap1ts <- ts(ssha2$ssha[ssha2$loc == 'a'], start = c(2000,1), frequency = 46)
sshap2ts <- ts(ssha2$ssha[ssha2$loc == 'b'], start = c(2000,1), frequency = 46)
sshap3ts <- ts(ssha2$ssha[ssha2$loc == 'c'], start = c(2000,1), frequency = 46)
chlp1ts <- ts(log(chl2$chl[chl2$loc == 'a']), start = c(2000,1), frequency = 46)
chlp2ts <- ts(log(chl2$chl[chl2$loc == 'b']), start = c(2000,1), frequency = 46)
chlp3ts <- ts(log(chl2$chl[chl2$loc == 'c']), start = c(2000,1), frequency = 46)
curlp1ts <- ts(curl2$curl[curl2$loc == 'a'], start = c(2000,1), frequency = 46)
curlp2ts <- ts(curl2$curl[curl2$loc == 'b'], start = c(2000,1), frequency = 46)
curlp3ts <- ts(curl2$curl[curl2$loc == 'c'], start = c(2000,1), frequency = 46)
scp1ts <- ts(sc2$sc[sc2$loc == 'a'], start = c(2000,1), frequency = 46)
scp2ts <- ts(sc2$sc[sc2$loc == 'b'], start = c(2000,1), frequency = 46)
scp3ts <- ts(sc2$sc[sc2$loc == 'c'], start = c(2000,1), frequency = 46)

# BFAST package -----------------------------------------------------------
# sstap1fit <- bfast(sstap1ts, season="harmonic", max.iter=1)
# sstap2fit <- bfast(sstap2ts, season="harmonic", max.iter=1)
# sstap3fit <- bfast(sstap3ts, season="harmonic", max.iter=1)
# sshap1fit <- bfast(sshap1ts, season="harmonic", max.iter=1)
# sshap2fit <- bfast(sshap2ts, season="harmonic", max.iter=1)
# sshap3fit <- bfast(sshap3ts, season="harmonic", max.iter=1)
# chlp1fit <- bfast(chlp1ts, season="harmonic", max.iter=1)
# chlp2fit <- bfast(chlp2ts, season="harmonic", max.iter=1) 
# chlp3fit <- bfast(chlp3ts, season="harmonic", max.iter=1)

# get extent filenames
fne1 <- paste(e1@xmin, e1@ymin)
fne2 <- paste(e2@xmin, e2@ymin)
fne3 <- paste(e3@xmin, e3@ymin)

# fit 8 day ts 
sstp1fit2 <- bfast(sstp1ts + min(sst2$sst), season="harmonic", max.iter=1)
sstp2fit2 <- bfast(sstp2ts + min(sst2$sst), season="harmonic", max.iter=1)
sstp3fit2 <- bfast(sstp3ts + min(sst2$sst), season="harmonic", max.iter=1)
sshap1fit2 <- bfast(sshap1ts + min(ssha2$ssha), season="harmonic", max.iter=1)
sshap2fit2 <- bfast(sshap2ts + min(ssha2$ssha), season="harmonic", max.iter=1)
sshap3fit2 <- bfast(sshap3ts + min(ssha2$ssha), season="harmonic", max.iter=1)
chlp1fit2 <- bfast(chlp1ts, season="harmonic", max.iter=1)
chlp2fit2 <- bfast(chlp2ts, season="harmonic", max.iter=1) 
chlp3fit2 <- bfast(chlp3ts, season="harmonic", max.iter=1)
curlp1fit2 <- bfast(curlp1ts + min(curl2$curl), season="harmonic", max.iter=1)
curlp2fit2 <- bfast(curlp2ts + min(curl2$curl), season="harmonic", max.iter=1) 
curlp3fit2 <- bfast(curlp3ts + min(curl2$curl), season="harmonic", max.iter=1)
scp1fit2 <- bfast(scp1ts, season="harmonic", max.iter=1)
scp2fit2 <- bfast(scp2ts, season="harmonic", max.iter=1) 
scp3fit2 <- bfast(scp3ts, season="harmonic", max.iter=1)



# plot all the different available plots

bfastggplot <- function(x, title = NULL) {
  niter <- length(x$output)
  out <- x$output[[niter]]
  Trend.bp <- !x$nobp$Vt
  Tt <- data.frame(x = index(out$Tt), y = coredata(out$Tt))
  St <- data.frame(x = index(out$St), y = coredata(out$St))
  noise <- data.frame(x = index(out$Nt), y = coredata(out$Nt))
  Yt <- data.frame(x = index(x$Yt), y = coredata(x$Yt))
  
  margin = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),
                 plot.margin = unit(c(0,0,-1,0.5), "lines"))
  
  p1 <- ggplot(data = Yt, aes(x = x, y = y))+
    geom_line() +
    labs(y = 'Yt', x = '')+ margin
  
  p2 <- ggplot(data = St, aes(x = x, y = y))+
    geom_line() +
    labs(y = 'St', x = '')
  
  if(!x$nobp$Wt){
    sbp <- data.frame(x = breakdates(out$bp.Wt), y = St$y[out$bp.Wt$breakpoints])
    p2 <- p2 + 
      geom_text(data = sbp, aes(x = x, y = y, label = round(x, digits = 2)),size = 2) + 
      geom_vline(xintercept = sbp$x, linetype = 'dashed') + margin
  } else( p2 <- p2 + margin)
  
  
  p3 <- ggplot(data = Tt, aes(x = x, y = y))+
    geom_line() +
    labs(y = 'Tt', x = '')
  
  if(Trend.bp == TRUE){
    bp <- data.frame(x = breakdates(out$bp.Vt), y = Tt$y[out$bp.Vt$breakpoints])
    p3 <- p3 + 
      geom_text(data = bp, aes(x = x, y = y, label = round(x, digits = 2)),  size = 2) + 
      geom_vline(xintercept = bp$x, linetype = 'dashed') + margin
  } else( p3 <- p3 + margin)
  
  p4 <- ggplot(data = noise, aes(x = x, y = y))+
    geom_line() +
    labs(y = 'Nt', x = '') +
    scale_x_continuous(breaks = noise$x[seq(1,nrow(noise),46)], labels = seq(00, 16, 1)) + 
    theme(plot.margin = unit(c(0,0,0,0.5), "lines"))
  
  # Left justify plots
  # Source: http://stackoverflow.com/a/13295880/496488
  pl <- list(p1,p2,p3,p4)
  gl <- lapply(pl, function(x){
    ggplotGrob(x)
  })
  
  maxWidth =  unit.pmax (gl[[1]]$widths[2:5], gl[[2]]$widths[2:5], gl[[3]]$widths[2:5], gl[[4]]$widths[2:5])
  
  gl <- lapply(gl, function(x) {
    x$widths[2:5] <- as.list(maxWidth)
    return(x)
  })
  
  gg <- grid.arrange(gl[[1]],gl[[2]],gl[[3]],gl[[4]], ncol=1, 
                     top = title)
  return(gg)
}

# plot for location 1
tn <- c('sst','ssha', 'curl', 'chl')
bl <- list(sstp1fit2, sshap1fit2, curlp1fit2, chlp1fit2)
pl <- list()
for(i in 1:length(tn)){
  p <- list(bfastggplot(bl[[i]], paste(tn[i], fne1)))
  pl <- c(pl, p)
}
tiff(filename = paste('bfast ', fne1, '.tiff', sep = ''),  width=11.5, height=7, units= "in", res = 300)
grid.arrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]], nrow = 2, ncol = 2)
dev.off()

# plot for location 2
tn <- c('sst','ssha', 'curl', 'chl')
bl <- list(sstp2fit2, sshap2fit2, curlp2fit2, chlp2fit2)
pl <- list()
for(i in 1:length(tn)){
  p <- list(bfastggplot(bl[[i]], paste(tn[i], fne2)))
  pl <- c(pl, p)
}
tiff(filename = paste('bfast ', fne2, '.tiff', sep = ''),  width=11.5, height=7, units= "in", res = 300)
grid.arrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]], nrow = 2, ncol = 2)
dev.off()

# plot for location 3
tn <- c('sst','ssha', 'curl', 'chl')
bl <- list(sstp3fit2, sshap3fit2, curlp3fit2, chlp3fit2)
pl <- list()
for(i in 1:length(tn)){
  p <- list(bfastggplot(bl[[i]], paste(tn[i], fne3)))
  pl <- c(pl, p)
}
tiff(filename = paste('bfast ', fne3, '.tiff', sep = ''),  width=11.5, height=7, units= "in", res = 300)
grid.arrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]], nrow = 2, ncol = 2)
dev.off()

# cross correlation  ------------------------------------------------------

# make list of bfast objects
bl <- list(sstp1fit2, 
           sstp2fit2, 
           sstp3fit2, 
           sshap1fit2,
           sshap2fit2,
           sshap3fit2,
           chlp1fit2 ,
           chlp2fit2 ,
           chlp3fit2 ,
           curlp1fit2,
           curlp2fit2,
           curlp3fit2,
           scp1fit2 ,
           scp2fit2,
           scp3fit2)
vn <- c('sst', 'ssha', 'chl', 'curl', 'sc')
pn <- c(1, 2, 3)
ln <- NULL
for(i in 1:length(vn)){
  for(j in 1:length(pn)){
    nn <- paste(vn[i], pn[j], sep = '')
    ln <- c(ln,nn)
  }
}
names(bl) <- ln

# give each list element id 
for(i in 1:length(ln)){
  bl[[i]]$id <- ln[i]
}

# check if random componenet is stationary 
lapply(bl, function(x){
  plot(x, main = x$id, type = 'noise')
})

# ljung-box test < 0.05 is stationary 
lapply(bl, function(x) {
  Box.test(x$output[[1]]$Nt, lag = 20, type = 'Ljung-Box')
})

dcl <- list(bl$curl1, bl$curl2, bl$curl3)
dcl <- lapply(dcl, function(x){
  x$output[[1]]$Nt <- diff(x$output[[1]]$Nt)
  return(x)
})

lapply(dcl, function(x) {
  Box.test(x$output[[1]]$Nt, lag = 20, type = 'Ljung-Box')
})

sl <- c(bl[-c(10,11,12)], dcl) # stationary list
names(sl)[13:15] <- c('dc1', 'dc2', 'dc3')

# save(sl, file = 'stationary_bfast_list.Rdata')

# do ccf for predictor variables with chl
save(sl, file = 'cwt_result.Rdata')

tiff(filename = paste('ccf vs chl p1', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(2,2))
lapply(sl[c(1,4,10,13)], function(x){
  ccf(x = x$output[[1]]$Nt, y = sl$chl1$output[[1]]$Nt, 
      main = paste(x$id, 'vs.', 'chl1', sep = ' '))})
dev.off()

tiff(filename = paste('ccf vs chl p2', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(2,2))
lapply(sl[c(2,5,11,14)], function(x){
  ccf(x = x$output[[1]]$Nt, y = sl$chl2$output[[1]]$Nt, 
      main = paste(x$id, 'vs.', 'chl2', sep = ' '))})
dev.off()

tiff(filename = paste('ccf vs chl p3', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(2,2))
lapply(sl[c(3,6,12,15)], function(x){
  ccf(x = x$output[[1]]$Nt, y = sl$chl3$output[[1]]$Nt, 
      main = paste(x$id, 'vs.', 'chl3', sep = ' '))})
dev.off()

# continuous wt analysis ------------------------------------------------------------
library(biwavelet)

dtl <- lapply(sl, function(x){
  o <- cbind(index(x$output[[1]]$Wt), coredata(x$output[[1]]$Wt))
  return(o)
})

wtl <- lapply(dtl, function(x){
  o <- wt(x, dt = 365/46)
  return(o)
})

nn <- names(wtl)
for(i in 1:length(nn)){
  wtl[[i]]$id <- nn[i]
}

# save(wtl, file = 'cwt_result.Rdata')

# Plot power
tiff(filename = paste('sst_cwt', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(3,1))
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1) # Make room to the right for the color bar
lapply(wtl[1:3], function(x){
  plot(x, plot.cb = TRUE, plot.phase = FALSE, main = x$id)
})
dev.off()

tiff(filename = paste('ssha_cwt', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(3,1))
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1) # Make room to the right for the color bar
lapply(wtl[4:6], function(x){
  plot(x, plot.cb = TRUE, plot.phase = FALSE, main = x$id)
})
dev.off()

tiff(filename = paste('chl_cwt', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(3,1))
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1) # Make room to the right for the color bar
lapply(wtl[7:9], function(x){
  plot(x, plot.cb = TRUE, plot.phase = FALSE, main = x$id)
})
dev.off()

tiff(filename = paste('sc_cwt', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(3,1))
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1) # Make room to the right for the color bar
lapply(wtl[10:12], function(x){
  plot(x, plot.cb = TRUE, plot.phase = FALSE, main = x$id)
})
dev.off()

tiff(filename = paste('curl_cwt', '.tiff', sep = ''),  width=7, height=7, units= "in", res = 300)
par(mfrow = c(3,1))
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1) # Make room to the right for the color bar
lapply(wtl[13:15], function(x){
  plot(x, plot.cb = TRUE, plot.phase = FALSE, main = x$id)
})
dev.off()


# cross wavelet transformation --------------------------------------------
# obtain de-trended ts
satdf3 <- tbl_df(data.frame(chl = coredata(chlp3fit2$output[[1]]$Wt), 
                           ssha = coredata(sshap3fit2$output[[1]]$Wt),
                           ssta = coredata(sstap3fit2$output[[1]]$Wt)))
satdf3$date <- as.Date(chl2$wkdate[chl2$loc == 'c'])

satdf1 <- tbl_df(data.frame(chl = coredata(chlp1fit2$output[[1]]$Wt), 
                           ssha = coredata(sshap1fit2$output[[1]]$Wt),
                           ssta = coredata(sstap1fit2$output[[1]]$Wt)))
satdf1$date <- as.Date(chl2$wkdate[chl2$loc == 'a'])


## computation of cross-wavelet power and wavelet coherence: use detrended ts
my.wc = analyze.coherency(satdf3, c("chl","ssha"), dt = 1/45.625, dj = 1/100, loess.span = 0)
my.wc2 = analyze.coherency(satdf1, c("chl","ssha"), dt = 1/45.625, dj = 1/100, loess.span = 0)

wc.image(my.wc, timelab="time (year)", periodlab="period (year)",
         main="cross-wavelet power",
         legend.params=list(lab="cross-wavelet power levels"),
         show.date = TRUE)

wc.image(my.wc2, timelab="time (year)", periodlab="period (year)",
         main="cross-wavelet power",
         legend.params=list(lab="cross-wavelet power levels"),
         show.date = TRUE)

wc.avg(my.wc)

wc.image(my.wc, which.image="wc", timelab="time (days)", periodlab="period (days)",
         main="wavelet coherence",
         legend.params=list(lab="wavelet coherence levels", lab.line=3.5, label.digits=3))
## plot of average coherence:
wc.avg(my.wc, which.avg="wc", legend.coords="topleft")




# biwavelet package -------------------------------------------------------

# obtain de-trended ts
satdf3 <- data.frame(chl = coredata(chlp3fit2$output[[1]]$Wt), 
                     ssha = coredata(sshap3fit2$output[[1]]$Wt),
                     ssta = coredata(sstap3fit2$output[[1]]$Wt))
satdf3$date <- as.Date(chl2$wkdate[chl2$loc == 'c'])

satdf1 <- data.frame(chl = coredata(chlp1fit2$output[[1]]$Wt), 
                     ssha = coredata(sshap1fit2$output[[1]]$Wt),
                     ssta = coredata(sstap1fit2$output[[1]]$Wt))
satdf1$date <- as.Date(chl2$wkdate[chl2$loc == 'a'])

satdf2 <- data.frame(chl = coredata(chlp2fit2$output[[1]]$Wt), 
                     ssha = coredata(sshap2fit2$output[[1]]$Wt),
                     ssta = coredata(sstap2fit2$output[[1]]$Wt))
satdf2$date <- as.Date(chl2$wkdate[chl2$loc == 'b'])


# XWT ---------------------------------------------------------------------
# Cross-wavelet 
# chl vs ssha
xwt1 <- xwt(subset(satdf1, select = c("date", "ssha")),
            subset(satdf1, select = c("date", "chl")))
# Make room to the right for the color bar
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(xwt1, plot.cb = TRUE, plot.phase = TRUE)

xwt2 <- xwt(subset(satdf2, select = c("date", "ssha")),
            subset(satdf2, select = c("date", "chl")))
# Make room to the right for the color bar
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(xwt2, plot.cb = TRUE, plot.phase = TRUE)

xwt3 <- xwt(subset(satdf3, select = c("date", "ssha")),
            subset(satdf3, select = c("date", "chl")))
# Make room to the right for the color bar
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(xwt3, plot.cb = TRUE, plot.phase = TRUE)



xwt4 <- xwt(subset(satdf1, select = c("date", "chl")),
            subset(satdf3, select = c("date", "chl")))
# Make room to the right for the color bar
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(xwt4, plot.cb = TRUE, plot.phase = TRUE)


# WTC ---------------------------------------------------------------------
# Specify the number of iterations. The more, the better (>1000).  For the
# purpose of this tutorial, we just set it = 10
nrands = 10
wtc.AB <- wtc(subset(satdf1, select = c("date", "ssha")),
              subset(satdf1, select = c("date", "chl")),
              nrands = nrands)
# Plotting a graph
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "Wavelet Coherence: A vs B")


wtc.AB <- wtc(subset(satdf3, select = c("date", "ssha")),
              subset(satdf3, select = c("date", "chl")),
              nrands = nrands)
# Plotting a graph
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Scale", xlab = "Period", 
     plot.cb = TRUE, main = "Wavelet Coherence: A vs B")


# Adding grid lines
n = length(t1[, 1])
abline(v = seq(1, n, 46.625), h = 1:16, col = "brown", lty = 1, lwd = 1)

# Defining x labels
axis(side = 3, at = c(seq(min(t1[,1]), max(t1[,1]), )), labels = c(seq(1997, 2016, 1)))


