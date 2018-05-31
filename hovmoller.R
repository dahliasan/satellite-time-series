# make hovmoller plots
library(raster)
library(rasterVis)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(latticeExtra)
library(lattice)

# set extent for each section
e1 <- extent(136, 137, -43, -36.5)
e2 <- extent(137 + 0.01, 138, -43, -36.5)
e3 <- extent(138 + 0.01, 139, -43, -36.5)
e4 <- extent(139 + 0.01, 140, -43, -36.5)
el <- list(e1, e2, e3, e4)

hovclimplot <- function(r, varname, cont = TRUE) {
  # r = full region raster
  # vn = variable name
  # put each location raster stack in a list  
  ll <- list()
  for(i in 1:length(el)){
    r2 <- crop(r, el[[i]])
    ll <- c(ll, r2)
  }
  
  # setZ for each stack in the list
  ll <- lapply(ll, function(x) {
    t <- names(x)
    x <- setZ(x, t)
  })
  
  # reverse month arrangement and rename to abbr letters
  ll <- lapply(ll, function(x){
    t <- factor(as.factor(names(x)), levels = month.abb, ordered = TRUE)
    x <- x[[order(t)]]
    x <- x[[12:1]]
  })
  
  # plot hovmoller diagram
  # find scale limits 
  tmp <- lapply(ll, function(x){
    df <- tbl_df(data.frame(rasterToPoints(x)))
    df2 <- df %>% group_by(y) %>% summarise_at(colnames(df)[-c(1:2)], mean)
    min <- min(as.matrix(df2[,-1]), na.rm = TRUE)
    max <- max(as.matrix(df2[,-1]), na.rm = TRUE)
    return(c(min, max))
  })
  
  smin <- min(unlist(tmp), na.rm = TRUE)
  smax <- max(unlist(tmp), na.rm = TRUE)
  
  # set plot parameters
  p.at <- seq(smin, smax, length.out = 16)
  p.ckey <- list(at = p.at, labels = list(cex = 0.7), space = 'left', height = 0.5 )
  p.xlab <- list(label = '')
  p.ylab <- list(label = '')
  p.scales <- list(cex=0.7, x = list( tck = c(1,1)),
                   y = list(labels = rev(c('J','F','M','A','M','J','J','A','S','O','N','D'))))
  
  # make hovmoller plot
    llp <- lapply(ll, function(x){
      p <- hovmoller(x,
                     interpolate = TRUE,
                     xscale.components = xscale.raster.subticks,
                     par.settings = BuRdTheme,
                     xlab = p.xlab,
                     ylab = p.ylab,
                     scales = p.scales,
                     colorkey = p.ckey,
                     at = p.at,
                     add.contour = cont)
      return(p)})
  
  # save multiple plots to single tiff
  tiff(file= paste(varname, "climatology hovmoller 136-140E 36.5-43S.tiff"), width=6, height=11.5, units="in", res=300)
  
  p <- c(llp[[4]], llp[[3]], llp[[2]], llp[[1]], 
         layout = c(1,4))
  p2 <- update(p, par.settings = BuRdTheme(),
         ylab.right = list(label = c('139-140 E', '138-139 E', '137-138 E', '136-137 E'),
                           cex = 0.75, rot = 270),
         scales = list(alternating = FALSE, y = list(rot = -90)))
  
  print(p2)
  
  dev.off()
  
  return(llp)
}


# CLIMATOLOGY -------------------------------------------------------------


# chl glob ----------------------------------------------------------------
r <- log(stack("climatology_chlglob_136_141_-36_-43.grd" ))
chlglobplot <- hovclimplot(r, varname = 'chl glob') 

# qscat curl --------------------------------------------------------------
r <- stack("climatology_qscurl_136_141_-36_-43.grd")
qscurlplot <- hovclimplot(r, 'qscurl', FALSE)

# ascat curl --------------------------------------------------------------
r <- stack("climatology_ascurl_136_141_-36_-43.grd")
qscurlplot <- hovclimplot(r, 'ascurl')

# ssha  ---------------------------------------------------------
r <- stack("climatology_ssha_136_141_-36.5_-43.grd")
sshaplot <- hovclimplot(r, varname = 'ssha test')

# sshaclim <- stack("climatology_ssha_136_141_-36.5_-43.grd")
# 
# # put each location raster stack in a list  
# sshacl <- list()
# for(i in 1:length(el)){
#   r <- crop(sshaclim, el[[i]])
#   sshacl <- c(sshacl, r)
# }
# 
# # setZ for each stack in the list
# sshacl <- lapply(sshacl, function(x) {
#   t <- names(x)
#   x <- setZ(x, t)
# })
# 
# # reverse month arrangement and rename to abbr letters
# sshacl <- lapply(sshacl, function(x){
#   x <- x[[12:1]]
# })
# 
# # plot hovmoller diagram
# 
# p.at <- seq(-0.06, 0.16, .02)
# p.ckey <- list(at = p.at, labels = list(cex = 0.7), space = 'left', height = 0.5 )
# p.xlab <- list(label = '')
# p.ylab <- list(label = '')
# p.scales <- list(cex=0.7, x = list( tck = c(1,1)),
#                  y = list(labels = rev(c('J','F','M','A','M','J','J','A','S','O','N','D'))))
# 
# 
# sshacpl <- lapply(sshacl, function(x){
#   p <- hovmoller(x,
#             interpolate = TRUE,
#             xscale.components = xscale.raster.subticks,
#             par.settings = BuRdTheme,
#             xlab = p.xlab,
#             ylab = p.ylab,
#             scales = p.scales,
#             colorkey = p.ckey,
#             at = p.at
#             )
#   return(p)
# })
# 
# 
# tiff(file= "ssha test 2 climatology hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)
# 
# p <- c(sshacpl[[4]], sshacpl[[3]], sshacpl[[2]], sshacpl[[1]], 
#        layout = c(1,4))
# update(p, par.settings = BuRdTheme(),
#        ylab.right = list(label = c('139-140 E', '138-139 E', '137-138 E', '136-137 E'),
#                          cex = 0.75, rot = 270),
#        scales = list(alternating = FALSE, y = list(rot = -90)))
# 
# dev.off()


# sst   ---------------------------------------------------------
r <- stack("climatology_sst_136_141_-36.5_-43.grd")
sstplot <- hovclimplot(r, varname = 'sst test')

# ### 
# sstclim <- stack("climatology_sst_136_141_-36.5_-43.grd")
# 
# # put each location raster stack in a list  
# sstcl <- list()
# for(i in 1:length(el)){
#   r <- crop(sstclim, el[[i]])
#   sstcl <- c(sstcl, r)
# }
# 
# # setZ for each stack in the list
# sstcl <- lapply(sstcl, function(x) {
#   t <- names(x)
#   x <- setZ(x, t)
# })
# 
# # reverse month arrangement and rename to abbr letters
# sstcl <- lapply(sstcl, function(x){
#   x <- x[[12:1]]
# })
# 
# 
# # plot hovmoller diagram
# 
# p.at <- seq(9, 20, 1)
# p.ckey <- list(at = p.at, labels = list(cex = 0.7), space = 'left', height = 0.5 )
# p.xlab <- list(label = '')
# p.ylab <- list(label = '')
# p.scales <- list(cex=0.7, x = list( tck = c(1,1)),
#                  y = list(labels = rev(c('J','F','M','A','M','J','J','A','S','O','N','D'))))
# 
# sstcpl <- lapply(sstcl, function(x){
#   p <- hovmoller(x,
#                  interpolate = TRUE,
#                  xscale.components = xscale.raster.subticks,
#                  par.settings = BuRdTheme,
#                  xlab = p.xlab,
#                  ylab = p.ylab,
#                  scales = p.scales,
#                  colorkey = p.ckey,
#                  at = p.at
#   )
#   return(p)
# })
# 
# 
# tiff(file= "sst climatology hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)
# 
# p <- c(sstcpl[[4]], sstcpl[[3]], sstcpl[[2]], sstcpl[[1]], 
#        layout = c(1,4))
# update(p, par.settings = BuRdTheme(),
#        ylab.right = list(label = c('139-140 E', '138-139 E', '137-138 E', '136-137 E'),
#                          cex = 0.75, rot = 270),
#        scales = list(alternating = FALSE, y = list(rot = -90)))
# 
# dev.off()


# CHL   ---------------------------------------------------------
r <- log(stack("climatology_chl_136_141_-36.5_-43.grd"))
chlplot <- hovclimplot(r, varname = 'chl test')

# chlclim <- stack( "climatology_chl_136_141_-36.5_-43.grd")
# # chlclim <- setZ(chlclim)
# # rearraging chlclim to start from Jan. 
# chlclim <- stack(chlclim[[5]], chlclim[[6]],chlclim[[7]],chlclim[[8]],chlclim[[9]],
#                   chlclim[[10]],chlclim[[11]],chlclim[[12]],chlclim[[1]],chlclim[[2]],
#                   chlclim[[3]],chlclim[[4]])
# 
# 
# # put each location raster stack in a list  
# chlcl <- list()
# for(i in 1:length(el)){
#   r <- crop(chlclim, el[[i]])
#   chlcl <- c(chlcl, r)
# }
# 
# # prepare data for plotting
# chlcl <- lapply(chlcl, function(x) {
#   # log chl values
#   x <- calc(x, log)
#   # setZ for each stack in the list
#   t <- names(x)
#   x <- setZ(x, t)
#   # reverse month arrangement and rename to abbr letters
#   x <- x[[12:1]]
# })
# 
# 
#   
# # plot hovmoller diagram
# p.at <- seq(-2, 1.2, .2)
# p.ckey <- list(at = p.at, labels = list(cex = 0.7), space = 'left', height = 0.5 )
# p.xlab <- list(label = '')
# p.ylab <- list(label = '')
# p.scales <- list(cex=0.7, x = list( tck = c(1,1)),
#                  y = list(labels = rev(c('J','F','M','A','M','J','J','A','S','O','N','D'))))
# 
# 
# chlcpl <- lapply(chlcl, function(x){
#   p <- hovmoller(x,
#                  interpolate = TRUE,
#                  xscale.components = xscale.raster.subticks,
#                  par.settings = BuRdTheme,
#                  xlab = p.xlab,
#                  ylab = p.ylab,
#                  scales = p.scales,
#                  colorkey = p.ckey,
#                  at = p.at,
#                  add.contour = TRUE
#   )
#   return(p)
# })
# 
# 
# tiff(file= "chl climatology hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)
# 
# p <- c(chlcpl[[4]], chlcpl[[3]], chlcpl[[2]], chlcpl[[1]], 
#        layout = c(1,4))
# update(p, par.settings = BuRdTheme(),
#        ylab.right = list(label = c('139-140 E', '138-139 E', '137-138 E', '136-137 E'),
#                          cex = 0.75, rot = 270),
#        scales = list(alternating = FALSE, y = list(rot = -90)))
# 
# dev.off()


# sea current strength  --------------------------------------------------------------
sc <- stack("climatology_oscar_5day_136_141_-36.5_-43.grd")
sc <- sc[[c(5,6,7,8,9,10,11,12,1,2,3,4)]]


# put each location raster stack in a list  
sccl <- list()
for(i in 1:length(el)){
  r <- crop(sc, el[[i]])
  sccl <- c(sccl, r)
}

# setZ for each stack in the list
sccl <- lapply(sccl, function(x) {
  t <- names(x)
  x <- setZ(x, t)
})

# reverse month arrangement and rename to abbr letters
sccl <- lapply(sccl, function(x){
  x <- x[[12:1]]
})

# plot hovmoller diagram
man <- max(unlist(lapply(sccl, maxValue)))
miin <-  min(unlist(lapply(sccl, minValue)))
p.at <- seq(0.04, 0.2, .02)
p.ckey <- list(at = p.at, labels = list(cex = 0.7), space = 'left', height = 0.5 )
p.xlab <- list(label = '')
p.ylab <- list(label = '')
p.scales <- list(cex=0.7, x = list( tck = c(1,1)),
                 y = list(labels = rev(c('J','F','M','A','M','J','J','A','S','O','N','D'))))


sccpl <- lapply(sccl, function(x){
  p <- hovmoller(x,
                 interpolate = TRUE,
                 xscale.components = xscale.raster.subticks,
                 par.settings = BuRdTheme,
                 xlab = p.xlab,
                 ylab = p.ylab,
                 scales = p.scales,
                 colorkey = p.ckey,
                 at = p.at
  )
  return(p)
})


tiff(file= "sc climatology hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)

p <- c(sccpl[[4]], sccpl[[3]], sccpl[[2]], sccpl[[1]], 
       layout = c(1,4))
update(p, par.settings = BuRdTheme(),
       ylab.right = list(label = c('139-140 E', '138-139 E', '137-138 E', '136-137 E'),
                         cex = 0.75, rot = 270),
       scales = list(alternating = FALSE, y = list(rot = -90)))

dev.off()




# wind stress curl --------------------------------------------------------
wc <- stack("climatology_windstresscurl_1day_136_141_-36.5_-43.grd")
wc <- wc[[c(7,8,9,10,11,12,1,2,3,4,5,6)]]


# put each location raster stack in a list  
wccl <- list()
for(i in 1:length(el)){
  r <- crop(wc, el[[i]])
  wccl <- c(wccl, r)
}

# setZ for each stack in the list
wccl <- lapply(wccl, function(x) {
  t <- names(x)
  x <- setZ(x, t)
})

# reverse month arrangement and rename to abbr letters
wccl <- lapply(wccl, function(x){
  x <- x[[12:1]]
})

# plot hovmoller diagram
man <- max(unlist(lapply(wccl, maxValue)))
miin <-  min(unlist(lapply(wccl, minValue)))
p.at <- seq(-6e-07, 8e-07, 1e-07)
p.ckey <- list(at = p.at, labels = list(at = p.at, cex = 0.7), space = 'left')
p.xlab <- list(label = '')
p.ylab <- list(label = '')
p.scales <- list(cex=0.7, x = list( tck = c(1,1)),
                 y = list(labels = rev(c('J','F','M','A','M','J','J','A','S','O','N','D'))))


wccpl <- lapply(wccl, function(x){
  p <- hovmoller(x,
                 interpolate = TRUE,
                 xscale.components = xscale.raster.subticks,
                 par.settings = BuRdTheme,
                 xlab = p.xlab,
                 ylab = p.ylab,
                 scales = p.scales,
                 colorkey = p.ckey,
                 at = p.at
  )
  return(p)
})


tiff(file= "wc climatology hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)

p <- c(wccpl[[4]], wccpl[[3]], wccpl[[2]], wccpl[[1]], 
       layout = c(1,4))
update(p, par.settings = BuRdTheme(),
       ylab.right = list(label = c('139-140 E', '138-139 E', '137-138 E', '136-137 E'),
                         cex = 0.75, rot = 270),
       scales = list(alternating = FALSE, y = list(rot = -90)))

dev.off()


# INTER-ANNUAL ------------------------------------------------------------
# ssha --------------------------------------------------------------------
dir(pattern = '.grd')
ssha <- stack("daily_ssha_136_141_-36.5_-43.grd")
tt <- seq(as.Date('1997-01-01'), length = nlayers(ssha), by='day')
ssha <- setZ(ssha, tt)

ssha1 <- crop(ssha, e1)
ssha2 <- crop(ssha, e2)
ssha3 <- crop(ssha, e3)
ssha4 <- crop(ssha, e4)

# plot hovmoller diagram
p.at <- seq(-0.25, 0.35, .05)
p.ckey <- list(labels = list(at = p.at, cex = 0.7), space = 'left', height = 0.5 )
p.xlab <- list(label = 'latitude', cex = 0.75)
p.ylab <- list(label = 'time', cex = 0.75)
p.scales <- list(cex=0.7)



p1 <- hovmoller(ssha1,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "SSHA 136 - 137 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

p2 <- hovmoller(ssha2,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "137 - 138 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey)

p3 <- hovmoller(ssha3,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "138 - 139 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

p4 <- hovmoller(ssha4,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "139 - 140 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

tiff(file= "ssha hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)
grid.arrange(p4, p3, p2, p1, ncol = 1, nrow = 4)
dev.off()


# sst ---------------------------------------------------------------------

sst <- stack("daily_sst_136_141_-36.5_-43.grd")
tt <- seq(as.Date('1997-01-01'), length = nlayers(sst), by='day')
sst <- setZ(sst, tt)

sst1 <- crop(sst, e1)
sst2 <- crop(sst, e2)
sst3 <- crop(sst, e3)
sst4 <- crop(sst, e4)
# plot hovmoller diagram

p.at <- seq(8, 23, 2)
p.ckey <- list(labels = list(at = p.at, space = 'left', height = 0.5 ))
p.xlab <- list(label = 'latitude', cex = 0.75)
p.ylab <- list(label = 'time', cex = 0.75)
p.scales <- list(cex=0.7)


p1 <- hovmoller(sst1,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "sst 136 - 137 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

p2 <- hovmoller(sst2,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "137 - 138 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

p3 <- hovmoller(sst3,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "138 - 139 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

p4 <- hovmoller(sst4,
                interpolate = TRUE,
                yscale.components = yscale.raster.subticks,
                par.settings = BuRdTheme,
                main = list(label = "139 - 140 E",
                            cex = 0.75),
                xlab = p.xlab,
                ylab = p.ylab,
                scales = p.scales,
                colorkey = p.ckey,
                at = p.at)

tiff(file= "sst hovmoller 136-140E 36.5-43S.tiff", width=6, height=11.5, units="in", res=300)
grid.arrange(p4, p3, p2, p1, ncol = 1, nrow = 4)
dev.off()

# CHL ---------------------------------------------------------------------
chl <- stack("weekly_chl_136_141_-36.5_-43.grd") 
tt <- seq(as.Date('1997-09-02'), length = nlayers(chl), by='8 day')
chl <- setZ(chl, tt)

chl1 <- crop(chl, e1)
chl2 <- crop(chl, e2)
chl3 <- crop(chl, e3)
chl4 <- crop(chl, e4)

# # make spatial resolution smaller
# chl_lowres <- aggregate(chl, fact=3)
# chl_lowres <- setZ(chl_lowres, tt)

#hovmoller

p.at <- seq(-6, 6, 1)
p.ckey <- list(at = p.at, labels = list(cex = 0.7), space = 'left', height = 0.5 )
p.xlab <- list(label = 'latitude', cex = 0.75)
p.ylab <- list(label = 'time', cex = 0.75)
p.scales <- list(cex=0.7)



hovmoller(chl1,
          interpolate = TRUE,
          yscale.components = yscale.raster.subticks,
          par.settings = BuRdTheme,
          xscale.components = xscale.raster.subticks,
          main = 'CHL')

hovmoller(chl2,
          interpolate = TRUE,
          yscale.components = yscale.raster.subticks,
          par.settings = BuRdTheme,
          xscale.components = xscale.raster.subticks,
          main = 'CHL')

hovmoller(chl3,
          interpolate = TRUE,
          yscale.components = yscale.raster.subticks,
          par.settings = BuRdTheme,
          xscale.components = xscale.raster.subticks,
          main = 'CHL')

hovmoller(chl4,
          interpolate = TRUE,
          yscale.components = yscale.raster.subticks,
          par.settings = BuRdTheme,
          xscale.components = xscale.raster.subticks,
          main = 'CHL')
