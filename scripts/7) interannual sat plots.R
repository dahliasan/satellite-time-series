# make interannual surface plots

library(dplyr)
library(ggplot2)
library(viridis)
library(maps)
library(mapdata)
library(marmap)
library(lubridate)
library(igraph)
library(scales)

# load list of dataframes
load('satellite data list (-36).Rdata')

# calculate annual monthly means
# ml <- lapply(ll, FUN = function(x) {
#   x %>% mutate(month = month(date, abbr = TRUE, label = TRUE), year = year(date)) %>% 
#     group_by(x, y, year, month) %>% 
#     filter(year > 1997 & year < 2017) %>% 
#     summarise(mean = mean(v1, na.rm = TRUE), sd = sd(v1, na.rm = TRUE))
# })
# save(ml, file = 'saellite data yearmon list (-36).Rdata')

# Run before plot -------------------------------------------
load('./extracted enviro data/saellite data yearmon list (-36).Rdata')

# log chl values
ml$chl <- ml$chl %>% mutate(mean = log(mean), sd = log(sd))
ml$glob <- ml$glob %>% mutate(mean = log(mean), sd = log(sd))

# merge scv and scu df
ml$scu$vm <- ml$scv$mean
colnames(ml$scu)[5] <- 'um'
ml$scu <- ml$scu %>%  mutate(mag = sqrt(um^2 + vm^2))

# get map 
map <- map_data("worldHires")

# get bathy data
bathy <- getNOAA.bathy(136,141,-43, -36.5,res=1, keep=TRUE)
bathy_map <- tbl_df(fortify(bathy))
b2000 <- bathy_map %>% filter(z >= -2000)
b2000 <- b2000 %>% group_by(x) %>% summarise(y = min(y))


# plot sst, ssha, chl type data -------------------------------------------

basicymonplot <- function(x, vn, st = 12) {
  # x = list of annual monthly dataframes with mean and sd variable
  # vn = vector of character data variable names
  # st = surface temperature for stf 
  # b2000 = dataframe for bathy up till 2000 
  
  out <- list()
  
  # first, for every month
  for(k in 1:12){
    print(k)
    mo <- month.abb[k]
    dl <- lapply(x, function(x){
      x <- x %>% filter(month == mo)
    })
    
    # get stf data
    stf <- ml$sst %>% filter(month == mo, mean >= st & mean < st+1)
    
    # second, plot each variable 
    for(i in 1:length(vn)){
      smin <- min(x[[i]]$mean, na.rm = TRUE)
      smax <- max(x[[i]]$mean, na.rm = TRUE)
      p <- ggplot(data = dl[[i]], aes(x, y)) +
        geom_raster(aes(fill = mean)) +
        geom_contour(aes(z = mean), color = 'white', size = 0.1, alpha = 1) +
        facet_wrap(~ year) +
        scale_fill_viridis(limits = c(smin, smax), name = vn) +
        geom_line(data = b2000, aes(x = x, y = y), color = 'grey', alpha = 0.7) +
        geom_map(map_id = "Australia", map = map) +
        geom_raster(data = stf, aes(x = x, y = y), fill = 'black', alpha = 0.2) +
        ylim(-43, -36) +
        xlim(136, 141) +
        labs(title = sprintf('%s (%s)', vn[i], mo), y = 'Lat', x = 'Lon') + 
        theme(text = element_text(size = 5)) + 
        theme_bw()
      
      # assign(vn[i],p)
      
      png(filename = paste0('./plots/annual plots/', sprintf('%s %s.%s', mo, vn[i], 'png')),  width=8, height=8, units= "in", res = 300)
      print(p)
      dev.off()
      
      p <- list(p)
      names(p) <- sprintf('%s (%s)', vn[i], mo)
      out <- c(out, p)
    }
  }
  return(out)
}
p <- basicymonplot(x = list(ml$ssha), vn = 'SSHA')

# plot surface currents ---------------------------------------------------

arrowymonplot <- function(x, vn, scalar = 4 , st = 12){
  # x = list of annual monthly dataframes with mean, sd, vector components (um, vm)
  # and mag (magnitude)
  # vn = vector of character data variable names
  # scaler = scaler for arrows 
  # st = surface temperature for stf 
  # b2000 = dataframe for bathy up till 2000 
  
  out <- list()
  
  # first, for every month
  for(k in 1:12){
    print(k)
    mo <- month.abb[k]
    dl <- lapply(x, function(x){
      x <- x %>% filter(month == mo)
    })
    
    # get stf data
    stf <- ml$sst %>% filter(month == mo, mean >= st & mean < st+1)
    
    #### plot surface current ####
    for(i in 1:length(vn)){
      p <- ggplot(data = dl[[i]],
                  aes(x = x, y = y)) +
        geom_raster(aes(fill = mag)) +
        geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.2, weight = 2)+
        geom_segment(arrow = arrow(length = unit(0.1, 'cm')), aes(xend = x + um * scaler, yend = y + vm * scaler), lwd = 0.4) +
        facet_wrap(~ year) +
        scale_fill_viridis() +
        geom_map(map_id = "Australia", map = map, colour = "grey") +
        labs(title = sprintf('%s (%s)', vn[i], mo)) +
        geom_tile(data = stf, aes(x = x, y = y), fill = 'black', alpha = 0.3, color = 'grey') +
        ylim(-43, -36) +
        xlim(136, 141)
      
      tiff(filename = sprintf('%s %s.%s', mo, vn[i], 'tiff'),  width=7, height=7, units= "in", res = 300)
      print(p)
      dev.off()
      p <- list(p)
      names(p) <- sprintf('%s (%s)', vn[i], mo)
      out <- c(out, p)
    }
  }
  return(out)
}
p2 <- arrowymonplot(x = list(ml$scu), vn = 'surface current')


# plot curl type data ----------------------------------------------------
curlymonplot <- function(x, vn, st = 12) {
  # x = list of annual monthly dataframes with mean and sd variable 
  # vn = vector of character data variable names
  # st = surface temperature for stf 
  # b2000 = dataframe for bathy up till 2000 
  
  out <- list()
  
  # first, for every month
  for(k in 1:12){
    print(k)
    mo <- month.abb[k]
    dl <- lapply(x, function(x){
      x <- x %>% filter(month == mo)
    })
    
    # get stf data
    stf <- ml$sst %>% filter(month == mo, mean >= st & mean < st+1)
    
    # second, plot each variable 
    for(i in 1:length(vn)){
      p <- ggplot(dl[[i]], aes(x = x, y = y)) +
        geom_tile(aes(fill = mean )) +
        geom_contour(aes(z = sd), color = 'white', size = 0.2) +
        geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.5, weight = 2) +
        facet_wrap(~ year) +
        ylim(-43, -36) +
        xlim(136, 141)+
        geom_map(map_id = "Australia", map = map, colour = "grey") +
        labs(title = sprintf('%s (%s)', vn[i], mo)) 
      
      mi <- (dl[[i]] %>%  filter(mean < 0))$mean
      mibreaks <- c( median(mi),min(mi))
      ma <- (dl[[i]] %>%  filter(mean > 0))$mean
      mabreaks <- c( max(ma),median(ma))
      my_colpal <- igraph::diverging_pal(11)[c(4:6,10:11)]
      brlab <- as.character(signif(c(mabreaks, 0, mibreaks)*10^5, 3))
      scname <- 'mean 10^-5'
      p <- p + scale_fill_gradientn(colours = my_colpal, 
                                       values = scales::rescale(c(mabreaks, 0, mibreaks)),
                                       guide = "legend",
                                       breaks = c(mabreaks, 0, mibreaks),
                                       labels = brlab,
                                       name = scname)
      
      tiff(filename = sprintf('%s %s.%s', mo, vn[i], 'tiff'),  width=7, height=7, units= "in", res = 300)
      print(p)
      dev.off()
      
      p <- list(p)
      names(p) <- sprintf('%s (%s)', vn[i], mo)
      out <- c(out, p)
    }
  }
  return(out)
}

p3 <- curlymonplot(x = list(ml$ascurl, ml$qscurl), vn = c('ascat curl', 'quikscat curl'))



# plot chl x surface currents ---------------------------------------------
x <- list(list(ml$chl, ml$scu))
CHLxSCymonplot <- function(x, vn, scaler = 3 , st = 12){
  # x = list of lists with CHL and SURFACE CURRENT
  # of annual monthly dataframes with mean, sd, vector components (um, vm)
  # and mag (magnitude)
  # vn = vector of character data variable names
  # scaler = scaler for arrows 
  # st = surface temperature for stf 
  # b2000 = dataframe for bathy up till 2000 
  
  out <- list()
  
  # first, for every month
  for(k in 1:12){
    mo <- month.abb[k]
    dl <- lapply(x[[1]], function(x){
      x <- x %>% filter(month == mo)
    })
    
    # get stf data
    stf <- ml$sst %>% filter(month == mo, mean >= st & mean < st+1)
    
    smin <- min(x[[1]][[1]]$mean, na.rm = TRUE)
    smax <- max(x[[1]][[1]]$mean, na.rm = TRUE)
    
    #### plot chl x surface current ####
    p <- ggplot(data = dl[[1]]) + 
      geom_raster(data = dl[[1]], aes(x = x, y = y, fill = mean)) +
      geom_segment(data = dl[[2]], 
                   arrow = arrow(length = unit(0.08, 'cm')), 
                   aes(x = x, y = y, xend = x + um * scaler, yend = y + vm * scaler), 
                   lwd = 0.13, colour = 'grey70') +
      facet_wrap(~ year) +
      scale_fill_viridis(limits = c(smin, smax), name = 'Chl') +
      geom_line(data = b2000, aes(x = x, y = y), color = 'grey', alpha = 0.7) +
      geom_map(map_id = "Australia", map = map) +
      geom_raster(data = stf, aes(x = x, y = y), fill = 'black', alpha = 0.2) +
      labs(title = sprintf('%s (%s)', vn, mo), x = 'Lon', y = 'Lat') +
      ylim(-43, -36) +
      xlim(136, 141) + 
      theme(text = element_text(size = 5)) + 
      theme_bw()
      
      png(filename = paste0('./plots/annual plots/',sprintf('%s %s.%s', mo, vn, 'png')),  width=8, height=8, units= "in", res = 300)
      print(p)
      dev.off()
      p <- list(p)
      names(p) <- sprintf('%s (%s)', vn, mo)
      out <- c(out, p)
      print(k)
  }
  return(out)
}
p4 <- CHLxSCymonplot(x = list(list(ml$chl, ml$scu)), 
                     vn = 'chl and surface current')


# combine wind stress curl and surface current plots ----------------
scaler = 6
p <- ggplot(data = ml$curl) + 
  geom_tile(data = ml$curl, aes(x = x, y = y, fill = mean < 0)) +
  geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.7) +
  geom_segment(data = ml$scu, 
               arrow = arrow(length = unit(0.1, 'cm')), 
               aes(x = x, y = y, alpha = mag, xend = x + um * scaler, yend = y + vm * scaler), 
               lwd = 0.5) +
  facet_wrap(~ month) + 
  ylim(-43, -36) + 
  geom_map(map_id = "Australia", map = map, colour = "grey") + 
  labs(title = 'wind stress curl and surface current') 

png(filename = paste('./plots/wind stress curl and surface current', '.png', sep = ''),  width=7, height=7, units= "in", res = 300)
p
dev.off()

# combine chl and wind stress curl ------------
p <- ggplot(data = ml$chl) + 
  geom_point(data = ml$curl, aes(x = x, y = y, alpha = mean < 0)) +
  geom_tile(data = ml$chl, aes(x = x, y = y, fill = mean), alpha = 0.9) +
  geom_contour(data = b2000, aes(x = x, y = y, z = z), color = 'grey', alpha = 0.7) +
  facet_wrap(~ month) + 
  ylim(-43, -36) + 
  geom_map(map_id = "Australia", map = map, colour = "grey") + 
  labs(title = 'wind stress curl and surface current') + 
  scale_fill_viridis()


# plot 1998/1999 chl and ssha during strong el nino ----------------------------
ssha <- ml$ssha %>% filter(year == 1998)
chl <- ml$chl %>% filter(year == 1998)
pal <- wesanderson::wes_palette('Zissou1')
p1 <- ggplot(data = chl, aes(x, y)) + 
  geom_raster(aes(fill = mean)) +
  scale_fill_viridis(name = 'Chl-a') +
  geom_line(data = b2000, aes(x = x, y = y), color = 'black', alpha = 0.7) +
  # geom_contour(data = ssha, aes(z = mean), size = 0.3, colour = 'white') + 
  stat_contour(data = ssha, aes(z = mean, colour = ..level..), size = 0.3) + 
  geom_map(map_id = "Australia", map = map) +
  labs( x = 'Lon', y = 'Lat') +
  facet_wrap(~ month) +
  ylim(-43, -36) +
  xlim(136, 141) + 
  theme(text = element_text(size = 5)) + 
  theme_bw() + 
  scale_colour_gradient2(high = pal[5], low = pal[1], name = 'SSHA')

png(filename = paste('./plots/annual plots/1998 chl vs ssha.png', '.png', sep = ''),  width=8, height=8, units= "in", res = 300)
p1
dev.off()


curr <- ml$scu %>% filter(year == 1998)
scaler = 3
ggplot(data = chl, aes(x, y)) + 
  geom_raster(aes(fill = mean)) +
  scale_fill_viridis(name = 'Chl-a') +
  geom_line(data = b2000, aes(x = x, y = y), color = 'black', alpha = 0.7) +
  geom_segment(data = curr, 
               arrow = arrow(length = unit(0.08, 'cm')), 
               aes(x = x, y = y, xend = x + um * scaler, yend = y + vm * scaler), 
               lwd = 0.13, colour = 'grey70') +
  geom_map(map_id = "Australia", map = map) +
  labs( x = 'Lon', y = 'Lat') +
  facet_wrap(~ month) +
  ylim(-43, -36) +
  xlim(136, 141) + 
  theme(text = element_text(size = 5)) + 
  theme_bw() + 
  scale_colour_gradient2(high = pal[5], low = pal[1], name = 'SSHA')
