### Plot BFAST objects ### 
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(bfast)
library(zoo)
library(broom)
# library(grid)
# library(gridExtra)

rm(list = ls())
## Load prepared BFAST objects
load('bfast zonal objects 6Sep2017.Rdata')
# load('bfast zonal objects 6Sep2017 h0.21.Rdata')
load('bfast soi 1998-2016.Rdata')
shf98bf <- bfl[[1]]
oc98bf <- bfl[[2]]
sam98bf <- bfl[[3]]
names(sam98bf) <- 'sam'


# ggplot bfast objects (Create Functions) ----------------------------------------------------
# source('plotbfast2.R')
# x = shf98bf[[1]]
BFASTggplot <- function(x, ylab1 = NULL, xlab1 = NULL, title = NULL, plotnoise = FALSE) {
  # ylab1 = NULL;xlab1 = NULL;title = NULL
  niter <- length(x$output)
  out <- x$output[[niter]]
  Trend.bp <- !x$nobp$Vt
  Tt <- data.frame(x = index(out$Tt), y = coredata(out$Tt))
  St <- data.frame(x = index(out$St), y = coredata(out$St))
  noise <- data.frame(x = index(out$Nt), y = coredata(out$Nt))
  Yt <- data.frame(x = index(x$Yt), y = coredata(x$Yt))
  startyr <- substr(first(noise$x),1,4)
  endyr <- substr(last(noise$x),1,4)
  freq <- frequency(out$Tt)
  
  # Set plot parameters
  
  margin <- theme(plot.margin = unit(c(0, 1, 0.1, 0.5), "lines"),
                  panel.grid = element_blank(),
                  text = element_text(size = 8),
                  plot.title = element_text(hjust = 0.5, margin = margin(5,0,1,0,'pt')),
                  axis.title.x = element_text(hjust = 0.5, margin = margin(1,0,1,0,'pt')),
                  panel.background = element_rect(fill = "white", colour = 'black'))
  
  
  # plot raw data
  p1 <- ggplot(data = Yt, aes(x = x, y = y))+
    geom_line(colour = 'grey60') +
    labs(y = ylab1, x = xlab1, title = title) + 
    margin + 
    scale_x_continuous(breaks = noise$x[seq(1,nrow(noise),freq)], labels = substr(seq(startyr, endyr, 1), 3, 4))
  
  # Fit ANOVA on segments
  if(Trend.bp == TRUE){
    ## ANOVA
    ft <- cbind(seasonal = out$St, trend = out$Tt, remainder = out$Nt)
    tsp(ft) <- tsp(x$Yt)
    ft <- list(time.series = ft)
    fit <- x
    niter <- length(fit$output) # nr of iterations
    out <- fit$output[[niter]]  # output of results of the final fitted seasonal and trend models and nr of breakpoints in both.
    out_ANOVA <- array()
    if (out$Vt.bp[1] > 0) {breaks <- length(out$Vt.bp) } else {breaks <- 0}  # number of breaks
    if (breaks > 0) {
      breakdates <- out$Vt.bp # breakdates
      coefs <- coef(out$bp.Vt) # output coefficients per segment
      sl <- coefs[,2] # slopes
    }
    out_breakdates <- data.frame(x1 = 1:(breaks+1), x2 = 1:(breaks+1))
    
    TS_anova <- fit$Yt - out$St   # time series Yt - St for ANOVA test
    dataframe <- data.frame(TIME=c(1:length(fit$Yt)),DATA=TS_anova)
    
    # determine segment startpoint and endpoint, calculate ANOVA
    for (m in 1:(breaks+1)) {
      startpoint <- if(m==1)  1 else breakdates[[m-1]]
      endpoint <- if(m==(breaks+1)) length(fit$Yt) else breakdates[m]-1
      out_breakdates$x1[m] <- startpoint 
      out_breakdates$x2[m] <- endpoint 
      df2 <- dataframe[startpoint:endpoint,]  # subset of dataframe (section)
      model <- lm(DATA~TIME, data=df2)        # linear model
      modelAnova <- anova(model)              # ANOVA
      out_ANOVA[m] <- modelAnova$Pr[1]        # save p-value
      if(breaks==0) {sl <- model$coefficients[2]}  ## JV updated -- this was causing problems !# slope Tt if breaks == 0
    }
    out_breakdates$group <- 1:nrow(out_breakdates)
    ## end ANOVA
    
    
    ## my version of ANOVA 
    # bpindex <- x$output[[1]]$Vt.bp
    # bpindex <- c(1, bpindex)
    # bpindex2 <- data.frame(x1 = 2:length(bpindex), x2 = 2:length(bpindex))
    # for(i in 2:length(bpindex)){
    #   bpindex2$x1[i-1] <- bpindex[i-1] + 1
    #   bpindex2$x2[i-1] <- bpindex[i]
    # }
    # lastrow <- data.frame(x1 = bpindex[length(bpindex)]+1, x2 = nrow(Tt))
    # bpindex2 <- bind_rows(bpindex2, lastrow)
    # bpindex2$x1[1] <- 1
    # bpindex2$group <- 1:nrow(bpindex2)
    # 
    # ## ANOVA
    # TS_anova <- x$Yt - out$St
    # dataframe <- data.frame(x=c(1:length(x$Yt)),y=TS_anova)
    # # mod <- map2(bpindex2$x1, bpindex2$x2, function(x,y) lm(y ~ x, data = dataframe[x:y,])) # linear model
    # # pval <- mod %>% map(., ~summary(.)$coefficients[2,4]) # save p-value
    # mod <- map2(bpindex2$x1, bpindex2$x2, function(x,y) lm(y ~ x, data = dataframe[x:y,])) # linear model
    # modAnova <- mod %>% map(., ~anova(.))
    # pval <- modAnova %>% map(., ~.$Pr[1]) # save p-value
    ## end my version of ANOVA
    
    ### add colour coded trendlines 
    
    Tt$signif <- NA
    Tt$group <- NA
    Tt2 <- map2(split(out_breakdates,1:nrow(out_breakdates)), out_ANOVA, function(x,y){
      a <- x$x1[1]
      b <- x$x2[1]
      out <- Tt[a:b,]
      out$signif <- y <= 0.05
      out$group <- x$group[1]
      return(out)
    }) %>% reduce(bind_rows)
    Tt2$signif <- factor(Tt2$signif, levels = c('TRUE', 'FALSE'))
  
  p3 <- p1 + geom_line(data = Tt2, aes(x = x, y = y, colour = signif, group = group)) + 
    scale_color_manual(values = c('black', 'darkgrey'), guide = FALSE)
  
  ### add confidence intervals for breakdates 
  bp <- data.frame(x = strucchange::breakdates(out$bp.Vt), y = Tt$y[out$bp.Vt$breakpoints])
  ci <- data.frame(x = strucchange::breakdates(out$ci.Vt)[,1], x2 = strucchange::breakdates(out$ci.Vt)[,3])
  ci$group <- 1:nrow(bp)
  p3 <- p3 + 
    # geom_text(data = bp, aes(x = x, y = y, label = round(x, digits = 2)),  size = 2) + 
    geom_vline(xintercept = bp$x, linetype = 'dashed') + 
    geom_segment(data = ci, aes(x = x, y = min(Yt), xend = x2, yend = min(Yt)), arrow = arrow(angle = 90, ends = 'both', length = unit(0.1, 'cm')), colour = 'red')
  
} else{
  p3 <- p1 + geom_line(data = Tt, aes(x = x, y = y), colour = 'black')}

  
  # p4 <- ggplot(data = noise, aes(x = x, y = y))+
  #   geom_line(colour = 'grey60') +
  #   labs(y = 'Noise', x = '') +
  #   scale_x_continuous(breaks = noise$x[seq(1,nrow(noise),freq)], labels = substr(seq(startyr, endyr, 1), 3, 4)) + 
  #   theme_bw() + 
  #   theme(plot.margin = unit(c(0, 1, 0, 0.5), "lines"),
  #                      panel.grid = element_blank(),
  #                      text = element_text(size = 8))
    # theme(plot.margin = unit(c(0,0,0,0.5), "lines"))
  
  # gg <- ggarrange(p3, p4, ncol = 1, nrow = 2, align = 'v')
  
  return(p3)
}



# Plot only BFAST trend ---------------------------------------------------
# plot shelf
p1 <- BFASTggplot(shf98bf$ssha, ylab1 = 'SSHA', title = '(a) Shelf')
p2 <- BFASTggplot(shf98bf$sst, ylab1 = 'SST')
p3 <- BFASTggplot(shf98bf$chl, ylab1 = 'CHL')
p4 <- BFASTggplot(sam98bf[[1]], ylab1 = 'SAM', title = '(c) SAM', xlab1 = 'Year')
pl <- list(p1,p2,p3,p4)
gg <- ggarrange(plotlist = pl, ncol = 1, nrow = 4, heights = c(1.1, 1,1,1.1), align = 'v')
gg
png(filename = paste('bfast ', 'shelf3', '.png', sep = ''),  width=7, height=7, units= "in", res = 300)
print(gg)
dev.off()

# plot oceanic
p1 <- BFASTggplot(oc98bf$ssha, ylab1 = 'SSHA', title = '(b) Oceanic')
p2 <- BFASTggplot(oc98bf$sst, ylab1 = 'SST')
p3 <- BFASTggplot(oc98bf$chl, ylab1 = 'CHL', xlab1 = 'Year')
p4 <- BFASTggplot(soi98bf, ylab1 = 'SOI', title = '(d) SOI', xlab1 = 'Year')
pl <- list(p1,p2,p3,p4)
# gg2 <- ggarrange(plotlist = pl, ncol = 1, nrow = 4, heights = c(0.9, 0.8 ,0.9 ,0.8), align = 'v')
gg2 <- ggarrange(plotlist = pl, ncol = 1, nrow = 4, heights = c(1.1, 1,1,1.1), align = 'v')
gg2
png(filename = paste('bfast ', 'ocean3', '.png', sep = ''),  width=7, height=7, units= "in", res = 300)
print(gg2)
dev.off()

# combine shelf and ocean
gg3 <- ggarrange(gg, gg2, ncol = 2, nrow = 1, align = 'hv')
gg3
png(filename = paste('bfast ', 'shelfocean3', '.png', sep = ''),  width=7, height=3.5, units= "in", res = 300)
print(gg3)
dev.off()



# Get table of slope and p values -----------------------------------------
BFASTanovatable <- function(x, term = NULL){
  niter <- length(x$output)
  out <- x$output[[niter]]
  Trend.bp <- !x$nobp$Vt
  Tt <- data.frame(x = index(out$Tt), y = coredata(out$Tt))
  St <- data.frame(x = index(out$St), y = coredata(out$St))
  noise <- data.frame(x = index(out$Nt), y = coredata(out$Nt))
  Yt <- data.frame(x = index(x$Yt), y = coredata(x$Yt))
  startyr <- substr(first(noise$x),1,4)
  endyr <- substr(last(noise$x),1,4)
  freq <- frequency(out$Tt)
  
  # Fit ANOVA on segments
  if(Trend.bp == TRUE){
    ## ANOVA
    ft <- cbind(seasonal = out$St, trend = out$Tt, remainder = out$Nt)
    tsp(ft) <- tsp(x$Yt)
    ft <- list(time.series = ft)
    fit <- x
    niter <- length(fit$output) # nr of iterations
    out <- fit$output[[niter]]  # output of results of the final fitted seasonal and trend models and nr of breakpoints in both.
    out_ANOVA <- array()
    if (out$Vt.bp[1] > 0) {breaks <- length(out$Vt.bp) } else {breaks <- 0}  # number of breaks
    if (breaks > 0) {
      breakdates <- out$Vt.bp # breakdates
      coefs <- coef(out$bp.Vt) # output coefficients per segment
      sl <- coefs[,2] # slopes
    }
    out_breakdates <- data.frame(x1 = 1:(breaks+1), x2 = 1:(breaks+1))
    # determine segment startpoint and endpoint, calculate ANOVA
    for (m in 1:(breaks+1)) {
      startpoint <- if(m==1)  1 else breakdates[[m-1]]
      endpoint <- if(m==(breaks+1)) length(fit$Yt) else breakdates[m]-1
      out_breakdates$x1[m] <- startpoint 
      out_breakdates$x2[m] <- endpoint 
    }
    
    TS_anova <- fit$Yt - out$St   # time series Yt - St for ANOVA test
    dataframe <- data.frame(TIME=c(1:length(fit$Yt)),DATA=TS_anova)
    
    mod <- map2(out_breakdates$x1, out_breakdates$x2, function(x,y) lm(DATA ~ TIME, data = dataframe[x:y,]))
    tidymod <- mod %>% map(., ~tidy(.)[2,]) %>% reduce(bind_rows)
    tidymod$period <- attributes(sl)$names
    tidymod$slope <- sl
    tidymod$n <- out_breakdates$x2 - out_breakdates$x1 +1
    tidymod$term <- term
    
  }
  return(tidymod)
}

names <- c('SSHA', 'SST', 'CHL','SAM', 'SSHA', 'SST', 'CHL')
atab <- map2(c(shf98bf, sam98bf, oc98bf), names, ~BFASTanovatable(.x, .y)) %>% 
  reduce(bind_rows)
write.csv(atab, file = 'BFAST_ANOVA_STATS.csv')



# # plot bfast with trend and season ----------------------------------------
# # plot shelf
# bl <- c(shf98bf[-4], sam98bf)
# tn <- c('ssha','sst', 'chl', 'sam')
# pl <- list()
# for(i in 1:length(tn)){
#   p <- list(BFASTggplot(bl[[i]], paste(tn[i], 'shelf')))
#   pl <- c(pl, p)
# }
# tiff(filename = paste('bfast ', 'shelf', '.tiff', sep = ''),  width=11.5, height=7, units= "in", res = 300)
# grid.arrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]], nrow = 2, ncol = 2)
# dev.off()
# 
# # plot oceanic
# bl <- c(oc98bf[-4], soi98bf)
# tn <- c('ssha','sst', 'chl', 'soi')
# pl <- list()
# for(i in 1:length(tn)){
#   p <- list(BFASTggplot(bl[[i]], paste(tn[i], 'oceanic')))
#   pl <- c(pl, p)
# }
# tiff(filename = paste('bfast ', 'oceanic', '.tiff', sep = ''),  width=11.5, height=7, units= "in", res = 300)
# grid.arrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]], nrow = 2, ncol = 2)
# dev.off()
# 
# # plot oni 
# tiff(filename = 'oni.tiff',  width=11.5, height=7, units= "in", res = 300)
# ggplot(oni98, aes(date, v1)) + 
#   geom_line() +
#   geom_hline(yintercept = -0.8, linetype = 'dashed') +
#   geom_hline(yintercept = 0.8, linetype = 'dashed') + 
#   scale_x_date(date_breaks = '1 year', date_labels = '%y')+
#   ylab('Ocean Nino Index')
# dev.off()
# 
# 
# 
# # Run all code above up till BFAST analysis. Plot each year seperately -----------------------------------------------
# ploteachyr <- function(d1, d2, vn) {
#   # d1 = shelf dataframe with date and varaible (as 'v1')
#   # d2 = oceanic dataframe 
#   # vn = variable name for labels
#   
#   tmp <- d1 %>%  mutate(week = week(date), year = year(date), 
#                         date2 = as.Date(strftime(date, format = '2000-%m-%d')))
#   tmpb <- d2 %>%  mutate(week = week(date), year = year(date),
#                          date2 = as.Date(strftime(date, format = '2000-%m-%d')))
#   
#   p <- ggplot(tmp) +
#     geom_line(data = tmp, aes(x = date2, y = v1))+
#     geom_line(data = tmpb, aes(x = date2, y = v1), colour = muted('blue'))+
#     facet_wrap(~ year) + 
#     ylab(vn) +
#     xlab('month') +
#     scale_x_date(date_breaks = '1 month', date_labels = '%m') + 
#     geom_vline(xintercept = c(as.numeric(as.Date('2000-04-01')),
#                               as.numeric(as.Date('2000-07-01')),
#                               as.numeric(as.Date('2000-10-01'))),
#                linetype = 'dashed', colour = muted('red'))
#   
#   tiff(filename = paste(vn, '_each_year','.tiff', sep =''),  width=11, height=7, units= "in", res = 300)
#   print(p)
#   dev.off()
#   
#   return(p)
#   
# }
# 
# vn <- c('ssha', 'sst', 'chl')
# pel <- list()
# for(i in seq_along(vn)){
#   p <- list(ploteachyr(shf98[[i]], oc98[[i]], vn[i]))
#   pel <- c(pel, p)
# }
# 
# # Plot differences between shelf and oceanic values for SST, SSHA, CHL --------
# ## Run all code above up till BFAST Analysis ---
# 
# ## calculate difference between shelf and oceanic values, sod = shelf ocean diff
# sod <- Map(function(x, y){
#   tbl_df(data.frame(date = x$date, v1 = x$v1 - y$v1))
# }, x = shf98[1:3], y = oc98[1:3])
# 
# sod <- lapply(sodiffs, function(x){
#   x$date2 <- as.Date(strftime(x$date, format = "2000-%m-%d"))
#   x$year <- year(x$date)
#   return(x)})
# 
# # subset ssha and chl only and bind into single df
# sod$ssha$var <- 'ssha'
# sod$chl$var <- 'chl'
# sodg <- rbind(sod$ssha, sod$chl)
# 
# # if you wanted to include sst
# sod$sst$var <- 'sst'
# sodg <- rbind(sod$ssha, sod$chl, sod$sst)
# 
# 
# tiff(filename = 'shelf oceans differences sst ssha sst interannual.tiff', width=7, height=7, units= "in", res = 300)
# ggplot() +
#   geom_line(data = sodg, aes(x = date2, y = v1, color = var))+
#   facet_wrap(~ year) + 
#   ylab('difference') +
#   xlab('month') +
#   scale_x_date(date_breaks = '1 month', date_labels = c('J','J', 'F', 'M', 'A','M','J','J','A','S','O','N','D')) + 
#   geom_vline(xintercept = c(as.numeric(as.Date('2000-03-01')),
#                             as.numeric(as.Date('2000-06-01')),
#                             as.numeric(as.Date('2000-09-01')),
#                             as.numeric(as.Date('2000-12-01'))),
#              linetype = 'dashed', size = 0.3) + 
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 0.3) +
#   theme(text = element_text(size = 8), legend.position = 'bottom', 
#         axis.title.x = element_blank())
# 
# dev.off()
# 
# 
# 
# 
# 

  

# OLD code ----------------------------------------------------------------

  # BFASTggplotNOSEASON <- function(x, title = NULL) {
  #   niter <- length(x$output)
  #   out <- x$output[[niter]]
  #   Trend.bp <- !x$nobp$Vt
  #   Tt <- data.frame(x = index(out$Tt), y = coredata(out$Tt))
  #   St <- data.frame(x = index(out$St), y = coredata(out$St))
  #   noise <- data.frame(x = index(out$Nt), y = coredata(out$Nt))
  #   Yt <- data.frame(x = index(x$Yt), y = coredata(x$Yt))
  #   startyr <- substr(first(noise$x),1,4)
  #   endyr <- substr(last(noise$x),1,4)
  #   freq <- frequency(out$Tt)
  #   
  #   margin = theme(axis.text.x = element_blank(),
  #                  plot.margin = unit(c(0,0,-1,0.5), "lines"), 
  #                  axis.text = element_text(size=6),
  #                  title = element_text(size = 7))
  #   
  #   p1 <- ggplot(data = Yt, aes(x = x, y = y))+
  #     geom_line() +
  #     labs(y = 'data', x = '')+ margin + 
  #     scale_x_continuous(breaks = noise$x[seq(1,nrow(noise),freq)], labels = NULL)
  #   
  #   p3 <- ggplot(data = Tt, aes(x = x, y = y))+
  #     geom_line() +
  #     labs(y = 'trend', x = '') + 
  #     scale_x_continuous(breaks = noise$x[seq(1,nrow(noise),freq)], labels = NULL)
  #   
  #   if(Trend.bp == TRUE){
  #     bp <- data.frame(x = strucchange::breakdates(out$bp.Vt), y = Tt$y[out$bp.Vt$breakpoints])
  #     p3 <- p3 + 
  #       geom_text(data = bp, aes(x = x, y = y, label = round(x, digits = 2)),  size = 2) + 
  #       geom_vline(xintercept = bp$x, linetype = 'dashed') + margin
  #   } else( p3 <- p3 + margin)
  #   
  #   p4 <- ggplot(data = noise, aes(x = x, y = y))+
  #     geom_line() +
  #     labs(y = 'noise', x = '') +
  #     scale_x_continuous(breaks = noise$x[seq(1,nrow(noise),freq)], labels = substr(seq(startyr, endyr, 1), 3, 4)) + 
  #     theme(plot.margin = unit(c(0,0,0,0.5), "lines"), 
  #           axis.text = element_text(size = 6),
  #           title = element_text(size = 7))
  #   
  #   # Left justify plots
  #   # Source: http://stackoverflow.com/a/13295880/496488
  #   pl <- list(p1,p3,p4)
  #   gl <- lapply(pl, function(x){
  #     ggplotGrob(x)
  #   })
  #   
  #   maxWidth =  unit.pmax (gl[[1]]$widths[2:5], gl[[2]]$widths[2:5], gl[[3]]$widths[2:5])
  #   
  #   gl <- lapply(gl, function(x) {
  #     x$widths[2:5] <- as.list(maxWidth)
  #     return(x)
  #   })
  #   
  #   gg <- grid.arrange(gl[[1]],gl[[2]],gl[[3]], ncol=1, 
  #                      top = title)
  #   return(gg)
  # }
  