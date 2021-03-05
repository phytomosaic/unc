######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- species *richness* bootstrap
#     bootstrap species richness, keep composition constant (deletion but not replacement)
#
#    Rob Smith, phytomosaic@gmail.com, 10 Nov 2020
#
##      GNU General Public License, Version 3.0    ###################

require(ecole)
require(labdsv)
require(scales)
require(boot)

### load and clean LICHEN SPECIES occurrence data
rm(list=ls())
`make_ratings_matrix` <- function(which_element='n', ...) {
  s  <- read.csv('./data_raw/MegaDbLICHEN_2020.05.09.csv',
                 stringsAsFactors=F, fileEncoding = 'latin1')
  names(s) <- ecole::clean_text(names(s), lower=T)
  s  <- s[s$growthform == 'Epiphytic macrolichen',] # keep only macrolichens
  s  <- s[,c('megadbid','sci_22chklst','fia_abun',
             'n_depmaxfreq','s_depmaxfreq')]
  s$sci_22chklst <- ecole::clean_text(s$sci_22chklst, TRUE)
  # get species ratings ('peak detection frequency' for N or S)
  #   this comes from HTR's 2018 spline regressions
  wa <- s[!duplicated(s$sci_22chklst),
          c('sci_22chklst','n_depmaxfreq','s_depmaxfreq')]
  wa <- wa[order(wa$sci_22chklst),]
  # reshape wide for abundance matrix
  spe <- labdsv::matrify(s[,c('megadbid','sci_22chklst','fia_abun')])
  # convert spe to trait (ratings to be sampled)
  x <- spe                                   # duplicate
  x[x > 0 ] <- 1                             # convert to binary
  # select which element
  if (which_element == 'n') {
    x <- sweep(x, 2, wa$n_depmaxfreq, '*')   # nitrogen ratings
  } else if (which_element == 's') {
    x <- sweep(x, 2, wa$s_depmaxfreq, '*')   # sulfur ratings
  }
  # filter out weak species
  x[is.na(x)] <- 0L    # zero out any NA (species that lacked ratings)
  j   <- colSums(x)==0 # 287 spp lacking valid spp ratings
  spe <- spe[,!j]
  x   <-   x[,!j]
  i   <- apply(x, 1, function(x) sum(x>0)) < 4 # species-poor (<4 rated species)
  spe <- spe[!i,]
  x   <-   x[!i,]
  x   <- as.matrix(x)
  return(list(spe = spe, x = x))
}
`match_site_data` <- function(x, spe, which_element='n', ...) {
  d <- read.csv('./data_raw/MegaDbPLOT_2020.05.09.csv', stringsAsFactors=F)
  names(d) <- ecole::clean_text(names(d), lower=T)
  d <- d[!(duplicated(d$megadbid)),]      # (!) one megadbid was duplicated
  d$megadbid  <- as.character(d$megadbid)
  rownames(d) <- d$megadbid
  d[,c('latusedd','longusedd')] <- NULL   # kill bad columns
  names(d)[names(d) %in% c('latusenad83','longusenad83')] <- c('lat','lon')
  d <- d[,c('lon','lat',names(d)[!names(d) %in% c('lon','lat')])] # reorder cols
  d <- d[!(is.na(d$lon) | is.na(d$lat)),] # keep only non-NA coords
  d <- d[!(d$lat > 49.5 & d$lat < 50),]   # rm one invalid location in Canada
  d <- d[!(d$lon > (-75) & d$lat < 38),]  # rm one invalid location in Atlantic
  is_ak <- d$lat > 50                     # assign Alaska data where none exist
  d$cmaq_n_3yroll[is_ak] <- d$n_lich_kghay[is_ak]
  d$cmaq_n_3yroll[is_ak & is.na(d$cmaq_n_3yroll)] <- 1.5
  d$cmaq_s_3yroll[is_ak] <- 1.5
  d$ubc_mwmt[d$ubc_mwmt < (-999)] <- NA
  d$meanmaxaugt5y_c[is_ak] <- d$ubc_mwmt[is_ak]
  d$ubc_map[d$ubc_map < 0] <- NA
  d$meanppt5y_cm[is_ak]    <- d$ubc_map[is_ak]
  inm <- intersect(dimnames(x)[[1]], d$megadbid)
  spe <- spe[dimnames(spe)[[1]] %in% inm, ]
  x   <- x[dimnames(x)[[1]] %in% inm, ]
  d   <- d[d$megadbid %in% inm, ]
  spe <- spe[order(dimnames(spe)[[1]]), ]
  x   <- x[order(dimnames(x)[[1]]),]
  d   <- d[order(d$megadbid),]
  # select the element
  if (which_element == 'n') {
    i   <- which(!is.na(d$cmaq_n_3yroll)) # omit if lacking CMAQ values
  } else if (which_element == 's') {
    i   <- which(!is.na(d$cmaq_s_3yroll)) # omit if lacking CMAQ values
  }
  spe <- spe[i, ] # omit if lacking CMAQ values
  x   <- x[i,]    # omit if lacking CMAQ values
  d   <- d[i,]    # omit if lacking CMAQ values
  # check that rownames and dimensions all equal
  stopifnot({
    identical(dimnames(spe)[[1]], d$megadbid) &
      identical(dimnames(x)[[1]], d$megadbid) &
      all.equal(dim(spe), dim(x)) &
      all.equal(NROW(spe), NROW(d))
  })
  # observed site values: airscore and richness
  d$scr_obs <- apply(x, 1, function(i) mean(i[i>0])) # obs airscore (this way)
  d$sr      <- apply(x, 1, function(i) sum(i>0))     # obs richness of RATED spp
  return(list(spe = spe, x = x, d = d))
}

######################################################################
### setup for N ratings
o   <- make_ratings_matrix('n')
o   <- match_site_data(o$x, o$spe, 'n')
x   <- o$x
d   <- o$d
# save(x, file='./res/x_n_scr.rda') # save ratings matrix for later bootstrapping
# save(d, file='./res/d_n_scr.rda') # save descriptor matrix for later plotting
pr  <- as.matrix(sweep(o$spe, 1, rowSums(o$spe), '/')) # species site occ probs
rm(o,d)
dim(x)  # ratings matrix to sample, 8875 sites, 346 species
dim(pr) # sampling probabilities (weights)
`fn`  <- function(x,i) { mean(x[i]) } # bootstrapping function to get the mean
B     <- 999  # number of bootstrap replicates
ncpus <- 11   # number of cores in parallel
### ! ! ! TIMEWARN ! ! !  ~19.0 min on 11 nodes, B=999, 8875 sites, 346 species
cat(paste0('start time: ', time_start <- Sys.time()), '\n')
b <- sapply(1:NROW(x), function(i) { # for each row in x
  if(i %% floor(NROW(x)/100) == 0) {
    cat(paste0(round(i/NROW(x),2)*100, '% '))
  }
  y    <-  x[i,]                 # subset row of species ratings (traits)
  wts  <- pr[i,]                 # species occ probabilities per row
  bb   <- boot::boot(y, fn, R=B,
                     weights = wts,
                     parallel='multicore', ncpus=ncpus)
  ci  <- boot::boot.ci(bb, type='perc', conf=0.95)$percent[4:5]
  if(is.null(ci)) { # if zero variance...
    ci <- c(NA,NA)
  }
  out <- c(mean(bb$t), ci)
  names(out) <- c('mean','lwr','upr')
  out
})
# save(b, file='./res/b_n_sr.rda')
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
######################################################################








######################################################################
### for sampling sufficiency: estimate SEM across varying sample sizes
### sample positive values (presences) from each row of x,
###   with probability equal to species abundances in spe
###   airscore = mean of these randomly sampled trait values
###   (inherently community-weighted by probabilities)
`bootscr` <- function(x, spe, i, B=999, n) {
  if(i %% floor(NROW(x)/100) == 0) {
    cat(paste0(round(i/NROW(x),2)*100, '% '))
  }
  x   <- x[i, ]                 # site row of species optima values
  pr  <- spe[i,] / sum(spe[i,]) # species probabilities per row
  if (missing(n)) {
    n  <- sum(pr > 0)    # sample size = observed richness
  }
  `bmean` <- function(i) { # mean of optima = airscore
    replicate(B, mean(sample(x, size=i, replace=T, prob=pr)))
  }
  res  <- sapply(n, bmean)  ### <--- TO FIX
  res # the bootstrapped airscores across varying sample sizes
}
# ### ! ! ! TIMEWARN ! ! ! 52 sec for 88 plots
# cat(paste0('start time: ', time_start <- Sys.time()), '\n')
# ### for varying sample sizes
# for (N in 2:20) {
#   cat(paste0('\n\t\t< < < < Sample size ',N,' of 20 > > > > \n'))
#   b <- sapply(1:NROW(spe),
#               FUN = bootscr,
#               x   = x,
#               spe = spe,
#               n   = N,
#               simplify='array')
#   dimnames(b)[[3]] <- dimnames(spe)[[1]]
#   save(b, file=paste0('./res/rich_boot/b_',N,'.rda'))
#   rm(b)
# }
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
999 * 8927 * 19 # size of FULL array, if using all plots...




####################3



### load all bootstrap results
### './res/old_richness_boots/'
### './res/rich_boot/'
fnm <- list.files('./res/rich_boot/', pattern='b_[0123456789]')
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) }) # 30 sec......
### calc for each list item
###     - bootstrapped SEM = SD of 999 means
###     - bootstrapped CI = 2.5--97.5 quantiles of 999 means
`f` <- function(b) {
  t(apply(b, 2, function(x) {
    c(sem = sd(x), ci_rng = unname(diff(quantile(x,c(0.025,0.975)))))
  }))
}
b <- lapply(lst, f) # 15 sec......
rm(lst)
b <- do.call(cbind, b)
b_sem <- t(b[,seq(1,NCOL(b),2)])
# b_ci  <- t(b[,seq(2,NCOL(b),2)])


hist(b[,'sem'])

### calculate sufficient species
mqo    <- 1.5 # kg ha y
load(file='./res/d_n_scr.rda', verbose=T)
sr     <- d$sr
# sr     <- apply(spe, 1, function(x) sum(x>0))
n_suff <- apply(b_sem, 2, function(x) min(which(x < mqo), na.rm=T)+1)
sr_deficit <- sr - n_suff # species needed meet MQO

### SEM across varying sample sizes (and MQO)
png('./fig/fig_00_SEM_curves.png',
    wid=7.5, hei=4, uni='in', res=700, bg='transparent')
set_par_mercury(2, mgp=c(1.8,0.2,0))
matplot(b_sem, type='l', xlab='Plot richness',
        ylab=expression('Standard error ('*kg~ha^-1~y^-1*')'),
        lty=1, col='#00000005')
add_text(0.98, 0.37, 'MQO', pos=2, col=2, cex=0.8)
abline(h=1.5, col=2, lwd=2)
hist(n_suff, breaks=1:20, xlab='Minimum richness at MQO',
     ylab='N plots', col='grey', main='', xaxs='i', yaxs='i')
box(bty='l')
dev.off()
rm(b_sem)
### sample size at which bootstrapped SEMs are below MQO threshold
###   (say, 1.5 kg ha-1 y-1)
png('./fig/fig_01_min_sufficient.png',
    wid=7.5, hei=4, uni='in', res=700, bg='transparent')
set_par_mercury(2, mgp=c(2.4,0.2,0), mar=c(4, 4, 0.5, 0.5))
`h` <- function(x,  ...) {
  hist(x, main='', xaxs='i', yaxs='i', ...)
  box(bty='l')
}
sr_deficit[sr_deficit > 60] <- NA # hack to soften extreme values
brk <- seq(-30,60,by=5)
mid <- which(brk==0)
n   <- length(brk)
col <- colorRampPalette(c('red','grey90','darkblue'))(n+mid)[-c(1:mid)]
stopifnot(length(brk) == length(col))
h(sr_deficit, col = col, breaks = brk,
  xlab='Richness sufficiency', ylab='N plots')
abline(v=0, col=2, lwd=2, lty=2)
add_text(0.02,0.95,'Deficit', pos=4)
add_text(0.98,0.95,'Surplus', pos=2)
### MAPS
`co` <- function (x, ...) {
  brk <- seq(-30,60,by=5)
  mid <- which(brk==0)
  n   <- length(brk)
  col <- colorRampPalette(c('red','grey90','darkblue'))(n+mid)[-c(1:mid)]
  col[as.numeric(cut(x, brk, include.lowest=T))]
}
# plot(d$lon, d$lat, col=co(sr_deficit), asp=1.6,
# main='', xlab='Richness sufficiency', ylab='')
`map_it` <- function(x, ...) {
  plot(d$lon, d$lat, type='n', asp=1.6, ylim=c(30,49),
       ylab='', xlab='Richness sufficiency', ...)
  maps::map('usa', add=T)
  points(d$lon, d$lat, pch=16, cex=0.5, col=co(sr_deficit))
  maps::map('state', add=T, lwd=0.5, col='#00000050')
}
map_it()
lbrk <- brk
lbrk[c(rep(T,3),F)] <- ''
legend('bottomright', legend=rev(lbrk), fill=rev(col), cex=0.5,
       y.intersp=0.3,
       bty='n')
dev.off()



#####################################################################
png('./fig/fig_01_map_sufficiency.png',
    wid=6.5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m', mar=c(1,1,0,0))
`h` <- function(x,  ...) {
  hist(x, main='', xaxs='i', yaxs='i', ...) ; box(bty='l')
}
sr_deficit[sr_deficit > 60] <- NA # hack to soften extreme values
brk <- seq(-30,60,by=5)
mid <- which(brk==0)
n   <- length(brk)
col <- colorRampPalette(c('red','grey90','darkblue'))(n+mid)[-c(1:mid)]
stopifnot(length(brk) == length(col))
### MAPS
`co` <- function (x, ...) {
  brk <- seq(-30,60,by=5)
  mid <- which(brk==0)
  n   <- length(brk)
  col <- colorRampPalette(c('red','grey90','darkblue'))(n+mid)[-c(1:mid)]
  col[as.numeric(cut(x, brk, include.lowest=T))]
}
plot(d$lon, d$lat, type='n', asp=1.6, ylim=c(30,49), xlim=c(-125,-50),
     ylab='', xlab='')
maps::map('usa', add=T)
points(d$lon, d$lat, pch=16, cex=0.5, col=co(sr_deficit))
maps::map('state', add=T, lwd=0.5, col='#00000050')
### INSET layout
ppar <- par('plt')
xmin <- ppar[2] - (ppar[2] - ppar[1]) * 0.32
xmax <- ppar[2] - (ppar[2] - ppar[1]) * 0.0001
ymin <- ppar[3] + (ppar[4] - ppar[3]) * 0.05
ymax <- ppar[3] + (ppar[4] - ppar[3]) * 0.5
op <- par(fig=c(xmin,xmax,ymin,ymax), mar=c(1.6,0,0,0),
          oma=c(0,0,0,0), mgp=c(0.5,-0.2,0), tcl=-0.05, new=T)
### HIST
h(sr_deficit, col = col, breaks = brk, yaxt='n',
  xlab='Richness sufficiency', ylab='', cex.lab=0.7, cex.axis=0.5)
abline(v=0, col=2, lwd=2, lty=2)
add_text(-0.02,0.93,'Deficit', pos=4, cex=0.7)
add_text(0.98,0.93,'Surplus', pos=2, cex=0.7)
dev.off()



#####################################################################
png('./fig/fig_01_map_sufficiency_rainbow.png',
    wid=6.5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m', mar=c(1,1,0,0))
`h` <- function(x,  ...) {
  hist(x, main='', xaxs='i', yaxs='i', ...) ; box(bty='l')
}
sr_deficit[sr_deficit > 60] <- NA # hack to soften extreme values
brk <- seq(-30,60,by=5)
mid <- which(brk==0)
n   <- length(brk)
# col <- colorRampPalette(c('red','grey90','darkblue'))(n+mid)[-c(1:mid)]
# stopifnot(length(brk) == length(col))
spectral <- rev(c(
  '#5E4FA2',
  '#4F61AA',
  '#4173B3',
  '#3386BC',
  '#4198B6',
  '#51ABAE',
  '#62BEA6',
  '#77C8A4','#8ED1A4', '#A4DAA4','#B8E2A1', '#CBEA9D',
  '#DEF199',
  '#EAF69F',
  '#FAFDB8',
  '#FEF0A5',
  '#FEE695',
  '#FDC978',
  '#FDB96A',
  '#FCA75E',
  '#F99254',
  '#F67D4A',
  '#F26943',
  '#D33C4E',
  '#C1284A',
  '#AF1446'))



col <- spectral[-c(1:mid+5)]
stopifnot(length(brk) == length(col))
### MAPS
`co` <- function (x, ...) {
  col[as.numeric(cut(x, brk, include.lowest=T))]
}
plot(d$lon, d$lat, type='n', asp=1.6, ylim=c(30,49), xlim=c(-125,-50),
     ylab='', xlab='')
maps::map('usa', add=T)
points(d$lon, d$lat, pch=16, cex=0.5, col=co(sr_deficit))
maps::map('state', add=T, lwd=0.5, col='#00000050')
### INSET layout
ppar <- par('plt')
xmin <- ppar[2] - (ppar[2] - ppar[1]) * 0.32
xmax <- ppar[2] - (ppar[2] - ppar[1]) * 0.0001
ymin <- ppar[3] + (ppar[4] - ppar[3]) * 0.05
ymax <- ppar[3] + (ppar[4] - ppar[3]) * 0.5
op <- par(fig=c(xmin,xmax,ymin,ymax), mar=c(1.6,0,0,0),
          oma=c(0,0,0,0), mgp=c(0.5,-0.2,0), tcl=-0.05, new=T)
### HIST
h(sr_deficit, col = col, breaks = brk, yaxt='n',
  xlab='Richness sufficiency', ylab='', cex.lab=0.7, cex.axis=0.5)
abline(v=0, col=2, lwd=2, lty=2)
add_text(-0.02,0.93,'Deficit', pos=4, cex=0.7)
add_text(0.98,0.93,'Surplus', pos=2, cex=0.7)
dev.off()

####    END    ####
