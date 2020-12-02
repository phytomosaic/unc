######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- species *compositions* bootstrap
#     bootstrap species that enter the lichen airscore, keep richness constant
#
#    Rob Smith, phytomosaic@gmail.com, 10 Nov 2020
#
##      GNU General Public License, Version 3.0    ###################
#
# Goal: estimate uncertainty in exceedance values.
#    EXCEEDANCE = OBS - CL
#    OBS = lichen airscore (or CMAQ)
#    CL = from Diversity paper nonlinear equation
#
#
# ####################################################################

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
save(d, file='./res/d_n_scr.rda') # save descriptor matrix for later plotting
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
save(b, file='./res/b_n_scr.rda')
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
######################################################################


######################################################################
### setup for S ratings
o   <- make_ratings_matrix('s')
o   <- match_site_data(o$x, o$spe, 's')
x   <- o$x
d   <- o$d
save(d, file='./res/d_s_scr.rda') # save descriptor matrix for later plotting
pr  <- as.matrix(sweep(o$spe, 1, rowSums(o$spe), '/'))  # species site occ probs
rm(o,d)
dim(x)  # ratings matrix to sample, 8918 sites, 324 species
dim(pr) # sampling probabilities (weights)
`fn`  <- function(x,i) { mean(x[i]) } # bootstrapping function to get the mean
B     <- 999  # number of bootstrap replicates
ncpus <- 11   # number of cores in parallel
### ! ! ! TIMEWARN ! ! !  ~20 min on 11 nodes, B=999, 8918 sites, 324 species
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
save(b, file='./res/b_s_scr.rda')
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
######################################################################




#############################################################################
###     load bootstrap results from file     ###########################
#############################################################################

rm(list=ls())

### load bootstrapped 95% CI for site airscores (each row = one site)
load(file='./res/b_n_scr.rda', verbose=T)
bn  <- t(b) # bootstrap CIs for N
load(file='./res/b_s_scr.rda', verbose=T)
bs  <- t(b) # bootstrap CIs for S
rm(b)

### load matching site-descriptor matrices
load(file='./res/d_n_scr.rda', verbose=T)
dn <- cbind(d,bn)
load(file='./res/d_s_scr.rda', verbose=T)
ds <- cbind(d,bs)
rm(d,bs,bn)

### remove degenerate solutions
i  <- is.na(dn$upr)
sum(i) #  37 degenerate solutions for N
dn <- dn[!i,]
i  <- is.na(ds$upr)
sum(i) # 480 degenerate solutions for S
ds <- ds[!i,]
rm(i)

### calc CI breadth = range of possible values = UNCERTAINTY
dn$ci_rng <- dn$upr - dn$lwr # calc breadth of CIs
ds$ci_rng <- ds$upr - ds$lwr # calc breadth of CIs

### calc site exceedances: bootstrapped mean minus the fixed CL
dn$exc <- dn$mean - 1.5   # N airscores CL = 1.5
ds$exc <- ds$mean - 2.7   # S airscores CL = 2.7

### histograms of uncertainties and exceedances
set_par_mercury(4)
hist(dn$exc, breaks=seq(-5, max(ds$exc)+0.5,by=0.5)) # N fairly tight
hist(ds$exc, breaks=seq(-5, max(ds$exc)+0.5,by=0.5)) # S ranges higher
hist(dn$ci_rng, breaks=seq(0, max(ds$ci_rng)+0.05,by=0.05)) # N fairly tight
hist(ds$ci_rng, breaks=seq(0, max(ds$ci_rng)+0.05,by=0.05)) # S ranges higher

### plot maps of exceedances for presentation
ca <- c('#5E4FA2', '#4F61AA','#4173B3', '#3386BC','#4198B6', '#51ABAE',
        '#62BEA6', '#77C8A4','#8ED1A4', '#A4DAA4','#B8E2A1', '#CBEA9D',
        '#DEF199', '#EAF69F','#F2FAAC', '#FAFDB8', '#FEFAB6', '#FEF0A5',
        '#FEE695', '#FDD985', '#FDC978', '#FDB96A','#FCA75E', '#F99254',
        '#F67D4A', '#F26943','#E85A47', '#DE4B4B', '#D33C4E', '#C1284A',
        '#AF1446', '#9E0142', rep('transparent',2))
`map_it` <- function(x, whichcol = 'exc', pal=ca, ...) {
  plot(x$lon, x$lat, type='n', asp=1.6, ylab='', xlab='', ...)
  maps::map('usa', fill=T, col='white', add=T)
  points(x$lon, x$lat, pch=16, cex=0.5, col=colvec(x[,whichcol], pal=ca))
  maps::map('state', add=T, lwd=0.5, col='#00000050')
  `colorbar` <- function(x, main='', pal, ...) {
    ppar <- par('plt')
    xmin <- ppar[2] - (ppar[2] - ppar[1]) * 0.15
    xmax <- ppar[2] - (ppar[2] - ppar[1]) * 0.05
    ymin <- ppar[3] + (ppar[4] - ppar[3]) * 0.05
    ymax <- ppar[3] + (ppar[4] - ppar[3]) * 0.35
    op <- par(fig=c(xmin,xmax,ymin,ymax), mar=c(0,2,2,0),
              oma=c(0,0,0,0), new=T)
    `bar` <- function(pal, min, max, nticks=5,
                      ticks=round(seq(min, max, len=nticks),0),
                      main=expression(kg~ha^-1~y^-1*'       ')) {
      scl <- (length(pal)-1)/(max-min)
      plot(c(0,10), c(min,max), type='n', bty='n',
           xaxt='n', xlab='', yaxt='n', ylab='', main='',
           font=1)
      title(main, adj = 0.5, line = 0.5, cex.main=0.5)
      axis(2, ticks, las=1, cex.axis=0.4,
           col='transparent', col.ticks='transparent')
      for (i in 1:(length(pal)-1)) {
        y <- (i-1) / scl + min
        rect(0, y, 10, y+1/scl, col=pal[i], border=NA)
      }
    }
    bar(pal=pal, min(x, na.rm=T), round(max(x, na.rm=T), 0))
    par(op)
  }
  colorbar(x[,whichcol], pal=pal)
}

### map exceedances
png('./fig/fig_01_map_n_exceedances.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(dn, 'exc', main='Airscore N exceedance') # ylim=c(30,49)
dev.off()
png('./fig/fig_01_map_s_exceedances.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(ds, 'exc', main='Airscore S exceedance') # ylim=c(30,49)
dev.off()
### map exceedance *uncertainty*
png('./fig/fig_01_map_n_exc_uncertainty.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(dn, 'ci_rng', main='CI breadth for N exceedance')
dev.off()
png('./fig/fig_01_map_s_exc_uncertainty.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(ds, 'ci_rng', main='CI breadth for S exceedance')
dev.off()





### DIAGNOSTICS

# correct NA climate values
dn$ubc_mat[dn$ubc_mat < -99] <- NA
ds$ubc_mat[ds$ubc_mat < -99] <- NA

### diagnostics N ........
png('./fig/fig_00_diagnostics_n.png',
    wid=12, hei=8, uni='in', res=700, bg='transparent')
set_par_mercury(6)
ca <- colvec(dn$cmaq_n_3yroll, alpha=0.9)
plot(jitter(dn$sr, factor=2),  dn$scr_obs, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Obsvd airscore')
plot(dn$cmaq_n_3yroll, dn$scr_obs, col=ca, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Obsvd airscore')
plot(dn$ubc_mat, dn$scr_obs, col=ca,
     xlab='Mean ann temp', ylab='Obsvd airscore')
plot(jitter(dn$sr, factor=2),  dn$ci_rng, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Uncertainty')
plot(dn$cmaq_n_3yroll, dn$ci_rng, col=ca, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Uncertainty')
plot(dn$ubc_mat, dn$ci_rng, col=ca,
     xlab='Mean ann temp', ylab='Uncertainty')
dev.off()
### diagnostics S ........
png('./fig/fig_00_diagnostics_s.png',
    wid=12, hei=8, uni='in', res=700, bg='transparent')
set_par_mercury(6)
ca <- colvec(ds$cmaq_s_3yroll, alpha=0.9)
plot(jitter(ds$sr, factor=2),  ds$scr_obs, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Obsvd airscore')
plot(ds$cmaq_s_3yroll, ds$scr_obs, col=ca, xlim=c(0,20),
     xlab='CMAQ S dep', ylab='Obsvd airscore')
plot(ds$ubc_mat, ds$scr_obs, col=ca,
     xlab='Mean ann temp', ylab='Obsvd airscore')
plot(jitter(ds$sr, factor=2),  ds$ci_rng, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Uncertainty')
plot(ds$cmaq_s_3yroll, ds$ci_rng, col=ca, xlim=c(0,20),
     xlab='CMAQ S dep', ylab='Uncertainty')
plot(ds$ubc_mat, ds$ci_rng, col=ca,
     xlab='Mean ann temp', ylab='Uncertainty')
dev.off()




### mean and CI are positively related to richness..... WHY???
`plot_ci` <- function (z, ...) {
  y   <- z$mean
  lwr <- z$lwr
  upr <- z$upr
  x   <- 1:NROW(z)
  ylm <- c(min(c(y,lwr),na.rm=T) - 0.4, max(c(y,upr),na.rm=T) + 0.1)
  plot(x, y, pch=16, cex=0.2, ylim=ylm, ylab='Airscore \u00B1 95% CI',
       xaxs='i', yaxs='r', xaxt='n', ...)
  segments(x0=x, x1=x, y0=lwr, y1=upr, lwd=0.1, lend='butt')
}

### plot bootstrapped airscores vs richness, then CMAQ
png('./fig/fig_01_CIs_richness-CMAQ-n.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(dn$sr, dn$cmaq_n_3yroll)
plot_ci(dn[o,], xlab='Sites, ordered by richness then CMAQ N dep')
csr <- cumsum(table(dn$sr))[1:27]
text(c(0,csr)+diff(c(0,csr,length(dn$sr)))*0.5,17,labels=c(4:30,'>30'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()
png('./fig/fig_01_CIs_richness-CMAQ-s.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(ds$sr, ds$cmaq_s_3yroll)
plot_ci(ds[o,], xlab='Sites, ordered by richness then CMAQ S dep')
csr <- cumsum(table(ds$sr))[1:27]
text(c(0,csr)+diff(c(0,csr,length(ds$sr)))*0.5,17,labels=c(4:30,'>30'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()


### plot bootstapped airscores vs observed CMAQ
png('./fig/fig_02_CIs_CMAQ-n.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(dn$cmaq_n_3yroll, dn$sr)
plot_ci(dn[o,], xlab='Sites, ordered by CMAQ N dep')
csr <- cumsum(table(round(dn$cmaq_n_3yroll[o])))[1:15]
text(c(0,csr)+diff(c(0,csr,length(dn$sr)))*0.5, 17, labels=c(1:15,'>15'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()
png('./fig/fig_02_CIs_CMAQ-s.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(ds$cmaq_n_3yroll, ds$sr)
plot_ci(ds[o,], xlab='Sites, ordered by CMAQ S dep')
csr <- cumsum(table(round(ds$cmaq_s_3yroll[o])))[1:15]
text(c(0,csr)+diff(c(0,csr,length(ds$sr)))*0.5, 17, labels=c(1:15,'>15'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()


# ### zoom-in
# png('./fig/fig_01_CIs_zoom_SR11.png',
#     wid=4.5, hei=4.5, uni='in', res=700, bg='transparent')
# set_par_mercury(1, pty='s')
# o   <- order(dn$cmaq_n_3yroll)
# i   <- which(dn$sr==11)
# set.seed(121)
# i   <- sort(sample(i,size=50))
# tci <- dn[o,]
# tci <- tci[i,]
# `plot_ci` <- function (z, ...) {
#   y   <- z$mean
#   lwr <- z$lwr
#   upr <- z$upr
#   x   <- 1:NROW(z)
#   ylm <- c(min(c(y,lwr),na.rm=T) - 0.4, max(c(y,upr),na.rm=T) + 0.1)
#   plot(x, y, pch=16, cex=0.5, ylim=ylm, ylab='Airscore \u00B1 95% CI',
#        xaxs='r', yaxs='r', xaxt='n', ...)
#   segments(x0=x, x1=x, y0=lwr, y1=upr, lwd=0.5, lend='butt')
# }
# plot_ci(tci, xlab='Sites, ordered by CMAQ N deposition')
# dev.off()

####    END    ####
