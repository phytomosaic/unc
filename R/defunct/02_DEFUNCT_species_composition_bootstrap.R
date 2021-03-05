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

rm(list=ls())
require(ecole)
require(scales)
# require(boot)


# ######################################################################
# ### bootstraps for N ratings
# # load(file='./res/d_n_scr.rda')   # descriptor matrix
# load(file='./res/x_n_scr.rda')   # ratings matrix
# load(file='./res/spe_n_scr.rda') # species abundance matrix
# ls()
# pr  <- as.matrix(sweep(spe, 1, rowSums(spe), '/')) # species site occ probs
# rm(spe)
# dim(x)  # ratings matrix to sample, 8875 sites, 346 species
# dim(pr) # sampling probabilities (weights)
# `fn`  <- function(x,i) { mean(x[i]) } # bootstrapping function to get the mean
# B     <- 999  # number of bootstrap replicates
# ncpus <- 11   # number of cores in parallel
# ### ! ! ! TIMEWARN ! ! !  ~19.0 min on 11 nodes, B=999, 8875 sites, 346 species
# cat(paste0('start time: ', time_start <- Sys.time()), '\n')
# b <- sapply(1:NROW(x), function(i) { # for each row in x
#   if(i %% floor(NROW(x)/100) == 0) {
#     cat(paste0(round(i/NROW(x),2)*100, '% '))
#   }
#   y    <-  x[i,]                 # subset row of species ratings (traits)
#   wts  <- pr[i,]                 # species occ probabilities per row
#   bb   <- boot::boot(y, fn, R=B,
#                      weights = wts,
#                      parallel='multicore',
#                      ncpus=ncpus)
#   ci  <- boot::boot.ci(bb, type='perc', conf=0.95)$percent[4:5]
#   if(is.null(ci)) { # if zero variance...
#     ci <- c(NA,NA)
#   }
#   out <- c(mean(bb$t), ci)
#   names(out) <- c('mean','lwr','upr')
#   out
# })
# save(b, file='./res/b_n_scr.rda')
# cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
# ######################################################################
#
#
# ######################################################################
# ### bootstraps for S ratings
# # load(file='./res/d_s_scr.rda')   # descriptor matrix
# load(file='./res/x_s_scr.rda')   # ratings matrix
# load(file='./res/spe_s_scr.rda') # species abundance matrix
# ls()
# pr  <- as.matrix(sweep(spe, 1, rowSums(spe), '/')) # species site occ probs
# rm(spe)
# dim(x)  # ratings matrix to sample, 8918 sites, 324 species
# dim(pr) # sampling probabilities (weights)
# `fn`  <- function(x,i) { mean(x[i]) } # bootstrapping function to get the mean
# B     <- 999  # number of bootstrap replicates
# ncpus <- 11   # number of cores in parallel
# ### ! ! ! TIMEWARN ! ! !  ~20 min on 11 nodes, B=999, 8918 sites, 324 species
# cat(paste0('start time: ', time_start <- Sys.time()), '\n')
# b <- sapply(1:NROW(x), function(i) { # for each row in x
#   if(i %% floor(NROW(x)/100) == 0) {
#     cat(paste0(round(i/NROW(x),2)*100, '% '))
#   }
#   y    <-  x[i,]                 # subset row of species ratings (traits)
#   wts  <- pr[i,]                 # species occ probabilities per row
#   bb   <- boot::boot(y, fn, R=B,
#                      weights = wts,
#                      parallel='multicore',
#                      ncpus=ncpus)
#   ci  <- boot::boot.ci(bb, type='perc', conf=0.95)$percent[4:5]
#   if(is.null(ci)) { # if zero variance...
#     ci <- c(NA,NA)
#   }
#   out <- c(mean(bb$t), ci)
#   names(out) <- c('mean','lwr','upr')
#   out
# })
# save(b, file='./res/b_s_scr.rda')
# cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
# ######################################################################




#############################################################################
###     load bootstrap results from file     ###########################
#############################################################################

# rm(list=ls())

# ### load bootstrapped 95% CI for site airscores (each row = one site)
# load(file='./res/b_n_scr.rda')   # bootstrap CIs
# bn  <- t(b)                      # bootstrap CIs for N
# load(file='./res/b_s_scr.rda')   # bootstrap CIs
# bs  <- t(b)                      # bootstrap CIs for S
# load(file='./res/d_n_scr.rda')   # descriptor matrix
# dn <- cbind(d,bn)
# load(file='./res/d_s_scr.rda')   # descriptor matrix
# ds <- cbind(d,bs)
# rm(b,d,bs,bn)




### load bootstrapped 95% CI for site airscores (each row = one site)
load(file='./res/b_n_scr.rda')   # bootstrap CIs
bn  <- b                      # bootstrap CIs for N
load(file='./res/b_s_scr.rda')   # bootstrap CIs
bs  <- b                      # bootstrap CIs for S
load(file='./res/d_n_scr.rda')   # descriptor matrix
dn <- cbind(d,bn)
load(file='./res/d_s_scr.rda')   # descriptor matrix
ds <- cbind(d,bs)
rm(b,d,bs,bn)

head(dn)
head(ds)




# ### remove degenerate solutions
# i  <- is.na(dn$upr)
# sum(i) #  37 degenerate solutions for N
# dn <- dn[!i,]
# i  <- is.na(ds$upr)
# sum(i) # 480 degenerate solutions for S
# ds <- ds[!i,]
# rm(i)

# ### calc CI breadth = range of possible values = UNCERTAINTY
# dn$ci_rng <- dn$upr - dn$lwr # calc breadth of CIs
# ds$ci_rng <- ds$upr - ds$lwr # calc breadth of CIs

### calc site exceedances: bootstrapped mean minus the fixed CL
dn$exc <- dn$med - 1.5   # N airscores CL = 1.5
ds$exc <- ds$med - 2.7   # S airscores CL = 2.7

### histograms of uncertainties and exceedances
set_par_mercury(4)
hist(dn$exc, breaks=seq(-5, max(ds$exc)+0.5,by=0.5)) # N fairly tight
hist(ds$exc, breaks=seq(-5, max(ds$exc)+0.5,by=0.5)) # S ranges higher
hist(dn$ci_rng, breaks=seq(-1, max(ds$ci_rng)+0.05,by=0.5)) # N fairly tight
hist(ds$ci_rng, breaks=seq(-1, max(ds$ci_rng)+0.05,by=0.5)) # S ranges higher


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
  y   <- z$med
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
o <- order(ds$cmaq_s_3yroll, ds$sr)
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
#   y   <- z$med
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
