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

### load and clean LICHEN SPECIES occurrence data
rm(list=ls())
s  <- read.csv('./data_raw/MegaDbLICHEN_2020.05.09.csv',
               stringsAsFactors=F, fileEncoding = 'latin1')
names(s) <- clean_text(names(s), lower=T)
s  <- s[s$growthform == 'Epiphytic macrolichen',] # keep only macros
s  <- s[,c('megadbid','sci_22chklst','fia_abun',
           'n_depmaxfreq','s_depmaxfreq')]
s$sci_22chklst <- clean_text(s$sci_22chklst, TRUE)

### get species ratings ('peak detection frequency' for N or S)
###   this comes from HTR's 2018 spline regressions
wa <- s[!duplicated(s$sci_22chklst),
        c('sci_22chklst','n_depmaxfreq','s_depmaxfreq')]
wa <- wa[order(wa$sci_22chklst),]

### reshape wide for abundance matrix
spe <- labdsv::matrify(s[,c('megadbid','sci_22chklst','fia_abun')])

### convert spe to trait (ratings to be sampled)
x <- spe                                   # duplicate
x[x > 0 ] <- 1                             # convert to binary
x <- sweep(x, 2, wa$n_depmaxfreq, '*')     # nitrogen ratings
# x <- sweep(x, 2, wa$s_depmaxfreq, '*')   # sulfur ratings
x[is.na(x)] <- 0L    # zero out any NA (species that lacked ratings)
j   <- colSums(x)==0 # 287 spp lacking valid spp ratings
spe <- spe[,!j]
x   <-   x[,!j]
i   <- apply(x, 1, function(x) sum(x>0)) < 4 # species-poor (< 4 rated species)
spe <- spe[!i,]
x   <-   x[!i,]
x   <- as.matrix(x)

### get descriptive data (CMAQ and original airscores)
d <- read.csv('./data_raw/MegaDbPLOT_2020.05.09.csv', stringsAsFactors=F)
names(d) <- ecole::clean_text(names(d), lower=T)
d <- d[!(duplicated(d$megadbid)),]      # (!) one megadbid was duplicated
d$megadbid  <- as.character(d$megadbid)
rownames(d) <- d$megadbid
d[,c('latusedd','longusedd')] <- NULL   # kill bad columns
names(d)[names(d) %in% c('latusenad83','longusenad83')] <- c('lat','lon')
d <- d[,c('lon','lat', names(d)[!names(d) %in% c('lon','lat')])] # reorder cols
d <- d[!(is.na(d$lon) | is.na(d$lat)),] # keep only non-NA coords
d <- d[!(d$lat>49.5 & d$lat < 50),]     # rm one invalid location in Canada
d <- d[!(d$lon > (-75) & d$lat < 38),]  # rm one invalid location in Atlantic
is_ak <- d$lat > 50    # assign Alaska data where none exist
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
i   <- which(!is.na(d$cmaq_n_3yroll)) # omit if lacking CMAQ values
spe <- spe[i, ] # omit if lacking CMAQ values
x   <- x[i,]    # omit if lacking CMAQ values
d   <- d[i,]    # omit if lacking CMAQ values

### checks
identical(dimnames(spe)[[1]], d$megadbid) # expect TRUE
identical(dimnames(x)[[1]], d$megadbid)   # expect TRUE
dim(d)   # 8875 sites, descriptor matrix
dim(spe) # 8875 sites, species abundance matrix
dim(x)   # 8875 sites, species ratings 'traits' matrix
d$scr_obs <- apply(x, 1, function(i) mean(i[i>0])) # obs airscore (this way)
d$sr      <- apply(x, 1, function(i) sum(i>0))     # obs richness of RATED spp
rm(i, j, s, wa, inm, is_ak) # cleanup





# ######################################################################
#    #### unwrap me at your peril...   ###################
# ######################################################################
# ### for CLAD WG-2: estimate bootstrap 95% CI of the estimated airscore
# ###   sample positive values (presences) from each row of x,
# ###      with probability equal to species abundances in spe;
# ###   airscore = mean of these randomly sampled trait values
# ###      (inherently community-weighted by probabilities)
# `boot_airscore` <- function(x, i, B=999) {
#     x    <- x[i, ]                 # subset row of species ratings (traits)
#     pr   <- spe[i,] / sum(spe[i,]) # species occ probabilities per row
#     n    <- sum(pr > 0)            # sample size = observed richness
#     res  <- replicate(B, mean(sample(x, size=n, replace=T, prob=pr)))
#     res                            # the bootstrapped airscores
# }
# ### ! ! ! TIMEWARN ! ! ! 10 min for 8927 plots (346 spp)
# cat(paste0('start time: ', time_start <- Sys.time()), '\n')
# b_scr <- sapply(1:NROW(spe), FUN=boot_airscore, x=x, simplify='array')
# dimnames(b_scr)[[2]] <- dimnames(spe)[[1]]
# cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
# save(b_scr, file='./res/b_scr.rda')
# ######################################################################
# ######################################################################
# ######################################################################






### load bootstraps and calc bootstrapped 95% CI for site scores
load(file='./res/b_scr.rda', verbose=T)
ci <- t(apply(b_scr, 2, function(x) quantile(x, c(0.025,0.50,0.975))))
# ci <- ci[!isna,]
ci_rng <- ci[,3] - ci[,1] # calc breadth of CIs
exc <- ci[,2] - 3.5       # calc species richness N exceedance (CL = 3.5)
exc <- exc - 2            # downweight fuzz...

### plot maps of exceedances for presentation
ca <- c('#5E4FA2', '#4F61AA','#4173B3', '#3386BC','#4198B6', '#51ABAE',
        '#62BEA6', '#77C8A4','#8ED1A4', '#A4DAA4','#B8E2A1', '#CBEA9D',
        '#DEF199', '#EAF69F','#F2FAAC', '#FAFDB8', '#FEFAB6', '#FEF0A5',
        '#FEE695', '#FDD985', '#FDC978', '#FDB96A','#FCA75E', '#F99254',
        '#F67D4A', '#F26943','#E85A47', '#DE4B4B', '#D33C4E', '#C1284A',
        '#AF1446', '#9E0142', rep('transparent',2))
`map_it` <- function(x, pal=ca, ...) {
  plot(d$lon, d$lat, type='n', asp=1.6, ylim=c(30,49),
       ylab='', xlab='', ...)
  maps::map('usa', fill=T, col='white', add=T)
  points(d$lon, d$lat, pch=16, cex=0.5, col=colvec(x, pal=ca))
  maps::map('state', add=T, lwd=0.5, col='#00000050')
  `colorbar` <- function(x, main='', pal, ...) {
    ppar <- par('plt')
    xmin <- ppar[2] - (ppar[2] - ppar[1]) * 0.19
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
  colorbar(x, pal=pal)
}
### exceedance
png('./fig/fig_01_map_exceedances.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(x=exc, main='Airscore N exceedance')
dev.off()
### exceedance uncertainty
png('./fig/fig_01_map_exceedance_unc.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(x=ci_rng, main='CI breadth for N exceedance')
dev.off()


### DIAGNOSTICS

### diagnostics........
png('./fig/fig_00_diagnostics.png',
    wid=12, hei=8, uni='in', res=700, bg='transparent')
set_par_mercury(6)
ca <- colvec(d$cmaq_n_3yroll, alpha=0.9)
plot(jitter(d$sr, factor=2),  d$scr_obs, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Obsvd airscore')
plot(d$cmaq_n_3yroll, d$scr_obs, col=ca, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Obsvd airscore')
plot(d$ubc_mat, d$scr_obs, col=ca,
     xlab='Mean ann temp', ylab='Obsvd airscore')
plot(jitter(d$sr, factor=2),  ci_rng, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Uncertainty')
plot(d$cmaq_n_3yroll, ci_rng, col=ca, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Uncertainty')
plot(d$ubc_mat, ci_rng, col=ca,
     xlab='Mean ann temp', ylab='Uncertainty')
dev.off()

### mean and CI are positively related to richness..... WHY???
`plot_ci` <- function (z, ...) {
  y   <- z[,2]
  lwr <- z[,1]
  upr <- z[,3]
  x   <- 1:NROW(z)
  ylm <- c(min(c(y,lwr),na.rm=T) - 0.4, max(c(y,upr),na.rm=T) + 0.1)
  plot(x, y, pch=16, cex=0.2, ylim=ylm, ylab='Airscore \u00B1 95% CI',
       xaxs='i', yaxs='r', xaxt='n', ...)
  segments(x0=x, x1=x, y0=lwr, y1=upr, lwd=0.1, lend='butt')
}

### plot bootstrapped airscores vs richness, then CMAQ
png('./fig/fig_01_CIs_richness-CMAQ.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(d$sr, d$cmaq_n_3yroll)
plot_ci(ci[o,], xlab='Sites, ordered by richness then CMAQ N dep')
csr <- cumsum(table(d$sr))[1:27]
text(c(0,csr)+diff(c(0,csr,length(d$sr)))*0.5,17,labels=c(4:30,'>30'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()

### plot bootstapped airscores vs observed CMAQ
png('./fig/fig_02_CIs_CMAQ.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(d$cmaq_n_3yroll, d$sr)
plot_ci(ci[o,], xlab='Sites, ordered by CMAQ N dep')
csr <- cumsum(table(round(d$cmaq_n_3yroll[o])))[1:15]
text(c(0,csr)+diff(c(0,csr,length(d$sr)))*0.5, 17, labels=c(1:15,'>15'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()

### zoom-in
png('./fig/fig_01_CIs_zoom_SR11.png',
    wid=4.5, hei=4.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='s')
o   <- order(d$cmaq_n_3yroll)
i   <- which(d$sr==11)
set.seed(121)
i   <- sort(sample(i,size=50))
tci <- ci[o,]
tci <- tci[i,]
`plot_ci` <- function (z, ...) {
  y   <- z[,2]
  lwr <- z[,1]
  upr <- z[,3]
  x   <- 1:NROW(z)
  ylm <- c(min(c(y,lwr),na.rm=T) - 0.4, max(c(y,upr),na.rm=T) + 0.1)
  plot(x, y, pch=16, cex=0.5, ylim=ylm, ylab='Airscore \u00B1 95% CI',
       xaxs='r', yaxs='r', xaxt='n', ...)
  segments(x0=x, x1=x, y0=lwr, y1=upr, lwd=0.5, lend='butt')
}
plot_ci(tci, xlab='Sites, ordered by CMAQ N deposition')
dev.off()

####    END    ####
