######################################################################
#
#  CLAD WG-2 -- estimating uncertainty
#
#    Rob Smith, phytomosaic@gmail.com, 12 June 2020
#
##      GNU General Public License, Version 3.0    ###################
#
#
# Goal: estimate uncertainty in exceedance values.
#
#    EXCEEDANCE = OBS - CL
#
#    OBS = lichen airscore (or CMAQ)
#
#    CL = from Diversity paper nonlinear equation
#
#
# Approach:
#
#    bootstrap species that enter the lichen airscore
#
#    sample from the SE of the nonlinear equation
#
# ####################################################################

require(ecole)
require(labdsv)
require(scales)

### load and clean LICHEN SPECIES occurrence data
rm(list=ls())
s  <- read.csv('./data_raw/megadb_lich.csv', stringsAsFactors=F)
names(s) <- clean_text(names(s), lower=T)
s  <- s[s$growthform == 'Epiphytic macrolichen',] # keep only macros
s  <- s[,c('megadbid','sci_22chklst','fia_abun',
           'n_depmaxfreq','s_depmaxfreq')]
s$sci_22chklst <- clean_text(s$sci_22chklst, TRUE)

### get species ratings ('peak detection frequency' for N or S)
wa <- s[!duplicated(s$sci_22chklst),
        c('sci_22chklst','n_depmaxfreq','s_depmaxfreq')]
wa <- wa[order(wa$sci_22chklst),]

### reshape wide for abundance matrix
spe <- labdsv::matrify(s[,c('megadbid','sci_22chklst','fia_abun')])
identical(wa$sci_22chklst, dimnames(spe)[[2]]) # expect TRUE
rm(s)

### convert spe to trait (ratings to be sampled)
x <- spe
x[x > 0 ] <- 1
x <- sweep(x, 2, wa$n_depmaxfreq, '*')    # nitrogen ratings
# x <- sweep(x, 2, wa$s_depmaxfreq, '*')   # sulfur ratings
x[is.na(x)] <- 0L # zero out any NA (species that lacked ratings)
j   <- colSums(x)==0 # 287 spp lacking valid spp ratings
spe <- spe[,!j]
x   <-   x[,!j]
sr  <- apply(x, 1, function(x) sum(x>0))
i   <- sr < 4 # 1706 species-poor sites (< 4 rated species)
spe <- spe[!i,]
x   <-   x[!i,]
rm(i,j,wa)
x   <- as.matrix(x)

### some summary statistics
scr_obs <- apply(x, 1, function(i) mean(i[i>0])) # obs airscore
sr <- apply(x, 1, function(x) sum(x>0)) # obs richness of RATED spp

### some descriptive data (CMAQ and original airscores)
d <- read.csv('./data_raw/MegaDbPLOT_2020.05.09.csv',
              stringsAsFactors=F)
names(d)  <- clean_text(names(d), lower=T)
i         <- match(dimnames(x)[[1]], d$megadbid)
naq       <- d$cmaq_n_3yroll[i]
saq       <- d$cmaq_s_3yroll[i]
nairscore <- d$n_airscore[i]
sairscore <- d$s_airscore[i]
mat       <- d$ubc_mat[i]
lon       <- d$longusedd[i]
lat       <- d$latusedd[i]
lon[lat > 49.01] <- NA
lat[lat > 49.01] <- NA
rm(d)

### CHOICE -- dont plot cases that lack CMAQ info?...
isna      <- is.na(naq)
naq       <- naq[!isna]
saq       <- saq[!isna]
nairscore <- nairscore[!isna]
sairscore <- sairscore[!isna]
scr_obs   <- scr_obs[!isna]
sr        <- sr[!isna]
lat       <- lat[!isna]
lon       <- lon[!isna]
mat       <- mat[!isna]

# ### get climate adjusted airscores
# load('./data/d.rda')
# nadj      <- d$nadj[match(dimnames(x)[[1]], d$megadbid)]
# nadj      <- nadj[!isna]
# rm(d)

# ######################################################################
# ######################################################################
# ######################################################################
# ### for CLAD WG-2: estimate bootstrap 95% CI of the estimated airscore
# ###   sample positive values (presences) from each row of x,
# ###      with probability equal to species abundances in spe;
# ###   airscore = mean of these randomly sampled trait values
# ###      (inherently community-weighted by probabilities)
# `bootscr` <- function(x, i, B=999) {
#     x    <- x[i, ]      # subset row
#     pr   <- spe[i,] / sum(spe[i,]) # probabilities per row
#     # B    <- 999         # number of bootstrap draws per row
#     n    <- sum(pr > 0) # sample size = observed richness
#     res  <- replicate(B, mean(sample(x,size=n,replace=T,prob=pr)))
#     res # the bootstrapped airscores
# }
# ### ! ! ! TIMEWARN ! ! ! 10 min for 8927 plots (346 spp)
# cat(paste0('start time: ', time_start <- Sys.time()), '\n')
# b_scr <- sapply(1:NROW(spe), FUN=bootscr, x=x, simplify='array')
# dimnames(b_scr)[[2]] <- dimnames(spe)[[1]]
# cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
# save(b_scr, file='./data/b_scr.rda')
# ######################################################################
# ######################################################################
# ######################################################################

### load bootstraps and calc bootstrapped 95% CI for site scores
load(file='./data/b_scr.rda', verbose=T)
ci <- t(apply(b_scr, 2, function(x) quantile(x, c(0.025,0.50,0.975))))
ci <- ci[!isna,]
ci_rng <- ci[,3] - ci[,1] # calc breadth of CIs
exc <- ci[,2] - 3.5 # calc species richness N exceedance (CL = 3.5)
exc <- exc - 2 # downweight fuzz...

### plot maps of exceedances for presentation
ca <- c('#5E4FA2', '#4F61AA','#4173B3', '#3386BC','#4198B6', '#51ABAE',
        '#62BEA6', '#77C8A4','#8ED1A4', '#A4DAA4','#B8E2A1', '#CBEA9D',
        '#DEF199', '#EAF69F','#F2FAAC', '#FAFDB8', '#FEFAB6', '#FEF0A5',
        '#FEE695', '#FDD985', '#FDC978', '#FDB96A','#FCA75E', '#F99254',
        '#F67D4A', '#F26943','#E85A47', '#DE4B4B', '#D33C4E', '#C1284A',
        '#AF1446', '#9E0142', rep('transparent',2))
`map_it` <- function(x, pal=ca, ...) {
  plot(lon, lat, type='n', asp=1.6, ylim=c(30,49),
       ylab='', xlab='', ...)
  maps::map('usa', fill=T, col='white', add=T)
  points(lon, lat, pch=16, cex=0.5, col=colvec(x, pal=ca))
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
png('C:/Users/Rob/Desktop/fig_01_map_exceedances.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(x=exc, main='Airscore N exceedance')
dev.off()
### exceedance uncertainty
png('C:/Users/Rob/Desktop/fig_01_map_exceedance_unc.png',
    wid=5, hei=3.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
map_it(x=ci_rng, main='CI breadth for N exceedance')
dev.off()


### DIAGNOSTICS

### diagnostics........
png('C:/Users/Rob/Desktop/fig_00_diagnostics.png',
    wid=12, hei=8, uni='in', res=700, bg='transparent')
set_par_mercury(6)
ca <- colvec(naq, alpha=0.9)
plot(jitter(sr, factor=2),  scr_obs, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Obsvd airscore')
plot(naq, scr_obs, col=ca, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Obsvd airscore')
plot(mat, scr_obs, col=ca,
     xlab='Mean ann temp', ylab='Obsvd airscore')
plot(jitter(sr, factor=2),  ci_rng, col=ca, xlim=c(0,40),
     xlab='Species richness', ylab='Uncertainty')
plot(naq, ci_rng, col=ca, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Uncertainty')
plot(mat, ci_rng, col=ca,
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
png('C:/Users/Rob/Desktop/fig_01_CIs_richness-CMAQ.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(sr, naq)
plot_ci(ci[o,], xlab='Sites, ordered by richness then CMAQ N dep')
csr <- cumsum(table(sr))[1:27]
text(c(0,csr)+diff(c(0,csr,length(sr)))*0.5,17,labels=c(4:30,'>30'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()

### plot bootstapped airscores vs observed CMAQ
png('C:/Users/Rob/Desktop/fig_02_CIs_CMAQ.png',
    wid=14.0, hei=5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='m')
o <- order(naq, sr)
plot_ci(ci[o,], xlab='Sites, ordered by CMAQ N dep')
csr <- cumsum(table(round(naq[o])))[1:15]
text(c(0,csr)+diff(c(0,csr,length(sr)))*0.5, 17, labels=c(1:15,'>15'),
     cex=0.5)
abline(v=csr, col=2)
dev.off()

### zoom-in
png('C:/Users/Rob/Desktop/fig_01_CIs_zoom_SR11.png',
    wid=4.5, hei=4.5, uni='in', res=700, bg='transparent')
set_par_mercury(1, pty='s')
o   <- order(naq)
i   <- which(sr==11)
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

ls()

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ### for sampling sufficiency: estimate SEM across varying sample sizes
# ### sample positive values (presences) from each row of x,
# ###   with probability equal to species abundances in spe
# ###   airscore = mean of these randomly sampled trait values
# ###   (inherently community-weighted by probabilities)
# ###
# ### bootstrap SEM per site
# `bootsem` <- function(x, i) {
#     x    <- x[i, ] # subset row
#     pr   <- spe[i,] / sum(spe[i,]) # probabilities per row
#     B    <- 999    # number of bootstrap draws per n
#     n    <- 2:30   # increasingly large sample sizes
#     sem  <- function(x) sd(x) # / sqrt(length(x)-1)
#     bsem <- function(ni) replicate(B,sem(sample(x,size=ni,replace=T,prob=pr)))
#     res  <- sapply(n, bsem)
#     dimnames(res)[[2]] <- paste0(n)
#     res
# }
`bootscr` <- function(x, spe, i, B=999, n) {
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
# for (N in 2:20){
#     b <- sapply(1:NROW(spe),
#                 FUN = bootscr,
#                 x   = x,
#                 spe = spe,
#                 n   = N,
#                 simplify='array')
#     dimnames(b)[[3]] <- dimnames(spe)[[1]]
#     save(b, file=paste0('./data/b_',N,'.rda'))
#     rm(b)
# }
# cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
999 * 8927 * 19 # size of FULL array, if using all plots...

### load all bootstrap results
fnm <- list.files('./data/', pattern='b_[0123456789]')
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./data/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) }) # 30 sec......
### calc for each list item
###     - bootstrapped SEM = SD of 999 means
###     - bootstrapped CI = 2.5--97.5 quantiles of 999 means
`f` <- function(b) {
  t(apply(b, 2, function(x) {
    c(sem = sd(x), ci_rng = unname(diff(quantile(x,c(0.025,0.975)))))
  }))
}
b <- lapply(lst, f)
rm(lst)
b <- do.call(cbind, b)
b_sem <- t(b[,seq(1,NCOL(b),2)])
# b_ci  <- t(b[,seq(2,NCOL(b),2)])

### calculate sufficient species
mqo    <- 1.5 # kg ha y
sr     <- apply(spe, 1, function(x) sum(x>0))
n_suff <- apply(b_sem, 2, function(x) min(which(x < mqo), na.rm=T)+1)
sr_deficit <- sr - n_suff # species needed meet MQO

### SEM across varying sample sizes (and MQO)
png('C:/Users/Rob/Desktop/fig_00_SEM_curves.png',
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
png('C:/Users/Rob/Desktop/fig_01_min_sufficient.png',
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
# plot(lon, lat, col=co(sr_deficit), asp=1.6,
# main='', xlab='Richness sufficiency', ylab='')
`map_it` <- function(x, ...) {
  plot(lon, lat, type='n', asp=1.6, ylim=c(30,49),
       ylab='', xlab='Richness sufficiency', ...)
  maps::map('usa', add=T)
  points(lon, lat, pch=16, cex=0.5, col=co(sr_deficit))
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
png('C:/Users/Rob/Desktop/fig_01_map_sufficiency.png',
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
plot(lon, lat, type='n', asp=1.6, ylim=c(30,49), xlim=c(-125,-50),
     ylab='', xlab='')
maps::map('usa', add=T)
points(lon, lat, pch=16, cex=0.5, col=co(sr_deficit))
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
png('C:/Users/Rob/Desktop/fig_01_map_sufficiency_rainbow.png',
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
plot(lon, lat, type='n', asp=1.6, ylim=c(30,49), xlim=c(-125,-50),
     ylab='', xlab='')
maps::map('usa', add=T)
points(lon, lat, pch=16, cex=0.5, col=co(sr_deficit))
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
