######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- species *composition* summary/plotting
#
#    Rob Smith, phytomosaic@gmail.com, 03 Mar 2021
#
##      GNU General Public License, Version 3.0    ###################

rm(list=ls())
require(ecole)
require(scales)

### load bootstrapped 95% CI for site airscores (each row = one site)
load(file='./res/b_n_scr.rda')   # bootstrap CIs
bn  <- b                         # bootstrap CIs for N
load(file='./res/b_s_scr.rda')   # bootstrap CIs
bs  <- b                         # bootstrap CIs for S
load(file='./res/d_n_scr.rda')   # descriptor matrix
dn <- cbind(d,bn)
load(file='./res/d_s_scr.rda')   # descriptor matrix
ds <- cbind(d,bs)
rm(b,d,bs,bn)
dn$ecoreg1 <- tolower(dn$ecoreg1)
ds$ecoreg1 <- tolower(ds$ecoreg1)
dn <- dn[dn$ecoreg1 != 'water',]
ds <- ds[ds$ecoreg1 != 'water',]

### calc site *exceedances*: bootstrapped median minus the fixed CL
dn$exc <- dn$med - 1.5   # N airscores CL = 1.5
ds$exc <- ds$med - 2.7   # S airscores CL = 2.7

# ### frequency distributions of uncertainties and exceedances (nationwide)
# set_par_mercury(4)
# hist(dn$exc, breaks=seq(-5, max(ds$exc)+0.5,by=0.5)) # N fairly tight
# hist(ds$exc, breaks=seq(-5, max(ds$exc)+0.5,by=0.5)) # S ranges higher
# hist(dn$ci_rng, breaks=seq(-1, max(ds$ci_rng)+0.05,by=0.5)) # N fairly tight
# hist(ds$ci_rng, breaks=seq(-1, max(ds$ci_rng)+0.05,by=0.5)) # S ranges higher

### summary table by ecoregion
a <- data.frame(
  aggregate(dn[,c('exc','ci_rng')], list(ecoregion=dn$ecoreg1), median),
  n = aggregate(dn[,'ci_rng'], list(ecoregion=dn$ecoreg1), length)[,2],
  aggregate(ds[,c('exc','ci_rng')], list(ecoregion=ds$ecoreg1), median)[,2:3],
  n = aggregate(ds[,'ci_rng'], list(ecoregion=ds$ecoreg1), length)[,2]
)
(a <- a[rev(order(a$ci_rng)),])
dn$ecoreg1 <- factor(dn$ecoreg1, levels=a$ecoregion)
ds$ecoreg1 <- factor(ds$ecoreg1, levels=a$ecoregion)

### boxplot uncertainties (regional)
grp <- a$ecoregion
k   <- length(grp)
u   <- viridis::inferno(k, begin=0.1, end=0.85, dir=-1)
`bxplt` <- function(x=dn, y=dn, yvar='exc', do_xaxt=TRUE, ...) {
  CEX <- 1
  plot(as.factor(x$ecoreg1), y[,yvar], outcex=0.4, # ylim=c(0,500),
       boxwex=0.3, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
       whisklty=1, staplewex=0,
       xlab='Ecoregion', xaxt='n', xaxs='i', pty='s', box.lty=1,
       mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L', cex=CEX,
       cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
  if(isTRUE(do_xaxt)) {
    text(1:k+0.1, y=par('usr')[3]-0.7, srt=0, adj=1, xpd=T, cex=CEX*0.9,
         labels=LETTERS[1:k])
  }
}
# png('./fig/fig_08_bxplt_unc_composition_by_region.png',
#     wid=6.5, hei=3.0, units='in', bg='transparent', res=700)
set_par_mercury(2, mar=c(3,4,0.5,0.5), oma=c(0.1,0.1,0,0))
nlab <- expression(Uncertainty~(kg~N~ha^-1~y^-1))
slab <- expression(Uncertainty~(kg~S~ha^-1~y^-1))
bxplt(dn, dn, 'ci_rng', T, ylab=nlab, ylim=c(0,16))
legend('topright', paste0(LETTERS[1:k], ' = ', grp),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
bxplt(ds, ds, 'ci_rng', T, ylab=slab, ylim=c(0,16))
# dev.off()







### plot maps of exceedances for presentation
u <- c('#5E4FA2', '#4F61AA','#4173B3', '#3386BC','#4198B6', '#51ABAE',
       '#62BEA6', '#77C8A4','#8ED1A4', '#A4DAA4','#B8E2A1', '#CBEA9D',
       '#DEF199', '#EAF69F','#F2FAAC', '#FAFDB8', '#FEFAB6', '#FEF0A5',
       '#FEE695', '#FDD985', '#FDC978', '#FDB96A','#FCA75E', '#F99254',
       '#F67D4A', '#F26943','#E85A47', '#DE4B4B', '#D33C4E', '#C1284A',
       '#AF1446', '#9E0142', rep('transparent',2))
u

`map_it` <- function(x, whichcol = 'exc', pal=u, ...) {
  plot(x$lon, x$lat, type='n', asp=1.6, ylab='', xlab='', ...)
  maps::map('usa', fill=T, col='white', add=T)
  points(x$lon, x$lat, pch=16, cex=0.5, col=colvec(x[,whichcol], pal=u))
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
           xaxt='n', xlab='Longitude\u00B0', yaxt='n',
           ylab='Longitude\u00B0', main='',
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
u <- colvec(dn$cmaq_n_3yroll, alpha=0.9)
plot(jitter(dn$sr, factor=2),  dn$scr_obs, col=u, xlim=c(0,40),
     xlab='Species richness', ylab='Obsvd airscore')
plot(dn$cmaq_n_3yroll, dn$scr_obs, col=u, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Obsvd airscore')
plot(dn$ubc_mat, dn$scr_obs, col=u,
     xlab='Mean ann temp', ylab='Obsvd airscore')
plot(jitter(dn$sr, factor=2),  dn$ci_rng, col=u, xlim=c(0,40),
     xlab='Species richness', ylab='Uncertainty')
plot(dn$cmaq_n_3yroll, dn$ci_rng, col=u, xlim=c(0,20),
     xlab='CMAQ N dep', ylab='Uncertainty')
plot(dn$ubc_mat, dn$ci_rng, col=u,
     xlab='Mean ann temp', ylab='Uncertainty')
dev.off()
### diagnostics S ........
png('./fig/fig_00_diagnostics_s.png',
    wid=12, hei=8, uni='in', res=700, bg='transparent')
set_par_mercury(6)
u <- colvec(ds$cmaq_s_3yroll, alpha=0.9)
plot(jitter(ds$sr, factor=2),  ds$scr_obs, col=u, xlim=c(0,40),
     xlab='Species richness', ylab='Obsvd airscore')
plot(ds$cmaq_s_3yroll, ds$scr_obs, col=u, xlim=c(0,20),
     xlab='CMAQ S dep', ylab='Obsvd airscore')
plot(ds$ubc_mat, ds$scr_obs, col=u,
     xlab='Mean ann temp', ylab='Obsvd airscore')
plot(jitter(ds$sr, factor=2),  ds$ci_rng, col=u, xlim=c(0,40),
     xlab='Species richness', ylab='Uncertainty')
plot(ds$cmaq_s_3yroll, ds$ci_rng, col=u, xlim=c(0,20),
     xlab='CMAQ S dep', ylab='Uncertainty')
plot(ds$ubc_mat, ds$ci_rng, col=u,
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
