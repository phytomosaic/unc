######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- species *richness* summary/plotting
#
#    Rob Smith, phytomosaic@gmail.com, 03 Mar 2021
#
##      GNU General Public License, Version 3.0    ###################


rm(list=ls())
require(ecole)
require(scales)


### load all bootstrap results
fnm <- list.files('./res/rich_boot/', pattern='b_[[:digit:]]')
fnm <- fnm[grep('n', fnm)] # grab nitrogen
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) }) # 30 sec......
### for each list item
###     - bootstrapped SEM = SD of 999 means
###     - bootstrapped CI = 2.5--97.5 quantiles of 999 means

# pick the 'sem' column only
b <- sapply(lst, '[', , 4, drop=TRUE) # SEM = 4th column, CI_RNG = 5th column
dimnames(b)[[2]] <- 2:20
rm(fnm, lst)

### calculate sufficient species
mqo    <- 1.5 # kg ha y
load(file='./res/d_n_scr.rda', verbose=T)
d$ecoreg1 <- tolower(d$ecoreg1)
b      <- b[d$ecoreg1 != 'water',]
d      <- d[d$ecoreg1 != 'water',]
sr     <- d$sr
n_suff <- apply(b, 1, function(x) min(which(x < mqo), na.rm=T)+1) # species needed meet MQO
sum(b < mqo) / prod(dim(b)) # 86% fall below MQO
sr_margin <- sr - n_suff   # species margin above/below whats needed meet MQO

# ### SEM across varying sample sizes (and MQO)
# png('./fig/fig_00_SEM_curves.png',
#     wid=7.5, hei=4, uni='in', res=700, bg='transparent')
set_par_mercury(2, mgp=c(1.8,0.2,0))
matplot(t(b), type='l', xlab='Plot richness',
        ylab=expression('Standard error ('*kg~ha^-1~y^-1*')'),
        lty=1, col='#00000005')
add_text(0.98, 0.37, 'MQO', pos=2, col=2, cex=0.8)
abline(h=1.5, col=2, lwd=2)
hist(n_suff, breaks=1:20, xlab='Minimum richness at MQO',
     ylab='N plots', col='grey', main='', xaxs='i', yaxs='i')
box(bty='l')
# dev.off()
# rm(b)

### summary table by ecoregion
a <- data.frame(
  aggregate(b[,c('20')], list(ecoregion=d$ecoreg1), median),
  n = c(table(d$ecoreg1))
)
(a <- a[rev(order(a$x)),])
d$ecoreg1 <- factor(d$ecoreg1, levels=a$ecoregion)


### boxplot uncertainties (regional)
grp  <- a$ecoregion
k    <- length(grp)
nlab <- expression(Uncertainty~(kg~N~ha^-1~y^-1))
# png('./fig/fig_08_bxplt_unc_richness_by_region.png',
#     wid=6.5, hei=3.0, units='in', bg='transparent', res=700)
`bxplt` <- function(x=d, y=b, yvar='20', do_xaxt=TRUE, ...) {
  CEX <- 1
  b[,'20']
  plot(as.factor(x$ecoreg1), y[,yvar], outcex=0.4, # ylim=c(0,500),
       boxwex=0.3, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
       whisklty=1, staplewex=0,
       xlab='Ecoregion', xaxt='n', xaxs='i', pty='s', box.lty=1,
       mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L', cex=CEX,
       cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
  if(isTRUE(do_xaxt)) {
    text(1:k+0.1, y=par('usr')[3]-0.2, srt=0, adj=1, xpd=T, cex=CEX*0.9,
         labels=LETTERS[1:k])
  }
}
set_par_mercury(3, mar=c(1,4,0.5,0.5), oma=c(5,0,0,0))
bxplt(d, b, '2', T, ylab=nlab, ylim=c(0,5))
add_label('A')
bxplt(d, b, '10', T, ylab=nlab, ylim=c(0,5))
add_label('B')
bxplt(d, b, '20', T, ylab=nlab, ylim=c(0,5))
add_label('C')
legend('topright', paste0(LETTERS[1:k], ' = ', grp),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
# dev.off()

# ### frequency distributions
# u <- viridis::inferno(NCOL(b), begin=0.1, end=0.85, dir=-1)
# set_par_mercury(1)
# plot(NA, xlim=c(0,4), ylim=c(0,2.7),
#      ylab='Probability density', xlab='Uncertainty (bootstrap CI width)')
# sapply(1:NCOL(b), function(j) {
#   d <- stats::density(b[,j])
#   # d$y <- d$y / (max(d$y))   # standardized 0-1
#   lines(d, col=u[j])
# })
# legend('topright', legend=dimnames(b)[[2]], fill=u, title='Richness',
#        bty='n', y.intersp=0.8, cex=0.75)
# abline(v=mqo, lty=2)


# # #####################################################################
# png('./fig/fig_09_map_sufficiency.png',
#     wid=6.5, hei=3.5, uni='in', res=700, bg='transparent')
`h` <- function(x,  ...) {
  hist(x, main='', xaxs='i', yaxs='i', ...) ; box(bty='l')
}
sr_margin[sr_margin > 34] <- NA # hack to soften extreme high values
brk <- seq(min(sr_margin, na.rm=T), max(sr_margin, na.rm=T), by=2)
mid <- which(brk==0)-1
n   <- length(brk)
col <- c(colorRampPalette(c('#f00000aa', '#ffffffaa'), alpha=T)(mid),
         colorRampPalette(c('#ffffffaa', '#00008aff'), alpha=T)(n-mid))
`co` <- function (x, col, brk, ...) {col[as.numeric(cut(x, brk, incl=T))]}
### MAPS
set_par_mercury(1, pty='m', mar=c(1.1,1.1,0,0))
plot(d$lon,d$lat,type='n',
     # asp=1.8,
     ylim=c(25,51.5),xlim=c(-125,-50),ylab='',xlab='')
maps::map('usa', add=T)
points(d$lon, d$lat, pch=16, cex=0.4, col=co(sr_margin, col=col, brk=brk))
maps::map('state', add=T, lwd=0.5, col='#00000050')
### INSET layout
ppar <- par('plt')
xmin <- ppar[2] - (ppar[2] - ppar[1]) * 0.32
xmax <- ppar[2] - (ppar[2] - ppar[1]) * 0.02
ymin <- ppar[3] + (ppar[4] - ppar[3]) * 0.05
ymax <- ppar[3] + (ppar[4] - ppar[3]) * 0.5
op <- par(fig=c(xmin,xmax,ymin,ymax), mar=c(1.6,0,0,0),
          oma=c(0,0,0,0), mgp=c(0.5,-0.2,0), tcl=-0.05, new=T)
h(sr_margin, col = col, breaks = brk, yaxt='n',
  xlab='Richness margin', ylab='', cex.lab=0.7, cex.axis=0.5)
abline(v=0, col=2, lwd=2, lty=2)
add_text(-0.04,0.93,'\u2190 Deficit', pos=4, cex=0.7)
add_text(0.90,0.93,'Surplus \u2192', pos=2, cex=0.7)
par(op)
# dev.off()


# ### create geographic regions -- DONT BOTHER, theyre pretty similar
# set.seed(97)
# k   <- 14
# grp <- stats::kmeans(d[,c('lon','lat')], k)$cluster


####    END    ####
