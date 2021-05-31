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
require(ggplot2)
require(gridExtra)
require(usmap)
require(viridis)

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

### labeling fixes
dn$ecoreg1 <- tolower(dn$ecoreg1)
dn$ecoreg2 <- tolower(dn$ecoreg2)
dn$ecoreg3 <- tolower(dn$ecoreg3)
ds$ecoreg1 <- tolower(ds$ecoreg1)
ds$ecoreg2 <- tolower(ds$ecoreg2)
ds$ecoreg3 <- tolower(ds$ecoreg3)
dn <- dn[dn$ecoreg1 != 'water',]
ds <- ds[ds$ecoreg1 != 'water',]
dn$ecoreg1[dn$ecoreg1 == 'northwestern forested mountains'] <- 'northwest montane forest'
ds$ecoreg1[ds$ecoreg1 == 'northwestern forested mountains'] <- 'northwest montane forest'
dn$ecoreg1[dn$ecoreg1 == 'southern semi-arid highlands'] <- 'semi-arid highlands'
ds$ecoreg1[ds$ecoreg1 == 'southern semi-arid highlands'] <- 'semi-arid highlands'

### *relative* range of variability (percent)
dn$rv <- dn$ci_rng / dn$med * 100
ds$rv <- ds$ci_rng / ds$med * 100

### setup plot labels (absolute and relative uncertainty)
nlaba <- expression(atop(NA, atop(textstyle('N absolute uncertainty'),
                                  textstyle((kg~ha^-1~y^-1)))))
slaba <- expression(atop(NA, atop(textstyle('S absolute uncertainty'),
                                  textstyle((kg~ha^-1~y^-1)))))
nlabr <- expression(atop(NA, atop(textstyle('N relative uncertainty (%)'),
                                  scriptscriptstyle(''))))
slabr <- expression(atop(NA, atop(textstyle('S relative uncertainty (%)'),
                                  scriptscriptstyle(''))))

### calc site *exceedances*: bootstrapped median minus the fixed CL
dn$exc <- dn$med - 1.5   # N airscores CL = 1.5
ds$exc <- ds$med - 2.7   # S airscores CL = 2.7

### split 'eastern temperate forest' ecoregion to north vs south
a <- c('mississippi alluvial and southeast usa coastal plains',
       'southeastern usa plains')
b <- c('southwestern appalachians','blue ridge','ozark highlands','ridge and valley')
dn$ecoreg1[dn$ecoreg2 %in% a | dn$ecoreg3 %in% b] <- 'southern temperate forests'
ds$ecoreg1[ds$ecoreg2 %in% a | ds$ecoreg3 %in% b] <- 'southern temperate forests'
rm(a,b)

### setup boxplot by region
a <- data.frame( # summary table
  aggregate(dn[,c('exc','ci_rng','rv')], list(ecoregion=dn$ecoreg1), median),
  n = aggregate(dn[,'ci_rng'], list(ecoregion=dn$ecoreg1), length)[,2],
  aggregate(ds[,c('exc','ci_rng','rv')], list(ecoregion=ds$ecoreg1), median)[,2:3],
  n = aggregate(ds[,'ci_rng'], list(ecoregion=ds$ecoreg1), length)[,2]
)
(a <- a[rev(order(a$rv)),]) # sort by relative range of variation
dn$ecoreg1 <- factor(dn$ecoreg1, levels=a$ecoregion)
ds$ecoreg1 <- factor(ds$ecoreg1, levels=a$ecoregion)
grp        <- a$ecoregion
k          <- length(grp)
dn$regionlab <- factor(dn$ecoreg1, labels=paste0(LETTERS[1:k],' = ',grp))
ds$regionlab <- factor(ds$ecoreg1, labels=paste0(LETTERS[1:k],' = ',grp))
`bxplt` <- function(x, y, yvar='20', CEX=0.7, ...) {
  plot(as.factor(x$ecoreg1), y[,yvar], outcex=0.4, boxwex=0.3, boxfill='#c1c1c1',
       outpch=16, outcol='#00000010', whisklty=1, staplewex=0,
       xlab='Ecoregion', xaxt='n', xaxs='i', pty='s', box.lty=1,
       mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L', cex=CEX,
       cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
  incr <- (par('usr')[4] - par('usr')[3]) * 0.04
  text(1:k+0.1, y=par('usr')[3]-incr, srt=0, adj=1, xpd=T, cex=CEX*0.9,
       labels=LETTERS[1:k])
}


### export CSVs for kriging in NCLAS tool, for Leah Charash, 26 May 2021
j <- c('lon','lat','ubc_mat','ubc_map','ubc_cmd','elevuse_m','ci_rng')
write.csv(dn[,j], file='./krig/n_for_nclas.csv')
write.csv(dn[,j], file='./krig/s_for_nclas.csv')


# ### --- Fig. xxx --- boxplots of *composition* uncertainties, by region
# png('./fig/fig_08_bxplt_unc_composition_by_region.png',
#     wid=6.5, hei=6.0, units='in', bg='transparent', res=700)
set_par_mercury(4, mar=c(3,4,0.5,0.5), oma=c(0.1,0.1,0,0))
bxplt(dn, dn, 'rv', T, ylab=nlabr, ylim=c(0,300), CEX=0.7)
legend('topleft', paste0(LETTERS[1:k], ' = ', grp),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
bxplt(ds, ds, 'rv', T, ylab=slabr, ylim=c(0,300), CEX=0.7)
bxplt(dn, dn, 'ci_rng', T, ylab=nlaba, ylim=c(0,20), CEX=0.7)
bxplt(ds, ds, 'ci_rng', T, ylab=slaba, ylim=c(0,20), CEX=0.7)
# dev.off()


### --- Fig. xxx --- map of *composition* uncertainties
`plot_map` <- function(d, zvar='rv', legtitle='Title here', tag='',
                       dir=1, brk, transf='identity', ...) {
  plot_usmap(fill='lightgrey') +
    geom_point(data=d, aes_string(x='lon.1', y='lat.1', color=zvar), size=0.5) +
    scale_color_viridis(name=legtitle, na.value='transparent', option='D',
                        direction=dir, trans=transf, breaks=brk, labels=brk) +
    guides(size='none',
           color=guide_colourbar(ticks=F, title.position='top',
                                 title.hjust=0.5, barheight=0.2, barwidth=5))+
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=4),
          legend.direction='horizontal', legend.position = c(0.53, 0.05),
          plot.background = element_blank(), panel.background = element_blank(),
          legend.background = element_blank(), legend.key  = element_blank(),
          title = element_text(size=14),
          plot.margin = margin(0,0,0,0),
          plot.tag.position = c(0.05, 1)) +
    labs(tag = tag) +
    theme(text=element_text(family='Routed Gothic', colour='black'))
}
dn <- dn[!(dn$lat>49.5 & dn$lon > -122),] # rm one invalid location in Canada
ds <- ds[!(ds$lat>49.5 & ds$lon > -122),] # rm one invalid location in Canada
dn <- dn[!(dn$rv > 200),] # rm 12 outliers
ds <- ds[!(ds$rv > 300),] # rm 40 outliers
dn <- usmap_transform(dn)                 # reproject
ds <- usmap_transform(ds)                 # reproject
png(file=paste0('./fig/fig_09_map_unc.png'),
    wid=6.5, hei=5.0, unit='in', bg='transparent', res=1000)
set_par_mercury(1)
grid.arrange(
  plot_map(dn, 'rv', legtitle=nlabr, brk=c(0,50,100,200)), #c(5,25,50,100,200,400),
  plot_map(ds, 'rv', legtitle=slabr,brk=c(0,50,100,200)),
  plot_map(dn, 'ci_rng', legtitle=nlaba, brk=seq(0,16,by=4)),
  plot_map(ds, 'ci_rng', legtitle=slaba, brk=seq(0,16,by=4)),
  ncol=2, widths=c(1,1))
dev.off()


### --- Fig. 1 --- maps of ecoregions
dn$regionlab<-factor(dn$regionlab,levels=levels(dn$regionlab)[k:1]) # no overplot
`plot_map_regions` <- function(d=dn,zvar='regionlab',legtitle='Ecoregions',tag='') {
  plot_usmap(fill='lightgrey') +
    geom_point(data=d, aes_string(x = 'lon.1', y = 'lat.1', color = zvar), size=0.5) +
    scale_color_manual(
      values = c('#de4e4e', '#40ada6','#e87878','#984ea3','#ff7f00',
                 '#ffeb94','#a65628','#377eb8','#479144','#75d7f4'),
      name=legtitle, na.value='transparent') +
    guides(size='none',
           color=guide_legend(title.position='top',title.hjust=0.5,
                              ncol=2, keyheight=0.4)) +
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=4),
          legend.direction='horizontal', legend.position = c(0.43, 0.03),
          plot.background = element_blank(), panel.background = element_blank(),
          legend.background = element_blank(), legend.key  = element_blank(),
          title = element_text(size=14), plot.margin = margin(0,0,0,0),
          plot.tag.position = c(0.05,1)) + labs(tag = tag) +
    theme(text=element_text(family='Routed Gothic', colour='black'))
}
png(file=paste0('./fig/fig_01_map_ecoregions.png'),
    wid=6.5, hei=5.0, unit='in', bg='transparent', res=1000)
plot_map_regions()
dev.off()





### --- NOTRUN --- DIAGNOSTICS ----------------------------------------------

# # correct NA climate values
# dn$ubc_mat[dn$ubc_mat < -99] <- NA
# ds$ubc_mat[ds$ubc_mat < -99] <- NA
#
# ### diagnostics N ........
# # png('./fig/fig_00_diagnostics_n.png',
# #     wid=12, hei=8, uni='in', res=700, bg='transparent')
# set_par_mercury(6)
# u <- colvec(dn$cmaq_n_3yroll, alpha=0.9)
# plot(jitter(dn$sr, factor=2),  dn$scr_obs, col=u, xlim=c(0,40),
#      xlab='Species richness', ylab='Obsvd airscore')
# plot(dn$cmaq_n_3yroll, dn$scr_obs, col=u, xlim=c(0,20),
#      xlab='CMAQ N dep', ylab='Obsvd airscore')
# plot(dn$ubc_mat, dn$scr_obs, col=u,
#      xlab='Mean ann temp', ylab='Obsvd airscore')
# plot(jitter(dn$sr, factor=2),  dn$ci_rng, col=u, xlim=c(0,40),
#      xlab='Species richness', ylab='Uncertainty')
# plot(dn$cmaq_n_3yroll, dn$ci_rng, col=u, xlim=c(0,20),
#      xlab='CMAQ N dep', ylab='Uncertainty')
# plot(dn$ubc_mat, dn$ci_rng, col=u,
#      xlab='Mean ann temp', ylab='Uncertainty')
# # dev.off()
# ### diagnostics S ........
# # png('./fig/fig_00_diagnostics_s.png',
# #     wid=12, hei=8, uni='in', res=700, bg='transparent')
# set_par_mercury(6)
# u <- colvec(ds$cmaq_s_3yroll, alpha=0.9)
# plot(jitter(ds$sr, factor=2),  ds$scr_obs, col=u, xlim=c(0,40),
#      xlab='Species richness', ylab='Obsvd airscore')
# plot(ds$cmaq_s_3yroll, ds$scr_obs, col=u, xlim=c(0,20),
#      xlab='CMAQ S dep', ylab='Obsvd airscore')
# plot(ds$ubc_mat, ds$scr_obs, col=u,
#      xlab='Mean ann temp', ylab='Obsvd airscore')
# plot(jitter(ds$sr, factor=2),  ds$ci_rng, col=u, xlim=c(0,40),
#      xlab='Species richness', ylab='Uncertainty')
# plot(ds$cmaq_s_3yroll, ds$ci_rng, col=u, xlim=c(0,20),
#      xlab='CMAQ S dep', ylab='Uncertainty')
# plot(ds$ubc_mat, ds$ci_rng, col=u,
#      xlab='Mean ann temp', ylab='Uncertainty')
# # dev.off()
#
# ### mean and CI are positively related to richness..... WHY???
# `plot_ci` <- function (z, ...) {
#   y   <- z$med
#   lwr <- z$lwr
#   upr <- z$upr
#   x   <- 1:NROW(z)
#   ylm <- c(min(c(y,lwr),na.rm=T) - 0.4, max(c(y,upr),na.rm=T) + 0.1)
#   plot(x, y, pch=16, cex=0.2, ylim=ylm, ylab='Airscore \u00B1 95% CI',
#        xaxs='i', yaxs='r', xaxt='n', ...)
#   segments(x0=x, x1=x, y0=lwr, y1=upr, lwd=0.1, lend='butt')
# }
#
# ### plot bootstrapped airscores vs richness, then CMAQ
# # png('./fig/fig_01_CIs_richness-CMAQ-n.png',
# #     wid=14.0, hei=5, uni='in', res=700, bg='transparent')
# set_par_mercury(1, pty='m')
# o <- order(dn$sr, dn$cmaq_n_3yroll)
# plot_ci(dn[o,], xlab='Sites, ordered by richness then CMAQ N dep')
# csr <- cumsum(table(dn$sr))[1:27]
# text(c(0,csr)+diff(c(0,csr,length(dn$sr)))*0.5,17,labels=c(4:30,'>30'),
#      cex=0.5)
# abline(v=csr, col=2)
# # dev.off()
# # png('./fig/fig_01_CIs_richness-CMAQ-s.png',
# #     wid=14.0, hei=5, uni='in', res=700, bg='transparent')
# set_par_mercury(1, pty='m')
# o <- order(ds$sr, ds$cmaq_s_3yroll)
# plot_ci(ds[o,], xlab='Sites, ordered by richness then CMAQ S dep')
# csr <- cumsum(table(ds$sr))[1:27]
# text(c(0,csr)+diff(c(0,csr,length(ds$sr)))*0.5,17,labels=c(4:30,'>30'),
#      cex=0.5)
# abline(v=csr, col=2)
# # dev.off()
#
# ### plot bootstapped airscores vs observed CMAQ
# # png('./fig/fig_02_CIs_CMAQ-n.png',
# #     wid=14.0, hei=5, uni='in', res=700, bg='transparent')
# set_par_mercury(1, pty='m')
# o <- order(dn$cmaq_n_3yroll, dn$sr)
# plot_ci(dn[o,], xlab='Sites, ordered by CMAQ N dep')
# csr <- cumsum(table(round(dn$cmaq_n_3yroll[o])))[1:15]
# text(c(0,csr)+diff(c(0,csr,length(dn$sr)))*0.5, 17, labels=c(1:15,'>15'),
#      cex=0.5)
# abline(v=csr, col=2)
# # dev.off()
# # png('./fig/fig_02_CIs_CMAQ-s.png',
# #     wid=14.0, hei=5, uni='in', res=700, bg='transparent')
# set_par_mercury(1, pty='m')
# o <- order(ds$cmaq_s_3yroll, ds$sr)
# plot_ci(ds[o,], xlab='Sites, ordered by CMAQ S dep')
# csr <- cumsum(table(round(ds$cmaq_s_3yroll[o])))[1:15]
# text(c(0,csr)+diff(c(0,csr,length(ds$sr)))*0.5, 17, labels=c(1:15,'>15'),
#      cex=0.5)
# abline(v=csr, col=2)
# # dev.off()
#
# ### zoom-in
# # png('./fig/fig_01_CIs_zoom_SR11.png',
# #     wid=4.5, hei=4.5, uni='in', res=700, bg='transparent')
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
# # dev.off()

####    END    ####
