######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- species *richness* summary/plotting
#
#    Rob Smith, phytomosaic@gmail.com, 03 Mar 2021
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(ecole)
require(scales)
require(ggplot2)
require(gridExtra)
require(usmap)
require(viridis)

### setup plot labels
nlab <- expression(N~uncertainty~(kg~ha^-1~y^-1))
slab <- expression(S~uncertainty~(kg~ha^-1~y^-1))

### load all bootstrap results --- nitrogen
fnm <- list.files('./res/rich_boot/', pattern='b_[[:digit:]]')
fnm <- fnm[grep('n_', fnm)] # grab nitrogen
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) })
ci_rng_n <- sapply(lst, '[', , 5, drop=T)  # uncertainty as CI95 range
se_n     <- sapply(lst, '[', , 4, drop=T)  # SEM (use only for MQO work)
dimnames(ci_rng_n)[[2]] <- dimnames(se_n)[[2]] <- 2:20
rm(lst, b, fnm)

### load all bootstrap results --- sulfur
fnm <- list.files('./res/rich_boot/', pattern='b_[[:digit:]]')
fnm <- fnm[grep('s_', fnm)] # grab sulfur
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) })
ci_rng_s <- sapply(lst, '[', , 5, drop=T)  # uncertainty as CI95 range
se_s     <- sapply(lst, '[', , 4, drop=T)  # SEM (use only for MQO work)
dimnames(ci_rng_s)[[2]] <- dimnames(se_s)[[2]] <- 2:20
rm(lst, b, fnm)

# ### SEM and range of variability are pretty correlated, as desired
# set_par_mercury(2)
# smoothScatter(se_n, ci_rng_n, colramp = viridis::viridis, col='white',
#               nrpoints = floor(prod(dim(ci_rng_n)) / 100)) # outliers = 1%
# smoothScatter(se_s, ci_rng_s, colramp = viridis::viridis, col='white',
#               nrpoints = floor(prod(dim(ci_rng_s)) / 100)) # outliers = 1%

### load descriptors and prepare
load(file='./res/d_n_scr.rda')   # descriptor matrix
dn <- d
load(file='./res/d_s_scr.rda')   # descriptor matrix
ds <- d
rm(d)
dn$ecoreg1    <- tolower(dn$ecoreg1)
dn$ci_rng_n02 <- ci_rng_n[,'2']
dn$ci_rng_n10 <- ci_rng_n[,'10']
dn$ci_rng_n20 <- ci_rng_n[,'20']
se_n          <- se_n[dn$ecoreg1 != 'water',]
dn            <- dn[dn$ecoreg1 != 'water',]
dn$ecoreg1[dn$ecoreg1 == 'northwestern forested mountains'] <- 'northwest montane forest'
dn$ecoreg1[dn$ecoreg1 == 'southern semi-arid highlands']    <- 'semi-arid highlands'
ds$ecoreg1    <- tolower(ds$ecoreg1)
ds$ci_rng_s02 <- ci_rng_s[,'2']
ds$ci_rng_s10 <- ci_rng_s[,'10']
ds$ci_rng_s20 <- ci_rng_s[,'20']
se_s          <- se_s[ds$ecoreg1 != 'water',]
ds            <- ds[ds$ecoreg1 != 'water',]
ds$ecoreg1[ds$ecoreg1 == 'northwestern forested mountains'] <- 'northwest montane forest'
ds$ecoreg1[ds$ecoreg1 == 'southern semi-arid highlands']    <- 'semi-arid highlands'

### calculate sufficient species
mqo       <- 1.5    # kg ha y
dn$n_suff <- apply(se_n, 1, function(x) min(which(x < mqo), na.rm=T)+1) # min suff
ds$n_suff <- apply(se_s, 1, function(x) min(which(x < mqo), na.rm=T)+1) # min suff
ds$n_suff[is.infinite(ds$n_suff)] <- 20L # some never fell below MQO
dn$sr_margin <- dn$sr - dn$n_suff  # species margin above/below needed meet MQO
ds$sr_margin <- ds$sr - ds$n_suff  # species margin above/below needed meet MQO
sum(se_n < mqo) / prod(dim(se_n))  # 86.6% fall below MQO - nitrogen
sum(se_s < mqo) / prod(dim(se_s))  # 84.4% fall below MQO - sulfur

### MEAN UNCERTAINTY to report
# mean(dn$ci_rng_n02)  # 6.92 kg N ha y (at least generous 2 species)
# mean(dn$ci_rng_n10)  # 3.40 kg N ha y (at moderate 10 species)
mean(dn$ci_rng_n20)    # 2.42 kg N ha y (at most generous 20 species) <-- USE THIS!
# mean(ds$ci_rng_s02)  # 7.01 kg S ha y (at least generous 2 species)
# mean(ds$ci_rng_s10)  # 3.45 kg S ha y (at moderate 10 species)
mean(ds$ci_rng_s20)    # 2.47 kg S ha y (at most generous 20 species) <-- USE THIS!



### --- Fig. 03 --- boxplots of uncertainty vs richness
png('./fig/fig_03_richness_boxplots.png',
    wid=6.5, hei=3.5, uni='in', res=1080, bg='transparent')
set_par_mercury(2, mgp=c(1.3,0.2,0), CEX=0.9)
boxplot(x = as.list(as.data.frame(ci_rng_n)), outcex=0.3, ylim=c(0,22),
        boxwex=0.4, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
        whisklty=1, staplewex=0, cex.axis=0.7,
        xlab='Simulated richness', ylab=nlab)
boxplot(x = as.list(as.data.frame(ci_rng_s)), outcex=0.3, ylim=c(0,22),
        boxwex=0.4, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
        whisklty=1, staplewex=0, cex.axis=0.7,
        xlab='Simulated richness', ylab=slab)
dev.off()


### --- Fig. S1 --- Measurement error: at what richness does SEM cross MQO?
png('./fig/fig_s1_MQO_richness.png',
    wid=6.5, hei=3.5, uni='in', res=1080, bg='transparent')
set_par_mercury(2, mgp=c(1.3,0.2,0))
matplot(t(se_n), type='l', xlab='Plot richness',
        ylab=expression('Standard error ('*kg~N~ha^-1~y^-1*')'),
        lty=1, col='#00000005', ylim=c(0,6))
add_text(0.98, 0.37, 'MQO', pos=2, col=2, cex=0.8)
abline(h=1.5, col=2, lwd=2)
matplot(t(se_s), type='l', xlab='Plot richness',
        ylab=expression('Standard error ('*kg~S~ha^-1~y^-1*')'),
        lty=1, col='#00000005', ylim=c(0,6))
add_text(0.98, 0.37, 'MQO', pos=2, col=2, cex=0.8)
abline(h=1.5, col=2, lwd=2)
dev.off()
rm(mqo, ci_rng_n, ci_rng_s, se_n, se_s) # cleanup

# ### setup for boxplot by ecoregion
# a <- data.frame(aggregate(dn$ci_rng_n20, list(ecoregion=dn$ecoreg1), median),
#                 n = c(table(dn$ecoreg1)))
# (a <- a[rev(order(a$x)),])
# dn$ecoreg1 <- factor(dn$ecoreg1, levels=a$ecoregion)
# grp  <- a$ecoregion
# k    <- length(grp)
# `bxplt` <- function(x, y, do_xaxt=TRUE, CEX=0.7, ...) {
#   plot(as.factor(x$ecoreg1), y, outcex=0.4, # ylim=c(0,500),
#        boxwex=0.3, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
#        whisklty=1, staplewex=0,
#        xlab='Ecoregion', xaxt='n', xaxs='i', pty='s', box.lty=1,
#        mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L', cex=CEX,
#        cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
#   if(isTRUE(do_xaxt)) {
#     incr <- (par('usr')[4] - par('usr')[3]) * 0.04
#     text(1:k+0.1, y=par('usr')[3]-incr, srt=0, adj=1, xpd=T, cex=CEX*0.9,
#          labels=LETTERS[1:k])
#   }
# }
# ### boxplot uncertainties per region
# png('./fig/fig_03_bxplt_unc_richness_by_region.png',
#     wid=6.5, hei=4.75, units='in', bg='transparent', res=1080)
# set_par_mercury(6, mar=c(3.1,3.1,0.5,0.5), oma=c(0,0,0,0))
# # nitrogen
# bxplt(dn, dn$ci_rng_n02, T, ylab='', ylim=c(0,15))
# title(ylab=nlab, line=0.75, cex.lab=0.9)
# add_label('    A    richness = 2', cex=0.8)
# bxplt(dn, dn$ci_rng_n10, T, ylab='', ylim=c(0,15))
# title(ylab=nlab, line=0.75, cex.lab=0.9)
# add_label('    B    richness = 10', cex=0.8)
# bxplt(dn, dn$ci_rng_n20, T, ylab='', ylim=c(0,15))
# title(ylab=nlab, line=0.75, cex.lab=0.9)
# add_label('    C    richness = 20', cex=0.8)
# legend('topright', paste0(LETTERS[1:k], ' = ', grp),
#        col='transparent', border=NA, bty='n', cex=0.46, ncol=2)
# # sulfur
# bxplt(ds, ds$ci_rng_s02, T, ylab='', ylim=c(0,15))
# title(ylab=slab, line=0.75, cex.lab=0.9)
# add_label('    D    richness = 2', cex=0.8)
# bxplt(ds, ds$ci_rng_s10, T, ylab='', ylim=c(0,15))
# title(ylab=slab, line=0.75, cex.lab=0.9)
# add_label('    E    richness = 10', cex=0.8)
# bxplt(ds, ds$ci_rng_s20, T, ylab='', ylim=c(0,15))
# title(ylab=slab, line=0.75, cex.lab=0.9)
# add_label('    F    richness = 20', cex=0.8)
# legend('topright', paste0(LETTERS[1:k], ' = ', grp),
#        col='transparent', border=NA, bty='n', cex=0.46, ncol=2)
# dev.off()
# rm(a, k, grp, bxplt) # cleanup



### --- Fig. xxx --- map of *richness* uncertainties
`plot_map` <- function(dn, zvar='ci_rng_n', legtitle='Title here', tag='',
                       dir=1, brk, transf='identity', ...) {
  plot_usmap(fill='lightgrey') +
    geom_point(data=dn, aes_string(x = 'lon.1', y = 'lat.1', color = zvar), size=0.5) +
    scale_color_viridis(name=legtitle, na.value='transparent',
                        option='D', direction=dir,
                        trans = transf, breaks = brk, labels = brk) +
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
dn <- dn[dn$state != 'Missouri',]         # rm another
dn <- dn[rev(order(dn$sr_margin)),]       # order by species richness margin
dn <- usmap_transform(dn)                 # reproject
ds <- ds[!(ds$lat>49.5 & ds$lon > -122),] # rm one invalid location in Canada
ds <- ds[ds$state != 'Missouri',]         # rm another
ds <- ds[rev(order(ds$sr_margin)),]       # order by species richness margin
ds <- usmap_transform(ds)                 # reproject
# png(file=paste0('./fig/fig_09_map_unc_richness.png'),
#     wid=6.5, hei=5.0, unit='in', bg='transparent', res=1000)
set_par_mercury(1)
grid.arrange(
  # plot_map(dn, 'ci_rng_n02', legtitle=nlab, brk=c(0:5)), #brk=c(0,50,100,200)), #c(5,25,50,100,200,400),
  plot_map(dn, 'ci_rng_n20', legtitle=nlab, brk=c(0:5)), #brk=c(0,50,100,200)), #c(5,25,50,100,200,400),
  plot_map(ds, 'ci_rng_s20', legtitle=slab, brk=c(0:10)),
  ncol=2, widths=c(1,1))
# dev.off()


### --- Fig. 04 --- maps of richness *margin* (surplus/deficit)
`plot_map_diverg` <- function(dn, zvar='sr_margin', legtitle='Richness margin',
                              tag='', dir=1, brk, transf='pseudo_log', ...) {
  plot_usmap(fill='grey90') +
    geom_point(data=dn, aes_string(x='lon.1',y='lat.1',color=zvar), size=0.25) +
    scale_colour_gradient2(name=legtitle, na.value='transparent',
                           trans=transf, breaks=brk, labels=brk, midpoint=0,
                           low='#FF0000', mid='#FFFFFF', high='#0000FF') +
    # # MUTED (poor look):
    # low='#832424', mid='#FFFFFF00', high='#3A3A98'
    guides(size='none',
           color=guide_colourbar(ticks=F, title.position='top',
                                 title.hjust=0.5, barheight=0.2, barwidth=5))+
    theme(legend.title = element_text(size=6),
          legend.text  = element_text(size=4),
          legend.direction='horizontal', legend.position = c(0.53, 0.05),
          plot.background = element_blank(), panel.background = element_blank(),
          legend.background = element_blank(), legend.key  = element_blank(),
          title = element_text(size=14, face='bold'),
          plot.margin = margin(0,0,0,0),
          plot.tag.position = c(0.05, 1)) + labs(tag = tag) +
    theme(text=element_text(family='Routed Gothic', colour='black'))
}
# dn <- dn[!(dn$lat>49.5 & dn$lon > -122),] # rm one invalid location in Canada
# dn <- dn[dn$state != 'Missouri',]        # rm another
dn <- dn[!(dn$sr_margin > 40),]          # rm 12 outliers
# dn <- dn[rev(order(dn$sr_margin)),]      # order by species richness margin
# dn <- usmap_transform(dn)               # reproject
png(file=paste0('./fig/fig_04_map_richnessmargin.png'),
    wid=6.5, hei=3.0, unit='in', bg='transparent', res=1080)
grid.arrange(
  plot_map_diverg(dn, 'sr_margin',
                  legtitle='Richness margin',
                  brk=c(-10,-5,0,5,10,40)),
  plot_map_diverg(ds, 'sr_margin',
                  legtitle='Richness margin',
                  brk=c(-10,-5,0,5,10,40)),
  ncol=2, widths=c(1,1))
dev.off()


####    END    ####
