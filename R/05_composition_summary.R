######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- species *composition* summary/plotting
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

### setup plot labels
nlab <- expression(N~uncertainty~(kg~ha^-1~y^-1))
slab <- expression(S~uncertainty~(kg~ha^-1~y^-1))

### MEAN UNCERTAINTY to report
mean(dn[,c('ci_rng')]) # 3.22 kg N ha y
mean(ds[,c('ci_rng')]) # 3.34 kg S ha y

### --- Fig. 05 --- boxplots of uncertainty vs richness
# png('./fig/fig_05_composition_boxplots.png',
#     wid=6.5, hei=3.5, uni='in', res=1080, bg='transparent')
tiff('./fig/fig_05_composition_boxplots.tif',
    wid=6.5, hei=3.5, uni='in', res=1080, bg='transparent')
dn$spprich_epimac[dn$spprich_epimac < 4] <- NA
set_par_mercury(2, mgp=c(1.3,0.2,0), CEX=0.9)
plot(x=factor(dn$spprich_epimac), y=dn$ci_rng,
     xlim=c(0,39), ylim=c(0,12), outcex=0.3,
     boxwex=0.3, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
     whisklty=1, staplewex=0, cex.axis=0.6,
     xlab='Observed richness', ylab=nlab)
plot(x=factor(ds$spprich_epimac), y=ds$ci_rng,
     xlim=c(0,39), ylim=c(0,12), outcex=0.3,
     boxwex=0.3, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
     whisklty=1, staplewex=0, cex.axis=0.6,
     xlab='Observed richness', ylab=slab)
dev.off()

# ### calc site *exceedances*: bootstrapped median minus the fixed CL
# dn$exc <- dn$med - 1.5   # N airscores CL = 1.5
# ds$exc <- ds$med - 2.7   # S airscores CL = 2.7

# ### split 'eastern temperate forest' ecoregion to north vs south
# a <- c('mississippi alluvial and southeast usa coastal plains',
#        'southeastern usa plains')
# b <- c('southwestern appalachians','blue ridge','ozark highlands','ridge and valley')
# dn$ecoreg1[dn$ecoreg2 %in% a | dn$ecoreg3 %in% b] <- 'southern temperate forests'
# ds$ecoreg1[ds$ecoreg2 %in% a | ds$ecoreg3 %in% b] <- 'southern temperate forests'
# rm(a,b)
#
# ### setup boxplot by region
# a <- data.frame( # summary table
#   aggregate(dn[,c('ci_rng')], list(ecoregion=dn$ecoreg1), median),
#   num_n = aggregate(dn[,'ci_rng'], list(ecoregion=dn$ecoreg1), length)[,2],
#   s_median_ci_rng = aggregate(ds[,c('ci_rng')], list(ecoregion=ds$ecoreg1), median)[,2],
#   num_s = aggregate(ds[,'ci_rng'], list(ecoregion=ds$ecoreg1), length)[,2]
# )
# (a <- a[rev(order(a[,2])),]) # sort by N uncertainty
# dn$ecoreg1 <- factor(dn$ecoreg1, levels=a$ecoregion)
# ds$ecoreg1 <- factor(ds$ecoreg1, levels=a$ecoregion)
# grp        <- a$ecoregion
# k          <- length(grp)
# dn$regionlab <- factor(dn$ecoreg1, labels=paste0(LETTERS[1:k],' = ',grp))
# ds$regionlab <- factor(ds$ecoreg1, labels=paste0(LETTERS[1:k],' = ',grp))
# `bxplt` <- function(x, y, yvar='20', CEX=0.7, ...) {
#   plot(as.factor(x$ecoreg1), y[,yvar], outcex=0.4, boxwex=0.3, boxfill='#c1c1c1',
#        outpch=16, outcol='#00000010', whisklty=1, staplewex=0,
#        xlab='Ecoregion', xaxt='n', xaxs='i', pty='s', box.lty=1,
#        mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L', cex=CEX,
#        cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
#   incr <- (par('usr')[4] - par('usr')[3]) * 0.04
#   text(1:k+0.1, y=par('usr')[3]-incr, srt=0, adj=1, xpd=T, cex=CEX*0.9,
#        labels=LETTERS[1:k])
# }
#
# # ### export CSVs for kriging in NCLAS tool, for Leah Charash, 26 May 2021
# # j <- c('lon','lat','ubc_mat','ubc_map','ubc_cmd','elevuse_m','ci_rng')
# # write.csv(dn[,j], file='./krig/n_for_nclas.csv')
# # write.csv(dn[,j], file='./krig/s_for_nclas.csv')
# # rm(j)
#
#
# ### --- Fig. 05 --- boxplots of *composition* uncertainties, by region
# png('./fig/fig_05_bxplt_unc_composition_by_region.png',
#     wid=6.5, hei=3, units='in', bg='transparent', res=1080)
# set_par_mercury(2, mar=c(3.1,3.1,0.1,0.1), oma=c(0,0,0,0))
# bxplt(dn, dn, 'ci_rng', T, ylab=nlab, ylim=c(0,20), CEX=0.7)
# legend('topleft', paste0(LETTERS[1:k], ' = ', grp),
#        col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
# bxplt(ds, ds, 'ci_rng', T, ylab=slab, ylim=c(0,20), CEX=0.7)
# dev.off()
# rm(a, grp, bxplt) # cleanup


### --- Fig. 06 --- map of *composition* uncertainties
`plot_map` <- function(d, zvar='ci_rng', legtitle='Title here', tag='',
                       dir=1, brk, transf='identity', ...) {
  plot_usmap(fill='lightgrey') +
    geom_point(data=d, aes_string(x='lon.1', y='lat.1', color=zvar), size=0.35) +
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
dn <- dn[dn$state != 'Missouri',]         # rm another
ds <- ds[ds$state != 'Missouri',]         # rm another
dn <- dn[!(dn$ci_rng > 9), ]              # rm 4 outliers, aids visualization
ds <- ds[!(ds$ci_rng > 12),]              # rm 52 outliers, aids visualization
dn <- usmap_transform(dn)                 # reproject
ds <- usmap_transform(ds)                 # reproject
# ### --- Fig. 06 --- map of *composition* uncertainties
# png(file=paste0('./fig/fig_06_map_unc_composition_HORIZONTAL.png'),
#     wid=6.5, hei=2.5, unit='in', bg='transparent', res=1080)
# set_par_mercury(1)
# grid.arrange(
#   plot_map(dn, 'ci_rng', legtitle=nlab, brk=seq(0,9,by=3)),
#   plot_map(ds, 'ci_rng', legtitle=slab, brk=seq(0,12,by=3)),
#   ncol=2, widths=c(1,1))
# dev.off()

### --- Fig. 06 --- map of *composition* uncertainties
png(file=paste0('./fig/fig_06_map_unc_composition.png'),
    wid=5.0, hei=6.0, unit='in', bg='transparent', res=1080)
tiff(file=paste0('./fig/fig_06_map_unc_composition.tif'),
    wid=5.0, hei=6.0, unit='in', bg='transparent', res=1080)
set_par_mercury(1)
grid.arrange(
  plot_map(dn, 'ci_rng', legtitle=nlab, brk=seq(0,9,by=3)),
  plot_map(ds, 'ci_rng', legtitle=slab, brk=seq(0,12,by=3)),
  nrow=2)
dev.off()


# ### --- Fig. 1 --- maps of ecoregions
# # dn$regionlab <- factor(dn$regionlab, levels=levels(dn$regionlab)[k:1])
# `map_regions` <- function(d=dn,zvar='regionlab',legtitle='Ecoregions',tag=''){
#   plot_usmap(fill='lightgrey') +
#     geom_point(data=d, aes_string(x = 'lon.1', y = 'lat.1', color = zvar), size=0.5) +
#     scale_color_manual(
#       values = c('#de4e4e',  # red  = n am deserts
#                  '#984ea3',  # purple = great plains
#                  '#40ada6',  # teal = temperate sierras
#                  '#e87878',  # rosepink = southern temperate forests
#                  '#ff7f00',  # orange = eastern temperate forests
#                  '#a65628',  # brown = semi-arid highlands
#                  '#ffeb94',  # yellow = northern forests
#                  '#377eb8',  # darkblue = marine west coast forest
#                  '#75d7f4',  # paleblue = mediterranean calif
#                  '#479144'), # green = northwest montane forest
#       name=legtitle, na.value='transparent') +
#     guides(size='none',
#            color=guide_legend(title.position='top',title.hjust=0.5,
#                               ncol=2, keyheight=0.4)) +
#     theme(legend.title = element_text(size=6),
#           legend.text = element_text(size=4),
#           legend.direction='horizontal', legend.position = c(0.43, 0.03),
#           plot.background = element_blank(), panel.background = element_blank(),
#           legend.background = element_blank(), legend.key  = element_blank(),
#           title = element_text(size=14), plot.margin = margin(0,0,0,0),
#           plot.tag.position = c(0.05,1)) + labs(tag = tag) +
#     theme(text=element_text(family='Routed Gothic', colour='black'))
# }
# png(file=paste0('./fig/fig_01_map_ecoregions.png'),
#     wid=6.5, hei=5.0, unit='in', bg='transparent', res=1080)
# map_regions()
# dev.off()

####    END    ####
