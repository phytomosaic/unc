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
require(ecole)
require(scales)
require(ggplot2)
require(gridExtra)
require(usmap)
require(viridis)

### load all bootstrap results
fnm <- list.files('./res/rich_boot/', pattern='b_[[:digit:]]')
fnm <- fnm[grep('n', fnm)] # grab nitrogen
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) }) # 30 sec......
### for each list item:
###     - bootstrapped SEM = SD of 999 means
###     - bootstrapped CI = 2.5--97.5 quantiles of 999 means
###     - ci_rng = range of variability
###     - ci_rng / med = *relative* range of variability (percent)

### *relative* range of variability (percent); CI_RNG = 5th col, MED = 2nd col
rv <- sapply(lst, '[', , 5, drop=T) / sapply(lst, '[', , 2, drop=T) * 100
dimnames(rv)[[2]] <- 2:20

### uncertainty (absolute terms); SEM = 4th column
se <- sapply(lst, '[', , 4, drop=T)
dimnames(se)[[2]] <- 2:20
rm(lst, b, fnm)

### setup plot labels (absolute and relative uncertainty)
nlaba <- expression(atop(NA, atop(textstyle('N absolute uncertainty'),
                                  textstyle((kg~ha^-1~y^-1)))))
slaba <- expression(atop(NA, atop(textstyle('S absolute uncertainty'),
                                  textstyle((kg~ha^-1~y^-1)))))
nlabr <- expression(atop(NA, atop(textstyle('N relative uncertainty (%)'),
                                  scriptscriptstyle(''))))
slabr <- expression(atop(NA, atop(textstyle('S relative uncertainty (%)'),
                                  scriptscriptstyle(''))))


# ### SEM and *relative* range of variability are pretty correlated, as desired
# set_par_mercury(1)
# smoothScatter(se, rv, colramp = viridis::viridis, col='white',
#               nrpoints = floor(prod(dim(se)) / 100)) # outliers = 1%

### calculate sufficient species
mqo    <- 1.5    # kg ha y
load(file='./res/d_n_scr.rda', verbose=T)
anyNA(d$ecoreg1) # expect FALSE
d$ecoreg1 <- tolower(d$ecoreg1)
rv     <- rv[d$ecoreg1 != 'water',]
se     <- se[d$ecoreg1 != 'water',]
d      <- d[d$ecoreg1 != 'water',]
sr     <- d$sr
n_suff <- apply(se, 1, function(x) min(which(x < mqo), na.rm=T)+1) # species needed meet MQO
sum(se < mqo) / prod(dim(se)) # 86% fall below MQO
d$sr_margin <- sr - n_suff    # species margin above/below whats needed meet MQO


# ### Measurement error: at what richness does SEM cross MQO?
# png('./fig/fig_10_MQO_richness.png',
#     wid=7.5, hei=4, uni='in', res=700, bg='transparent')
set_par_mercury(2, mgp=c(1.8,0.2,0))
matplot(t(se), type='l', xlab='Plot richness',
        ylab=expression('Standard error ('*kg~ha^-1~y^-1*')'),
        lty=1, col='#00000005')
add_text(0.98, 0.37, 'MQO', pos=2, col=2, cex=0.8)
abline(h=1.5, col=2, lwd=2)
hist(n_suff, breaks=1:19+0.5, xlab='Minimum richness at MQO',
     ylab='N plots', col='grey', main='', xaxs='i', yaxs='i')
box(bty='l')
# dev.off()

### setup for boxplot by ecoregion
a <- data.frame(aggregate(rv[,c('20')], list(ecoregion=d$ecoreg1), median),
                n = c(table(d$ecoreg1)))
(a <- a[rev(order(a$x)),])
d$ecoreg1 <- factor(d$ecoreg1, levels=a$ecoregion)
grp  <- a$ecoregion
k    <- length(grp)
`bxplt` <- function(x, y, yvar='20', do_xaxt=TRUE, CEX=1, ...) {
  plot(as.factor(x$ecoreg1), y[,yvar], outcex=0.4, # ylim=c(0,500),
       boxwex=0.3, boxfill='#c1c1c1', outpch=16, outcol='#00000010',
       whisklty=1, staplewex=0,
       xlab='Ecoregion', xaxt='n', xaxs='i', pty='s', box.lty=1,
       mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L', cex=CEX,
       cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
  if(isTRUE(do_xaxt)) {
    incr <- (par('usr')[4] - par('usr')[3]) * 0.04
    text(1:k+0.1, y=par('usr')[3]-incr, srt=0, adj=1, xpd=T, cex=CEX*0.9,
         labels=LETTERS[1:k])
  }
}

### boxplot RRV uncertainties per region
png('./fig/fig_09_bxplt_unc_richness_by_region.png',
    wid=9.5, hei=6.0, units='in', bg='transparent', res=700)
# relative:
set_par_mercury(6, mar=c(3,4,1.5,0.5), oma=c(5,0,0,0))
bxplt(d, rv, '2', T, ylim=c(0,750), ylab=nlabr)
add_label('A     richness = 2')
bxplt(d, rv, '10', T, ylim=c(0,750), ylab=nlabr)
add_label('B     richness = 10')
bxplt(d, rv, '20', T, ylim=c(0,750), ylab=nlabr)
add_label('C     richness = 20')
legend('topright', paste0(LETTERS[1:k], ' = ', grp),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
# absolute:
bxplt(d, se, '2', T, ylab='', ylim=c(0,5)) ; title(ylab=nlaba, line=0.75)
add_label('D     richness = 2')
bxplt(d, se, '10', T, ylab='', ylim=c(0,5)) ; title(ylab=nlaba, line=0.75)
add_label('E     richness = 10')
bxplt(d, se, '20', T, ylab='', ylim=c(0,5)) ; title(ylab=nlaba, line=0.75)
add_label('F     richness = 20')
legend('topright', paste0(LETTERS[1:k], ' = ', grp),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
dev.off()

# ### NOTRUN: frequency distributions
# u <- viridis::inferno(NCOL(rv), begin=0.1, end=0.85, dir=-1)
# set_par_mercury(1)
# plot(NA, xlim=c(0,300), ylim=c(0,0.02),
#      ylab='Probability density', xlab='Uncertainty (bootstrap CI width)')
# sapply(1:NCOL(rv), function(j) lines(stats::density(rv[,j]), col=u[j]) )
# legend('topright', legend=dimnames(rv)[[2]], fill=u, title='Richness',
#        bty='n', y.intersp=0.8, cex=0.75)




d$se <- se[,dimnames(se)[[2]] == '20'] # hold richness constant at 20 (most generous)
d$rv  <- rv[,dimnames(se)[[2]] == '20'] # hold richness constant at 20 (most generous)
head(d)

### --- Fig. xxx --- map of *richness* uncertainties
`plot_map` <- function(d, zvar='rv', legtitle='Title here', tag='',
                       dir=1, brk, transf='identity', ...) {
  plot_usmap(fill='lightgrey') +
    geom_point(data=d, aes_string(x = 'lon.1', y = 'lat.1', color = zvar), size=0.5) +
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
d <- d[!(d$lat>49.5 & d$lon > -122),] # rm one invalid location in Canada
d <- d[d$state != 'Missouri',]        # rm another
# d <- d[!(d$sr_margin > 40),]          # rm 12 outliers
d <- d[rev(order(d$sr_margin)),]      # order by species richness margin
d <- usmap_transform(d)               # reproject
# png(file=paste0('./fig/fig_09_map_unc_richness.png'),
#     wid=6.5, hei=5.0, unit='in', bg='transparent', res=1000)
set_par_mercury(1)
grid.arrange(
  plot_map(d, 'rv', legtitle=nlabr, brk=c(0,50,100,200)), #c(5,25,50,100,200,400),
  plot_map(d, 'se', legtitle=nlaba, brk=seq(0,16,by=4)),
  ncol=2, widths=c(1,1))
# dev.off()


# # png(file=paste0('./fig/fig_10_map_richnessmargin.png'),
# #     wid=4.5, hei=3.0, unit='in', bg='transparent', res=1000)
# plot_map_diverg(d, 'sr_margin', legtitle=lab, brk=c(-10,-5,0,5,10,40))
# # dev.off()













### maps of richness *margin* (surplus/deficit)
`plot_map_diverg` <- function(d, zvar='sr_margin', legtitle='Richness margin',
                              tag='', dir=1, brk, transf='pseudo_log', ...) {
  plot_usmap(fill='grey90') +
    geom_point(data=d, aes_string(x='lon.1',y='lat.1',color=zvar), size=0.25) +
    scale_colour_gradient2(name=legtitle, na.value='transparent',
                           trans = transf, breaks = brk, labels = brk,
                           # low = '#FF0000',
                           # mid = '#FFFFFF00',
                           # high = '#0000FF',
                           low =  '#832424',
                           mid = '#FFFFFF00',
                           high = '#3A3A98',
                           # high='#472D7BFF',
                           # mid='#FFFFFF00',
                           # low='#AADC32FF',
                           midpoint = 0) +
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
d <- d[!(d$lat>49.5 & d$lon > -122),] # rm one invalid location in Canada
d <- d[d$state != 'Missouri',]        # rm another
d <- d[!(d$sr_margin > 40),]          # rm 12 outliers
d <- d[rev(order(d$sr_margin)),]      # order by species richness margin
d <- usmap_transform(d)               # reproject
lab <- 'Richness margin'
# png(file=paste0('./fig/fig_10_map_richnessmargin.png'),
#     wid=4.5, hei=3.0, unit='in', bg='transparent', res=1000)
plot_map_diverg(d, 'sr_margin', legtitle=lab, brk=c(-10,-5,0,5,10,40))
# dev.off()






####    END    ####
