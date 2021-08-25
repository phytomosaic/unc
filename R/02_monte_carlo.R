######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- Monte Carlo parameter resampling
#
#    Rob Smith, phytomosaic@gmail.com, 16 Jun 2021
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(ecole)        # for plotting and name handling
require(quantreg)     # for 90th percentile regressions
load('./data/d.rda')  # site data

# ### NB: PRN originally used N from calibration eqn? here we use CMAQ instead
# hist(d$cmaq_n_3yroll, breaks=seq(0,30,by=0.5), col='#00000050')
# hist(d$n_lich_kghay, breaks=seq(0,30,by=0.5), add=T, col='#00000050')
# ### PRN used thresholds for each model: N=12 and S=20 (kg ha y)

### rename columns
onm <- c('spprich_total','spprich_oligo','spprich_s_sens','abun_cyano',
         'abun_forage', 'cmaq_n_3yroll', 'cmaq_s_3yroll')
nnm <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano',
         'abun_forage', 'N', 'S')
names(d)[match(onm,names(d))] <- nnm
nnm <- c(nnm,'ubc_td','lat','lon') # keep if has lat/lon and climate info

### fill NA to 0 in diversity measures
j <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano','abun_forage')
d[,j] <- sapply(d[,j], function(x) { x[is.na(x)] <- 0 ; return(x) })

# ### pick by PRN's megadbid
# u <- paste0('https://raw.githubusercontent.com/phytomosaic/Lichen_CL_quantile',
#             '_regression/master/output/LichDb_sFncGrpSens_All_Abun.csv')
# fnm <- tempfile()
# download.file(u, fnm, method = 'curl')
# f <- read.csv(fnm)
# d <- d[d$megadbid %in% f$megadbid, nnm] # keep only megadb's from PRN

### remove nonstandard and species-poor plots
d <- d[d$plottype == 'Standard',]
d <- d[d$fia_prot==1,]
d <- d[,nnm]              # keep only used columns
d <- na.omit(d)           # keep only complete rows (omits w/o climate)
d <- d[d$spp_rich > 4, ]  # drop species-poor plots
range(d$lat)              # AK and CONUS
dim(d)                    # 6024 sites included

### Monte Carlo resampling of parametric coefficients
`monte_carlo` <- function(fmla,n=99999,do_plot=T,...) {
   cat('now doing', format(fmla), '...\n')
   xvar   <- substr(gsub('.*poly\\(', '', format(fmla)), 1,1)
   sek    <- seq(0.1,  max(d[,xvar], na.rm=T), by=0.01)
   mod    <- quantreg::rq(fmla,tau=0.90,data=d) # focal model
   crit   <- 0.20                               # critical decline = 20%
   newdat <- data.frame(sek)                    # sequence
   colnames(newdat) <- if(grepl('poly\\(S', paste0(fmla)[3])) 'S' else 'N'
   pr  <- predict(mod,newdat,type='none',interval='conf') # *predicted* vals
   `nadir` <- function(x) { # trim to parabola nadir (forbid increasing curve)
      is_decr <- sapply(2:length(x), function(i) (x[i] < x[i-1]) * 1)
      if(all(is_decr == 1)) length(x) else which.min(is_decr)
   }
   i   <- nadir(pr[,'fit'])                     # index the nadir
   sek <- sek[1:i]                              # trimmed sequence
   pr  <- pr[1:i,]                              # trimmed predicted values
   p   <-  (1 - pr[,'fit'] / max(pr[,'fit']))   # percent decline from max
   CL  <- sek[which.min(abs(p - crit))]         # CL is N dep that's closest
   cx  <- coefficients(mod)                     # model coefficients
   se  <- summary(mod)$coefficients[,'Std. Error'] * 1.0 # plug-in parameter SE
   set.seed(1926)                               # sample from parameter priors
   b   <- data.frame(b0=rnorm(n, cx[1], se[1]), # intercept
                     b1=rnorm(n, cx[2], se[2]), # slope
                     b2=rnorm(n, cx[3], se[3])) # curvature
   `f` <- function(b0, b1, b2, N=sek) { b0 + b1*N + b2*N^2 } # polynomial fn
   CLs <- mapply(function(b0, b1, b2) {         # fit polynomial function
      yhat <- f(b0, b1, b2, N=sek)              # fitted values for richness
      p    <-  (1 - yhat / max(yhat))           # percent decline from max
      return(sek[which.min(abs(p - crit))])     # CL is N dep that's closest
   }, b[,1], b[,2], b[,3])                      # random parameters
   ci   <- round(quantile(CLs, probs=c(0.025,0.975)),2)
   pval <- (sum(CLs > CL) + 1) / (length(CLs) + 1)
   if(do_plot) {
      xvar <- substr(gsub('.*poly\\(', '', format(fmla)), 1,1)
      yvar <- gsub('\\ ~.*', '', format(fmla))
      plot(d[,xvar],d[,yvar],pch=16,cex=0.7,col='#00000050',ylab=yvar,xlab=xvar)
      lines(sek, f(cx[1], cx[2], cx[3]), col='cyan', lwd=2)
      abline(v=CL, col='red')
      abline(h=pr[,'fit'][which.min(abs(sek - CL))], col='red')
      points(CL, pr[,'fit'][which.min(abs(sek - CL))],pch=21,bg='gold',cex=1.5)
      ecole::add_text(0.60, 0.90, paste0('CL = ', round(CL,3)))
   }
   return(list(stats=c(CL_fitted=CL, ci, pval=pval, coefs=cx), CLs=CLs))
}
ys <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano','abun_forage')
ys <- ys[c(1,2,5,4,1,3,5,4)]
xs <- c(rep('N',4), rep('S',4))
fmla_lst <- lapply(paste0(ys, '~ poly(', xs, ', 2, raw=T)'), as.formula)

### do the Monte Carlo resampling
# png('./fig/fig_99_monte_carlo_fitlines_trace.png',
#     wid=7.5, hei=3.85, uni='in', res=1080, bg='transparent')
set_par_mercury(8, mfrow=c(2,4))
mc <- lapply(fmla_lst, function(i) { monte_carlo(i, n=99999) }) # ! TIMEWARN ! ! !
# dev.off()

### summary table
(tab <- round(sapply(mc, `[[`, 1), 3))
colnames(tab) <- paste0(xs, '_vs_', ys)
ci_rng <- tab[3,] - tab[2,] ### RANGE of UNCERTAINTY to report
# 95% variability band range as a percentage relative to each point estimate:
vb_rng_pct <- round((tab[3,] - tab[2,]) / tab[1,] * 100, 1)
(tab <- rbind(tab, ci_rng, vb_rng_pct))
write.csv(tab, file='./fig/tab_02_monte_carlo_output.csv')

### boxplot of Monte Carlo CLs, grouped by CL model
s       <- stack(data.frame(sapply(mc, `[[`, 2)))
s$ind   <- factor(s$ind, labels=LETTERS[1:length(mc)])
`bxplt` <- function(CEX = 0.7, ...) {
   plot(s$ind, s$values, outcex=0.4, boxwex=0.3, boxfill='#c1c1c1', outpch=16,
        outcol=NA, whisklty=1, staplewex=0, xlab='CL Model', xaxs='i',
        pty='s', box.lty=1, mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L',
        cex=CEX, cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
}
# png('./fig/fig_02_bxplt_monte_carlo.png',
#     wid=4, hei=4, units='in', bg='transparent', res=1080)
tiff('./fig/fig_02_bxplt_monte_carlo.tif',
     wid=4, hei=4, units='in', bg='transparent', res=700, compr='lzw')
set_par_mercury(1, mar=c(3,4,0.5,0.5), oma=c(0.1,0.1,0,0))
ylab <- expression(Randomization~CLs~(kg~ha^-1~y^-1))
bxplt(ylab=ylab, ylim=c(0,9.0), medlwd=1)
yys <- rep(c('Total spp. richness','Sensitive spp. richness',
             'Cyanolichen abundance','Forage lichen abundance'), 2)
legend('topleft', paste0(unique(s$ind), ' = ', paste0(xs, ' vs ', yys)),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
points(1:8, tab[1,], pch=23, bg='white', cex=0.5) # fitted midpoint values
# cl_2019 <- c(3.5, 3.1, 1.9, 1.3, 6.0, 2.5, 2.6, 2.3) # add CLs from 2019 ms?
# points(1:8, cl_2019, pch=23, bg='red', cex=0.6)      # add CLs from 2019 ms?
# points(1:8, tab[3,], pch=23, bg='blue', cex=0.4)     # add upr CI?
# points(1:8, tab[2,], pch=23, bg='blue', cex=0.4)     # add lwr CI?
dev.off()



### --- Fig. 1 --- one Monte Carlo run, for explanatory purposes
`monte_carlo_plot` <- function(fmla,n=999,do_plot=T,...) {
   cat('now doing', format(fmla), '...\n')
   xvar   <- substr(gsub('.*poly\\(', '', format(fmla)), 1,1)
   sek    <- seq(0.1,  max(d[,xvar], na.rm=T), by=0.01)
   mod    <- quantreg::rq(fmla,tau=0.90,data=d) # focal model
   crit   <- 0.20                               # critical decline = 20%
   newdat <- data.frame(sek)                    # sequence
   colnames(newdat) <- if(grepl('poly\\(S', paste0(fmla)[3])) 'S' else 'N'
   pr  <- predict(mod,newdat,type='none',interval='conf') # *predicted* vals
   `nadir` <- function(x) { # trim to parabola nadir (forbid increasing curve)
      is_decr <- sapply(2:length(x), function(i) (x[i] < x[i-1]) * 1)
      if(all(is_decr == 1)) length(x) else which.min(is_decr)
   }
   i   <- nadir(pr[,'fit'])                     # index the nadir
   sek <- sek[1:i]                              # trimmed sequence
   pr  <- pr[1:i,]                              # trimmed predicted values
   p   <-  (1 - pr[,'fit'] / max(pr[,'fit']))   # percent decline from max
   CL  <- sek[which.min(abs(p - crit))]         # CL is N dep that's closest
   cx  <- coefficients(mod)                     # model coefficients
   se  <- summary(mod)$coefficients[,'Std. Error'] * 1.0 # plug-in parameter SE
   set.seed(1926)                               # sample from parameter priors
   b   <- data.frame(b0=rnorm(n, cx[1], se[1]), # intercept
                     b1=rnorm(n, cx[2], se[2]), # slope
                     b2=rnorm(n, cx[3], se[3])) # curvature
   `f` <- function(b0, b1, b2, N=sek) { b0 + b1*N + b2*N^2 } # polynomial fn
   # setup plot, add original CL
   xvar <- substr(gsub('.*poly\\(', '', format(fmla)), 1,1)
   yvar <- gsub('\\ ~.*', '', format(fmla))
   plot(d[,xvar], d[,yvar], pch=NA, xlim=c(0,20), ylab='Species richness',
        xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})))
   # iterate for random parameter values
   CLs <- mapply(function(b0, b1, b2) {        # fit polynomial function
      yhat <- f(b0, b1, b2, N=sek)              # fitted values for richness
      p    <-  (1 - yhat / max(yhat))           # percent decline from max
      lines(sek, yhat, col='#00000010')         ###  <<----- simulated CLs
      return(sek[which.min(abs(p - crit))])     # CL is N dep that's closest
   }, b[,1], b[,2], b[,3])                     # random parameters
   ci   <- round(quantile(CLs, probs=c(0.025,0.975)),2)
   pval <- (sum(CLs > CL) + 1) / (length(CLs) + 1)
   # plotting the lines
   lines(sek, f(cx[1], cx[2], cx[3]), col='cyan', lwd=2) # original CL
   abline(v=CLs, col='#00000010') # simulated CLs
   abline(v=CL,  col='cyan', lwd=2)      # original CL
   # output
   return(list(stats=c(CL_fitted=CL, ci, pval=pval, coefs=cx), CLs=CLs))
}
# png('./fig/fig_01_monte_carlo_explanatory.png',
#     wid=6.5, hei=3.5, uni='in', res=1080, bg='transparent')
tiff('./fig/fig_01_monte_carlo_explanatory.tif',
     wid=6.5, hei=3.5, uni='in', bg='transparent', res=700, compr='lzw')
set_par_mercury(2, CEX=0.8)
# monte carlo fitlines
mc <- monte_carlo_plot(spprich_n_sens ~ poly(N, 2, raw = T), n=999) # ! TIMEWARN ! ! !
add_label('A')
# monte carlo histogram
hist(mc$CLs, breaks=seq(0,10,by=0.03), main='', prob=F, ylab='Frequency',
     xlab=expression(Critical ~ Loads ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     xlim=c(2,4.2), col='grey', yaxs='i') ; box(bty='L')
add_label('B')
dev.off()



####    END    ####
