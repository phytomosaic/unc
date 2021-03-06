######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- Monte Carlo parameter resampling
#
#    Rob Smith, phytomosaic@gmail.com, 12 July 2020
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(ecole)        # for plotting and name handling
require(quantreg)     # for 90th percentile regressions
load('./data/d.rda')  # site data

# ### NB: PRN originally used N from calibration eqn, here we use CMAQ instead
# hist(d$cmaq_n_3yroll, breaks=seq(0,30,by=0.5), col='#00000050')
# hist(d$n_lich_kghay, breaks=seq(0,30,by=0.5), add=T, col='#00000050')

### rename columns
onm <- c('spprich_epimac','spprich_oligo','spprich_s_sens','abun_cyano',
         'abun_forage', 'cmaq_n_3yroll', 'cmaq_s_3yroll')
nnm <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano',
         'abun_forage', 'N', 'S')
names(d)[match(onm,names(d))] <- nnm
nnm <- c(nnm, 'meanmaxaugt5y_c',
         'meanmindect5y_c',
         'meanppt5y_cm',
         'ubc_td',
         'ubc_cmd',
         'lat','lon') # keep climate columns

### fill NA to 0 in diversity measures
j <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano','abun_forage')
d[,j] <- sapply(d[,j], function(x) { x[is.na(x)] <- 0 ; return(x) })

# ### trim a few extremely high outlying deposition values
# d$N[d$N > 20] <- NA
# d$S[d$S > 45] <- NA

### remove nonstandard and species-poor plots
d <- d[d$plottype == 'Standard',]
d <- d[d$fia_prot==1,]
d <- d[,nnm]              # keep only used columns
d <- na.omit(d)           # keep only complete rows (omits AK without climate)
d <- d[d$spp_rich > 4, ]  # drop species-poor plots ! ! ! CAUTION ! ! !
range(d$lat)              # effectively CONUS
dim(d)                    # 5422 sites included

### Monte Carlo resampling of parametric coefficients
`monte_carlo` <- function(fmla, sek = seq(0.1, 30, by=0.01), n=99999, ...) {
   cat('now doing', format(fmla), '...\n')
   mod    <- quantreg::rq(fmla, tau=0.90, data=d) # focal model
   crit   <- 0.20                               # critical species decline = 20%
   newdat <- data.frame(sek)
   colnames(newdat) <- if(grepl('poly\\(S', paste0(fmla)[3])) 'S' else 'N'
   pr  <- predict(mod,newdat,type='none',interval='conf') # *predicted* vals
   p   <-  (1 - pr[,'fit'] / max(pr[,'fit']))  # percent decline from max
   CL  <- sek[which.min(abs(p - crit))]        # CL is N dep that's closest
   cx  <- coefficients(mod)  # OLD --> 32.08192918  -1.19453724  0.01107993
   se  <- summary(mod)$coefficients[,'Std. Error'] * 1.0 # plug-in parameter SE
   set.seed(1926)             # sample from parameter prior distributions = b
   b   <- data.frame(b0=rnorm(n, cx[1], se[1]), # 5
                      b1=rnorm(n, cx[2], se[2]), # 0.1
                      b2=rnorm(n, cx[3], se[3])) # 0.01
   `f` <- function(b0, b1, b2, N=sek) { b0 + b1*N + b2*N^2 } # polynomial fn
   CLs <- mapply(function(b0, b1, b2) {    # fit polynomial w random parameters
      yhat <- f(b0, b1, b2, N=sek)          # fitted values for richness
      p    <-  (1 - yhat / max(yhat))       # percent decline from max
      return(sek[which.min(abs(p - crit))]) # CL is N dep that's closest
   }, b[,1], b[,2], b[,3])
   ci   <- round(quantile(CLs, probs=c(0.025,0.975)),2)
   pval <- (sum(CLs > CL) + 1) / (length(CLs) + 1)
   return(list(stats=c(CL_fitted=CL, ci, pval=pval), CLs=CLs))
}
ys <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano','abun_forage')
ys <- ys[c(1,2,4,5,1,3,4,5)]
xs <- c(rep('N',4), rep('S',4))
fmla_lst <- lapply(paste0(ys, '~ poly(', xs, ', 2, raw=T)'), as.formula)
mc <- lapply(fmla_lst, function(i) { monte_carlo(i) })

### summary table
(tab <- sapply(mc, `[[`, 1))
colnames(tab) <- paste0(xs, ' vs ', ys)
write.csv(tab, file='./fig/tab_02_monte_carlo_output.csv')

### boxplot of Monte Carlo CLs
s     <- stack(data.frame(sapply(mc, `[[`, 2)))
s$ind <- factor(s$ind, labels=LETTERS[1:length(mc)])
`bxplt` <- function(CEX = 0.7, ...) {
   plot(s$ind, s$values, outcex=0.4, boxwex=0.3, boxfill='#c1c1c1', outpch=16,
        outcol=NA, whisklty=1, staplewex=0, xlab='CL Model', xaxs='i',
        pty='s', box.lty=1, mgp=c(CEX+1.4,0.4,0), tcl=-0.2, las=1, bty='L',
        cex=CEX, cex.lab=CEX*1.4, cex.axis=CEX*1.1, cex.main=CEX*2.4, ...)
}
png('./fig/fig_02_bxplt_monte_carlo.png',
    wid=4, hei=4, units='in', bg='transparent', res=700)
set_par_mercury(1, mar=c(3,4,0.5,0.5), oma=c(0.1,0.1,0,0))
ylab <- expression(Randomization~CLs~(kg~N~ha^-1~y^-1))
bxplt(ylab=ylab, ylim=c(0,12.75), medlwd=1)
yys <- rep(c('Total spp. richness','Sensitive spp. richness',
             'Cyanolichen abundance','Forage lichen abundance'), 2)
legend('topleft', paste0(unique(s$ind), ' = ', paste0(xs, ' vs ', yys)),
       col='transparent', border=NA, bty='n', cex=0.5, ncol=2)
points(1:8, tab[1,], pch=23, bg='white', cex=0.5)
dev.off()

####    END    ####
