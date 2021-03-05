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

# ### NB: Peter originally used N from calibration eqn, here we use CMAQ instead
# ecole::set_par_mercury(3)
# plot(d$cmaq_n_3yroll, d$n_lich_kghay) ; abline(0,1)
# hist(d$cmaq_n_3yroll, breaks=seq(0,30,by=0.5))
# hist(d$n_lich_kghay, breaks=seq(0,30,by=0.5))

### rename columns
onm <- c('spprich_epimac',
         'spprich_oligo', 'spprich_s_sens', 'abun_cyano', 'abun_forage',
         'cmaq_n_3yroll', 'cmaq_s_3yroll')
nnm <- c('spp_rich',
         'spprich_n_sens', 'spprich_s_sens', 'abun_cyano', 'abun_forage',
         'N', 'S')
names(d)[match(onm,names(d))] <- nnm
nnm <- c(nnm, 'meanmaxaugt5y_c','meanmindect5y_c','meanppt5y_cm',
         'ubc_td','ubc_cmd','lat','lon') # keep climate columns
range(d$lat) # 30.06993 63.80750

### fill NA to 0 in diversity measures
j <- c('spp_rich','spprich_n_sens','spprich_s_sens','abun_cyano','abun_forage')
d[,j] <- sapply(d[,j], function(x) { x[is.na(x)] <- 0 ; return(x) })

### remove nonstandard and species-poor plots
d <- d[d$plottype == 'Standard',]
d <- d[d$fia_prot==1,]
d <- d[,nnm]              # keep only used columns
d <- na.omit(d)           # keep only complete rows
d <- d[d$spp_rich > 4, ]  # drop species-poor plots ! ! ! CAUTION ! ! !
range(d$lat) # 30.71974 49.19305

### fit 90th quantile regression models
fmla <- as.formula('spp_rich ~ poly(N, 2, raw=T)') # formula
m    <- quantreg::rq(fmla, tau=0.90, data=d)       # focal model
mval <- max(m$fitted.values)          # maximum species richness
v    <- c(0, 0.10, 0.20, 0.50, 0.80)  # 0,5,10,20,50,80% thresholds
sr_decline <- mval - (mval * v)       # incremental richness decline
names(sr_decline) <- v                # clean up names
ndep <- seq(0.1, 30, by=0.01)           # sequence of possible N dep
ci <- predict(m, data.frame(N=ndep), type='none', interval='conf') # *frequentist* CI
# set_par_mercury(1)  ;  plot(d$N, d$spp_rich)  ;  points(d$N, fitted(m), col=3)
# ndep where species decline matches fit/lwr/upr (to nearest 0.01 kghay)
CL  <- sapply(sr_decline, function(x) ndep[which.min(abs(c(ci[,'fit']) - x))])
lwr <- sapply(sr_decline, function(x) ndep[which.min(abs(c(ci[,'lower']) - x))])
upr <- sapply(sr_decline, function(x) ndep[which.min(abs(c(ci[,'higher']) - x))])
data.frame(perc_decline=v, sr_decline, CL, lwr, upr) # <-- final CRITICAL LOADS

### uncertainty of CLs - Monte Carlo resampling of parametric coefficients
crit <- 0.20               # critical species decline = 20%
ndep <- seq(0,25,by=0.05)  # vector of possible N dep
(cx  <- coefficients(m))   # OLD --> 32.08192918  -1.19453724  0.01107993
(se  <- summary(m)$coefficients[,'Std. Error']) * 0.50 # plug-in SE of parameters?
n    <- 999                # n replicates
set.seed(1926)             # sample from parameter prior distributions = b
# b <- data.frame(b0=rnorm(n, cx[1], 5),
#                 b1=rnorm(n, cx[2], 0.1),
#                 b2=rnorm(n, cx[3], 0.01))
b <- data.frame(b0=rnorm(n, cx[1], se[1]),
                b1=rnorm(n, cx[2], se[2]),
                b2=rnorm(n, cx[3], se[3]))
`f` <- function(b0, b1, b2, N=ndep) { b0 + b1*N + b2*N^2 } # polynomial function
CLs <- mapply(function(b0, b1, b2) {   # fit polynomial w random parameters
   yhat <- f(b0, b1, b2, N=ndep)       # fitted values for richness
   mval <- max(yhat)                   # maximum species richness
   v    <- crit                        # 20% change at CL
   p    <-  (1 - yhat / mval)          # percent decline along N dep
   CL   <- ndep[which.min(abs(p - v))] # CL is N dep that's closest
   return(CL) }, b[,1], b[,2], b[,3])
# histogram and Monte Carlo p-value
hist(CLs, breaks=seq(0,99,by=0.1), main='', prob=F, ylab='Frequency',
     xlab=expression(Critical ~ Loads ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     xlim=c(0,25), col='grey', yaxs='i') ; box(bty='L')
abline(v=CL[names(CL)==crit], col='red', lwd=3)
pval  <- paste0('p = ', (sum(CLs > CL[names(CL)==crit]) + 1) / (length(CLs) + 1))
seval <- paste0('SE = ', round(sd(CLs),2)) # SEM = SD of Monte Carlo CLs
(cirng <- round(quantile(CLs, probs=c(0.025,0.975)),2))
cival <- paste0('CI = ', cirng[1], ' - ', cirng[2]) #
# mnval <- paste0('Mean = ', round(mean(CLs),2)) # mean of Monte Carlo CLs
add_text(0.7, 0.85, labels=pval)  # p-value
# add_text(0.7, 0.78, labels=mnval) # mean
add_text(0.7, 0.78, labels=seval) # SE-value
add_text(0.7, 0.71, labels=cival) # SE-value
# dev.off()


### TODO: repeat for all CLs...


### TODO: break out exceedances by region






####    END    ####
