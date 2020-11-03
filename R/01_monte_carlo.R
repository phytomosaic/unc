######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- Monte Carlo parameter resampling
#
#    Rob Smith, phytomosaic@gmail.com, 12 July 2020
#
##      GNU General Public License, Version 3.0    ###################


### TODO: repeat for all CLs...


### preamble
rm(list=ls())
require(ecole)     # for plotting and name handling
require(quantreg)  # for 90th percentile regressions

### load data
load('./data/d.rda')

# ### Peter originally used N from calibration eqn, here we use CMAQ instead
# ecole::set_par_mercury(3)
# plot(d$cmaq_n_3yroll, d$n_lich_kghay)
# hist(d$cmaq_n_3yroll, breaks=seq(0,30,by=0.5))
# hist(d$n_lich_kghay, breaks=seq(0,30,by=0.5))

### rename columns
onm <- c('spprich_epimac','cmaq_n_3yroll','meanmaxaugt5y_c',
         'meanmindect5y_c','meanppt5y_cm','ubc_td','ubc_cmd')
nnm <- c('spp_rich','N','maxaug_c','mindec_c',
         'precip_cm','continen','CMD')
names(d)[match(onm,names(d))] <- nnm
rm(onm, nnm)

### -- remove nonstandard plots
d <- d[d$plottype == 'Standard',]
d <- d[d$fia_prot==1,]
j <- c('spp_rich','N','maxaug_c','mindec_c','precip_cm','continen','CMD')
d <- d[,j]                # keep only used columns
d <- na.omit(d)           # keep only complete rows
round(cor(d), 2)          # correlation matrix
d <- d[d$spp_rich > 4, ]  # drop species-poor plots ! ! ! CAUTION ! ! !

### list of four models with different combinations of predictors
XVAR <- c('N',
          'poly(N, 2, raw=T)',
          'N + maxaug_c + mindec_c + precip_cm + CMD',
          'poly(N, 2, raw=T) + maxaug_c + mindec_c + precip_cm + CMD')
XVAR_names <- c('N', 'N_poly', 'N_clim', 'N_poly_clim')
# fit 90th quantile regression models
fmla_lst <- lapply(paste0('spp_rich ~ ', XVAR), as.formula)
mod_lst  <- lapply(fmla_lst, function(x) quantreg::rq(x, tau=0.90, data=d))
names(mod_lst) <- XVAR_names
# null model
mod_null <- rq(formula=spp_rich~1, tau=0.90, data=d, model=T)
# R1 goodness-of-fit, but the intercept is hard coded
(mod_R1 <- sapply(1:length(XVAR),
                  function(y) 1 - mod_lst[[y]]$rho / mod_null$rho))
# AIC for each model and combine R1 statistic
(mod_compare <- data.frame(nm  = XVAR_names,
                           AIC = sapply(mod_lst, AIC),
                           R1  = mod_R1,
                           fmla= paste(fmla_lst)))
# write.csv(mod_compare,'./N_mod_compare.csv')

# ### percent decline in lichen SR for some select N dep increments
# m    <- mod_lst$N_poly                # selected model
# mfit <- m$fitted.values               # selected model fitted values
# mval <- max(m$fitted.values)          # maximum species richness
# ndep <- data.frame(N=c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20))
# pr   <- predict(m, ndep)              # absolute decline
# 100 * (1 - pr / mval)                 # percent decline from max value

### calc Critical Loads (CLs) at 0,5,10,20,50,80% changes
m    <- mod_lst$N_poly                # focal model
mval <- max(m$fitted.values)          # maximum species richness
v    <- c(0, 0.10, 0.20, 0.50, 0.80)  # 0,5,10,20,50,80% changes
sr_decline <- mval - (mval * v)       # incremental richness decline
names(sr_decline) <- v                # clean up names
ndep <- seq(0.1,30,by=0.01)           # sequence of possible N dep
# fit and *frequentist* CI across WHOLE range of N
ci <- data.frame(predict(m, data.frame(N=ndep),
                         type='none', interval='confidence', level=0.95))
# # fit and *bootstrap* CIs (IGNORE for now) ! ! ! TIMEWARN boot CIs ! ! !
# ci <- as.data.frame(mapply(function(x) {
#    as.data.frame(predict(x, d, interval='confidence', level=0.95,
#                          type='percentile', se='boot', N=999, se='nid'))},
#    mod_lst))
# index NDEP where species decline matches fitted/lwr/upr (to nearest 0.1 kghay)
CL  <- sapply(sr_decline,
              function(x) ndep[which.min(abs(c(ci[,'fit']) - x))])
lwr <- sapply(sr_decline,
              function(x) ndep[which.min(abs(c(ci[,'lower']) - x))])
upr <- sapply(sr_decline,
              function(x) ndep[which.min(abs(c(ci[,'higher']) - x))])
data.frame(perc_decline=v, sr_decline, CL, lwr, upr) # <--- final CRITICAL LOADS


# ### plot the regression and resulting CLs
# png('./fig/CL.png', wid=4.5, hei=4.5, uni='in', bg='transparent', res=700)
# ecole::set_par_mercury(1)
# plot(d$N, d$spp_rich, cex=0.7, xlim=c(0,28), # col='#00000050',
#      col='#00000000',
#      xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
#      ylab='Species richness')
# # fit line
# points(m$x[,2], m$fitted.values, col=4, cex=0.2, pch=16)
# # CL points
# points(CL, sr_decline, col=2, cex=1.1, pch=16)
# # text for percent decline
# text(CL, sr_decline + 3, labels=v, col=2, cex=0.8, pch=16, font=2)
# # droplines for CLs
# arrows(CL, sr_decline, CL, rep(0,5), col=2, angle=12, len=0.1)
# # CIs for CLs
# segments(CL, sr_decline, lwr, sr_decline, col=4)
# segments(CL, sr_decline, upr, sr_decline, col=5)
# dev.off()


### uncertainty of CLs - Monte Carlo resampling of parametric coefficients
(cx <- coefficients(m))    # OLD --> 32.08192918  -1.19453724  0.01107993
ndep <- seq(0,25,by=0.05)  # vector of N dep

`f` <- function(b0, b1, b2, N=ndep) { b0 + b1*N + b2*N^2 } # polynomial function
summary(m)
(se <- summary(m)$coefficients[,'Std. Error']) # SE of parameters - plug-in?
n <- 999
# set up prior sampling distributions for each parameter:
set.seed(1926)
b <- data.frame(b0=rnorm(n, cx[1], 5),
                b1=rnorm(n, cx[2], 0.1),
                b2=rnorm(n, cx[3], 0.01))
png('./fig/CL_monte_carlo.png', wid=8.75, hei=4.5, uni='in',
    bg='transparent', res=700)
ecole::set_par_mercury(2)
plot(d$N, d$spp_rich, xlim=c(0,28), type='n',
     xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     ylab='Species richness')
CLs <- mapply(function(b0,b1,b2) {
   # ndep <- seq(0,25,by=0.05)           # vector of N dep
   yhat <- f(b0, b1, b2, N=ndep)       # fitted values for richness
   mval <- max(yhat)                   # maximum species richness
   v    <- 0.10                        # 10% change at CL
   p    <-  (1 - yhat / mval)          # percent decline along N dep
   CL   <- ndep[which.min(abs(p - v))] # CL is N dep that's closest
   lines(ndep, yhat, col='#00000040')
   return(CL) }, b[,1], b[,2], b[,3])
# lines(ndep, f(32.0819, -1.19453724, 0.01107993), col=2, lwd=3)
lines(ndep, f(cx[1], cx[2], cx[3]), col=2, lwd=3) # observed regression fit
abline(v=CLs, col='#00000005')                 # resampled CLs
abline(v=CL[names(CL)==0.1], col='red', lwd=3) # observed CL
# histogram and Monte Carlo p-value
hist(CLs, breaks=seq(0,99,by=0.1), main='', prob=F, ylab='Frequency',
     xlab=expression(Critical ~ Loads ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     xlim=c(0,6), col='grey', yaxs='i') ; box(bty='L')
abline(v=CL[names(CL)==0.1], col='red', lwd=3)
pval  <- paste0('p = ', (sum(CLs > CL[names(CL)==0.1]) + 1) / (length(CLs) + 1))
seval <- paste0('SE = ', round(sd(CLs),2)) # SEM = SD of Monte Carlo CLs
# mnval <- paste0('Mean = ', round(mean(CLs),2)) # mean of Monte Carlo CLs
add_text(0.7, 0.85, labels=pval)  # p-value
# add_text(0.7, 0.78, labels=mnval) # mean
add_text(0.7, 0.78, labels=seval) # SE-value
dev.off()


### TODO: repeat for all CLs...




####    END    ####
