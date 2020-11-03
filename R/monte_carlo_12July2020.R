######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- Monte Carlo parameter resampling
#
#    Rob Smith, phytomosaic@gmail.com, 12 July 2020
#
##      GNU General Public License, Version 3.0    ###################


### Root/Nelson CL models -- Lichen species richness vs Nitrogen

### preamble
rm(list=ls())
require(ecole)
require(quantreg)

### load and clean data
# dd <- read.csv(url(paste0('https://raw.githubusercontent.com/nelsopet/',
#                          'Lichen_CL_quantile_regression/master/output/',
#                          'LichDb_sFncGrpSens_All_Abun.csv')))
# m <- dd$megadbid  ;  rm(dd)# keep exact same plots Peter used???

### load MEGADB data
d <- read.csv('./data_raw/MegaDbPLOT_2020.05.09.csv', stringsAsFactors=F)

# ### Peter's used N from calibration equation
# set_par_mercury(6)
# plot(d$cmaq_n_3yroll, d$n_lich_kghay)
# hist(d$cmaq_n_3yroll, breaks=seq(0,30,by=0.5))
# hist(d$n_lich_kghay, breaks=seq(0,30,by=0.5))
# plot(d$latusedd, d$n_lich_kghay)
# plot(d$latusedd, d$cmaq_n_3yroll)
# # hist(d$n_lich_kghay, breaks=seq(0,30,by=0.5))

### CMAQ doesnt exist for Alaska, substitute N-calibration estimates
is_ak <- !is.na(d$latusedd) & (d$latusedd > 49.00001)
d$cmaq_n_3yroll[is_ak] <- d$n_lich_kghay[is_ak]
rm(is_ak)

### rename columnns
onm <- c('spprich_epimac','cmaq_n_3yroll','meanmaxaugt5y_C',
         'meanmindect5y_c','meanppt5y_cm','ubc_td','ubc_cmd')
nnm <- c('spp_rich','N','maxaug_c','mindec_c',
         'precip_cm','continen','CMD')
names(d)[match(onm,names(d))] <- nnm
rm(onm, nnm)

# ### OPTION 1 -- try matching Peter's plots
# d <- d[d$megadbid %in% m,]  ;  rm(m)

### OPTION 2 -- remove nonstandard plots
d <- d[d$plottype == 'Standard',]
d <- d[d$fia_prot==1,]
j <- c('spp_rich','N','maxaug_c','mindec_c','precip_cm','continen','CMD')
d <- d[,j]       # keep only used columns
d <- na.omit(d)  # only complete rows
round(cor(d), 2) # correlations

### drop species-poor plots ! ! ! CAUTION ! ! !
d <- d[d$spp_rich > 4, ]

### list of models with different combinations of predictors
XVAR <- c('N',
          'poly(N, 2, raw=T)',
          'N + maxaug_c + mindec_c + precip_cm + CMD',
          'poly(N, 2, raw=T) + maxaug_c + mindec_c + precip_cm + CMD')
XVAR_names <- c('N', 'N_poly', 'N_clim', 'N_poly_clim')

### fit 90th quantile regression models
fmla_lst <- lapply(paste0('spp_rich ~ ', XVAR), as.formula)
mod_lst  <- lapply(fmla_lst, function(x) rq(x, tau=0.90, data=d))
names(mod_lst) <- XVAR_names
# null model
mod_null <- rq(formula=spp_rich~1, tau=0.90, data=d, model=T)
# R1 goodness-of-fit, but the intercept is hard coded
mod_R1 <- sapply(1:length(XVAR),
                 function(y) 1 - mod_lst[[y]]$rho / mod_null$rho)
# AIC for each model and combine R1 statistic
(mod_compare <- data.frame(nm  = XVAR_names,
                           AIC = sapply(mod_lst, AIC),
                           R1  = mod_R1,
                           fmla= paste(fmla_lst)))
# write.csv(mod_compare,'./N_mod_compare.csv')

### fitted values and confidence intervals
`pr` <- function (x) {
   as.data.frame(predict(x, d, interval='confidence', level=0.95,
                         # type='percentile', se='boot', N=9999,
                         type='none', se='nid'))
}
mod_pred <- as.data.frame(mapply(pr,mod_lst)) # ! ! ! TIME WARN ! ! !

# ### percent decline in lichen SR for some select N dep increments
# m <- mod_lst$N_poly      # selected model
# mfit <- m$fitted.values  # selected model fitted values
# ndep <- data.frame(N=c(1,1.5,2,2.5,3,5,7.5,10,12.5,15,17.5,20))
# predict(m, ndep)                   # absolute decline
# 100 * (1 - (predict(m, ndep)) / mval) # percent decline
# range(d$N, na.rm=T)

### calc Critical Loads (CLs) at 0,5,10,20,50,80% changes
m    <- mod_lst$N_poly                # focal model
mval <- max(m$fitted.values)          # maximum species richness
v    <- c(0, 0.10, 0.20, 0.50, 0.80)  # 0,5,10,20,50,80% changes
sr_decline <- mval - (mval * v)       # incremental richness decline
names(sr_decline) <- v                # clean up names
ndep <- seq(0.1,30,by=0.01)           # sequence of possible N dep
# fit and CI across WHOLE range of N (can bootstrap CIs instead)
ci <- data.frame(predict(m, data.frame(N=ndep),
                         type='none',
                         interval='confidence',
                         level=0.95))
# index NDEP where species decline (nearly) matches fitted/lwr/upr
CL  <- sapply(sr_decline,
              function(x) ndep[which.min(abs(c(ci[,'fit']) - x))])
lwr <- sapply(sr_decline,
              function(x) ndep[which.min(abs(c(ci[,'lower']) - x))])
upr <- sapply(sr_decline,
              function(x) ndep[which.min(abs(c(ci[,'higher']) - x))])
### final CRITICAL LOADS
data.frame(perc_decline=v, sr_decline, CL, lwr, upr)

### plot
tiff('C:/Users/Rob/Desktop/CL.tiff', wid=4.5, hei=4.5, uni='in',
     bg='transparent', res=700)
set_par_mercury(1)
plot(d$N, d$spp_rich, cex=0.7, xlim=c(0,28), # col='#00000050',
     col='#00000000',
     xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     ylab='Species richness')
# fit line
points(m$x[,2], m$fitted.values, col=4, cex=0.2, pch=16)
# CL points
points(CL, sr_decline, col=2, cex=1.1, pch=16)
# text for percent decline
text(CL, sr_decline + 3, labels=v, col=2, cex=0.8, pch=16, font=2)
# droplines for CLs
arrows(CL, sr_decline, CL, rep(0,5), col=2, angle=12, len=0.1)
# CIs for CLs
segments(CL, sr_decline, lwr, sr_decline, col=4)
segments(CL, sr_decline, upr, sr_decline, col=5)
# arrows(CL, sr_decline, lwr, sr_decline, col=4, angle=12, len=0.1)
# arrows(CL, sr_decline, upr, sr_decline, col=5, angle=12, len=0.1)
dev.off()


### uncertainty of CLs - Monte Carlo resampling of parametric coefficients
(cx <- coefficients(m)) # 32.08192918  -1.19453724  0.01107993
ndep <- seq(0,25,by=0.05)         # vector of N dep
`fn` <- function(b0, b1, b2, N=ndep) {
   b0 + b1*N + b2*N^2 # polynomial function
}
summary(m)
(se <- summary(m)$coefficients[,'Std. Error']) # SE of parameters - plug-in?
n <- 999
b <- data.frame(b0=rnorm(n, cx[1], 5),
                b1=rnorm(n, cx[2], 0.1),
                b2=rnorm(n, cx[3], 0.01))
tiff('C:/Users/Rob/Desktop/CL_monte_carlo.tiff', wid=8.75, hei=4.5, uni='in',
     bg='transparent', res=700)
set_par_mercury(2)
plot(d$N, d$spp_rich, xlim=c(0,28), type='n',
     xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     ylab='Species richness')
ndep <- seq(0,25,by=0.05)           # vector of N dep
CLs <- mapply(function(b0,b1,b2) {
   ndep <- seq(0,25,by=0.05)        # vector of N dep
   yhat <- fn(b0, b1, b2, N=ndep)   # fitted values for richness
   mval <- max(yhat)                # maximum species richness
   v    <- 0.10                     # 10% change at CL
   p    <-  (1 - yhat / mval)       # percent decline along N dep
   CL   <- ndep[which.min(abs(p - v))]  # CL is N dep that's closest
   lines(ndep, yhat, col='#00000040')
   return(CL) }, b[,1], b[,2], b[,3])
lines(ndep, fn(32.0819, -1.19453724, 0.01107993), col=2, lwd=3)
abline(v=CLs, col='#00000005')
abline(v=CL[names(CL)==0.1], col='red', lwd=3)
### histogram and Monte Carlo p-value
hist(CLs, breaks=seq(0,99,by=0.1), main='', prob=F, ylab='Frequency',
     xlab=expression(Critical ~ Loads ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
     xlim=c(0,6), col='grey', yaxs='i') ; box(bty='L')
# lines(density(CLs), col=4)
# plot(ecdf(CLs), vert=T, do.points=F, add=T, col=4, col.01line=NULL)
abline(v=CL[names(CL)==0.1], col='red', lwd=3)
pval <- paste0('p = ', (sum(CLs > CL[names(CL)==0.1]) + 1) / (length(CLs) + 1))
seval <- paste0('SE = ', round(sd(CLs),2)) # SEM = SD of Monte Carlo CLs
# mnval <- paste0('Mean = ', round(mean(CLs),2)) # mean of Monte Carlo CLs
add_text(0.7, 0.85, labels=pval)  # p-value
# add_text(0.7, 0.78, labels=mnval) # mean
add_text(0.7, 0.78, labels=seval) # SE-value
dev.off()




#' # upper/lower bounds for CI of incremental % decline
#' # 0 % 32.79772 30.98566 34.62121
#' # 5 % 31.14400 29.78724 32.55391
#' # 10% 29.46393 28.26562 30.44257
#' # 20% 26.21598 25.56544 26.81200
#' # 50% 21.68773 20.96845 22.41112
#' # 80% 17.08599 16.36759 17.73913
#'
#' # Use same heuristic solving approach as above to find deposition
#' # associated with upper and lower confidence intervals of incremental % decline
#' ## Heuristically solved for X given a Y ...	34.58681
#' ## 0 % (max spp_rich) Heuristically solved for X given a Y 10.41...   undef-0.08-0.94  N kg/ha/yr
#' ## 5 % loss of N. Cya Heuristically solved for X given a Y 9.89...	  0.19 -0.86 -1.54 N kg/ha/yr
#' ## 10% loss of N. Cya Heuristically solved for X given a Y 9.37...	  1.21 -1.70 -2.34 N kg/ha/yr
#' ## 20% loss of N. Cya Heuristically solved for X given a Y 8.324995...3.15 -3.50 -3.9 N kg/ha/yr
#' ## 50% loss of N. Cya Heuristically solved for X given a Y 5.203122...n/a  N kg/ha/yr
#' ## 80% loss of N. Cya Heuristically solved for X given a Y 2.081249...n/a kg/ha/yr
#'
#' ### plots with fitted lines, CIs and percent decline point estimates
#' # Get coeffecients to plot in figures
#' round(coefficients(m),2)
#' #32.97                -2.19                 0.07
#' # richness = 32.97 - 2.19*Ndep + 0.07*Ndep^2
#' #
#' # tiff(paste('SR_vs_N_percentiles','_',
#' #            format(Sys.time(),'%y%d%m'),'.tiff',sep=''),
#' #      wid=1200, hei=1200, units='px', pointsize=24)
#' names(mod_pred)
#' x         <- unlist(mod_pred[2], recursive=F, use.names=T)
#' clim      <- unlist(mod_pred[4],recursive=F, use.names=T)
#' clim_pred <- unlist(clim[1], recursive=F, use.names=T)
#'
#' fit  <- unlist(x[1], recursive = F, use.names =T)
#' low  <- unlist(x[2], recursive = F, use.names =T)
#' high <- unlist(x[3], recursive = F, use.names =T)
#' tmp  <- cbind(d, fit, low, high, clim_pred)  ### <--- need `match`
#' tmp  <- subset(tmp, N < 12.5)
#' plot(d$spp_rich ~ d$N,
#'      xlab=expression(Nitrogen ~ (kg ~ N ~ ha^{-1} ~ y^{-1})),
#'      ylab='Lichen Species Richness',
#'      # col=c('darkorange','darkgreen')[Area],
#'      cex.main = 1.6, cex.lab  = 1.7, cex.axis = 1.7)
#' points(tmp$N, tmp$clim_pred, col=alpha('grey',0.25), pch= 19)
#' points(tmp$N, tmp$fit, col='blue', pch= 19)
#' points(tmp$N, tmp$low, col='black', pch= 4, cex=0.25)
#' points(tmp$N, tmp$high, col='black', pch= 4, cex=0.25)
#' #Red dot coordinates are determined by heuristically solving for
#' #incremental percent decline (right value) to determine assoc
#' #  deposition (left value)
#' points(0.08,32.80, col='red',     pch=20,cex=2) #0% decline
#' #points(0.86,31.16,col='red',     pch=20,cex=2) #5% decline
#' points(1.70,29.52, col='red',     pch=20,cex=2) #10% decline
#' points(3.50,26.24, col='red',     pch=20,cex=2) #20% decline
#' points(6.65,5.203122, col='red', pch=20,cex=2) #50% decline
#' points(12.8,2.081249, col='red', pch=20,cex=2) #80% decline
#' text(max(d$N)*0.6,max(d$spp_rich)*0.65,
#'      labels=bquote(atop(.('Lichen Species richness = 32.97 -'),
#'                         .('2.19*N + 0.07*')*N^2)), cex=1.5)
#' #text(max(d$N)*0.6,max(d$spp_rich)*0.9,
#' #labels=bquote(atop(.('Lichen Species richness = 32.97 -'),.('2.19*N + 0.07*')*N^2)), cex=1.5)
#' legend('topright', legend=c('Raw Data East','Raw Data West', 'Fitted Values Deposition only', '95% Bootstrapped Confidence Interval Deposition only','Fitted Values Climate + Deposition','0/10/20/50/80% loss'), col=c('darkorange','darkgreen','blue', 'black','grey','red'), pch=16, cex=1.47)
#' #legend('topright', legend=c('Raw Data East','Raw Data West',
#' #'Fitted Values Deposition only', '95% Bootstrapped CI Deposition only',
#' #''Fitted Values Climate + Deposition','0/10/20/50/80% loss'),
#' #'col=c('darkorange','darkgreen','blue', 'black','grey','red'), pch=16, cex=1)
#' dev.off()
#'
#' ####    END    ####
