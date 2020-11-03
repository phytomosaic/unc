######################################################################
#
#  uncertainty -- data processing
#
#    Rob Smith, phytomosaic@gmail.com, 17 Aug 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
# remotes::install_github('phytomosaic/ecole') # <<----- if needed
require(ecole)
require(crs)

### load and clean data
d <- read.csv('./data_raw/MegaDbPLOT_2020.05.09.csv', header=T, stringsAsFactors=F)
names(d) <- ecole::clean_text(names(d), lower=T)
d <- d[!(duplicated(d$megadbid)),]      # (!) one megadbid was duplicated
rownames(d) <- d$megadbid
d <- d[!(duplicated(d$plot)),]          # keep 1 instance of ea plot location
d[,c('latusedd','longusedd')] <- NULL   # kill bad columns
names(d)[names(d) %in% c('latusenad83','longusenad83')] <- c('lat','lon')
d <- d[,c('lon','lat', names(d)[!names(d) %in% c('lon','lat')])] # reorder cols
d <- d[!(is.na(d$lon) | is.na(d$lat)),] # keep only non-NA coords
d <- d[!(d$lat>49.5 & d$lat < 50),]     # rm one invalid location in Canada
d <- d[!(d$lon > (-75) & d$lat < 38),]  # rm one invalid location in Atlantic
is_ak <- d$lat > 50    # assign Alaska data where none exist
d$cmaq_n_3yroll[is_ak] <- d$n_lich_kghay[is_ak]
d$cmaq_n_3yroll[is_ak & is.na(d$cmaq_n_3yroll)] <- 1.5
d$cmaq_s_3yroll[is_ak] <- 1.5
d$ubc_mwmt[d$ubc_mwmt < (-999)] <- NA
d$meanmaxaugt5y_c[is_ak] <- d$ubc_mwmt[is_ak]
d$ubc_map[d$ubc_map < 0]  <- NA
d$meanppt5y_cm[is_ak]    <- d$ubc_map[is_ak]


##   detrending climate  #########################################
x <- d                              # temporary data.frame
x <- x[!is.na(x$n_airscore),]       # remove if lacking airscore
x <- x[!is.na(x$meanmaxaugt5y_c),]  # remove if lacking climate data
x <- x[!is.na(x$meanppt5y_cm),]     # remove if lacking climate data
x <- x[!is.na(x$cmaq_n_3yroll ),]   # remove if lacking CMAQ data
set.seed(88) ; x$meanmaxaugt5y_c <- jitter(x$meanmaxaugt5y_c,fac=2.5)
### bin 2 climate variables into quantile groups (E vs W separately)
k <- 5 # number of climate groups (5 = quintiles)
`bin` <- function(x) {
  x2  <- as.numeric(ecole::cut_quantile(x$meanppt5y_cm, n=k))
  x3  <- as.numeric(ecole::cut_quantile(x$meanmaxaugt5y_c, n=k))
  return(interaction(x2, x3, sep='_'))
}
iswest <- which(x$us_loc == 'west')
iseast <- which(x$us_loc == 'east')
x$grp  <- factor(
  NA, levels=interaction(rep(1:k,ea=k),rep(1:k,time=k),sep='_'))
x$grp[iseast] <- bin(x[iseast,])
x$grp[iswest] <- bin(x[iswest,])
table(x$grp[iseast])
table(x$grp[iswest])
### CRS model, using smooth for each discrete climate 'grp'
`f` <- function(x, yvar='n_airscore', xvar='cmaq_n_3yroll', ...) {
  # x <- x[iseast,]
  x$y   <- x[,yvar]
  x$x1  <- x[,xvar]
  m <- crs::crs(y ~ x1 + grp, degree=4, segments=1, lambda=0.1,
                cv='none', kernel=TRUE, basis='tensor', data=x)
  f <- fitted(m) # fitted raw values
  # get offset coefficients for groups
  o <- plot(m, mean=T, ci=T, plot.behavior='data')[[2]]
  # offset by group mean (plus 'intercept')
  a <- f - o$mean[match(x$grp, o$grp)] + mean(o$mean)
  return(list(m = m, fit = f, adj = a))
}
ne <- f(x[iseast,], 'n_airscore', 'cmaq_n_3yroll')
nw <- f(x[iswest,], 'n_airscore', 'cmaq_n_3yroll')
se <- f(x[iseast,], 's_airscore', 'cmaq_s_3yroll')
sw <- f(x[iswest,], 's_airscore', 'cmaq_s_3yroll')
x$nfit <- x$nadj <- x$sfit <- x$sadj <- NA
x$nfit[iswest]  <- nw$fit
x$nadj[iswest]  <- nw$adj
x$nfit[iseast]  <- ne$fit
x$nadj[iseast]  <- ne$adj
x$sfit[iswest]  <- sw$fit
x$sadj[iswest]  <- sw$adj
x$sfit[iseast]  <- se$fit
x$sadj[iseast]  <- se$adj
### match to original full dataframe
d$nadj <- x$nadj[match(d$megadbid, x$megadbid)]
d$sadj <- x$sadj[match(d$megadbid, x$megadbid)]
d$nfit <- x$nfit[match(d$megadbid, x$megadbid)]
d$sfit <- x$sfit[match(d$megadbid, x$megadbid)]
rm(x,iswest,iseast,nw,ne,sw,se,f)
###   END detrending climate  #####################################



###   calc EXCEEDANCES   ##########################################
# lookup table CLs from Diversity paper (are invariant to geography)
tab <- data.frame(metric = rep(c('Species richness',
                                 'Sensitive species richness',
                                 'Forage lichen abundance',
                                 'Cyanolichen abundance',
                                 'Comm composition shifts'),len=10),
                  element = rep(c('N','S'),ea=5),
                  deposition = c(3.5, 3.1, 1.9, 1.3, 1.5,
                                 6.0, 2.5, 2.6, 2.3, 2.7))
lab <- c('ex_sr','ex_sr_sens','ex_sr_forag','ex_sr_cyano','ex_compn')
### based on CMAQ
n_exceed <- sweep(x = d[,rep('cmaq_n_3yroll',ea=5)],
                  MARGIN = 2, STATS  = tab$deposition[1:5])
s_exceed <- sweep(x = d[,rep('cmaq_s_3yroll',ea=5)],
                  MARGIN = 2, STATS  = tab$deposition[6:10])
names(n_exceed) <- paste0('n_cmaq_', lab)
names(s_exceed) <- paste0('s_cmaq_', lab)
d <- cbind(d, n_exceed, s_exceed)  ; rm(n_exceed, s_exceed)
### based on raw airscores
n_exceed <- sweep(x = d[,rep('n_airscore',ea=5)],
                  MARGIN = 2, STATS  = tab$deposition[1:5])
s_exceed <- sweep(x = d[,rep('s_airscore',ea=5)],
                  MARGIN = 2, STATS  = tab$deposition[6:10])
names(n_exceed) <- paste0('n_raw_', lab)
names(s_exceed) <- paste0('s_raw_', lab)
d <- cbind(d, n_exceed, s_exceed)  ; rm(n_exceed, s_exceed)
### based on adjusted airscores
n_exceed <- sweep(x = d[,rep('nadj',ea=5)],
                  MARGIN = 2, STATS  = tab$deposition[1:5])
s_exceed <- sweep(x = d[,rep('sadj',ea=5)],
                  MARGIN = 2, STATS  = tab$deposition[6:10])
names(n_exceed) <- paste0('n_adj_', lab)
names(s_exceed) <- paste0('s_adj_', lab)
d <- cbind(d, n_exceed, s_exceed)  ; rm(n_exceed, s_exceed)
###   END exceedances   ###########################################


### sort by coordinates
d <- d[order(d$lat, d$lon),]


### save
save(d, file='./data/d.rda')







####    END    ####
