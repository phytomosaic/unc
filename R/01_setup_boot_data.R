######################################################################
#
#  CLAD WG-2 -- setup for ALL bootstraps
#
#    Rob Smith, phytomosaic@gmail.com, 16 Feb 2021
#
##      GNU General Public License, Version 3.0    ###################

rm(list=ls())
require(ecole)
require(labdsv)

# removes if:
#   species lacks a rating
#   site has < 4 rated species
#   site lacks CMAQ data
#   site is apparently outside USA

### define functions
`make_ratings_matrix` <- function(which_element='n', ...) {
  s  <- read.csv('./data_raw/MegaDbLICHEN_2020.05.09.csv',
                 stringsAsFactors=F, fileEncoding = 'latin1')
  names(s) <- ecole::clean_text(names(s), lower=T)
  s  <- s[s$growthform == 'Epiphytic macrolichen',] # keep only macrolichens
  s  <- s[,c('megadbid','sci_22chklst','fia_abun',
             'n_depmaxfreq','s_depmaxfreq')]
  s$sci_22chklst <- ecole::clean_text(s$sci_22chklst, TRUE)
  # get species ratings ('peak detection frequency' for N or S)
  #   this comes from HTR's 2018 spline regressions
  wa <- s[!duplicated(s$sci_22chklst),
          c('sci_22chklst','n_depmaxfreq','s_depmaxfreq')]
  wa <- wa[order(wa$sci_22chklst),]
  # reshape wide for abundance matrix
  spe <- labdsv::matrify(s[,c('megadbid','sci_22chklst','fia_abun')])
  # convert spe to trait (ratings to be sampled)
  x <- spe                                   # duplicate
  x[x > 0 ] <- 1                             # convert to binary
  # select which element
  if (which_element == 'n') {
    x <- sweep(x, 2, wa$n_depmaxfreq, '*')   # nitrogen ratings
  } else if (which_element == 's') {
    x <- sweep(x, 2, wa$s_depmaxfreq, '*')   # sulfur ratings
  }
  # filter out weak species
  x[is.na(x)] <- 0L    # zero out any NA (species that lacked ratings)
  j   <- colSums(x)==0 # 287 spp lacking valid spp ratings
  spe <- spe[,!j]
  x   <-   x[,!j]
  i   <- apply(x, 1, function(x) sum(x>0)) < 4 # species-poor (<4 rated species)
  spe <- spe[!i,]
  x   <-   x[!i,]
  x   <- as.matrix(x)
  return(list(spe = spe, x = x))
}
`match_site_data` <- function(x, spe, which_element='n', ...) {
  d <- read.csv('./data_raw/MegaDbPLOT_2020.05.09.csv', stringsAsFactors=F)
  names(d) <- ecole::clean_text(names(d), lower=T)
  d <- d[!(duplicated(d$megadbid)),]      # (!) one megadbid was duplicated
  d$megadbid  <- as.character(d$megadbid)
  rownames(d) <- d$megadbid
  d[,c('latusedd','longusedd')] <- NULL   # kill bad columns
  names(d)[names(d) %in% c('latusenad83','longusenad83')] <- c('lat','lon')
  d <- d[,c('lon','lat',names(d)[!names(d) %in% c('lon','lat')])] # reorder cols
  d <- d[!(is.na(d$lon) | is.na(d$lat)),] # keep only non-NA coords
  d <- d[!(d$lat > 49.5 & d$lat < 50),]   # rm one invalid location in Canada
  d <- d[!(d$lon > (-75) & d$lat < 38),]  # rm one invalid location in Atlantic
  is_ak <- d$lat > 50                     # assign Alaska data where none exist
  d$cmaq_n_3yroll[is_ak] <- d$n_lich_kghay[is_ak]
  d$cmaq_n_3yroll[is_ak & is.na(d$cmaq_n_3yroll)] <- 1.5
  d$cmaq_s_3yroll[is_ak] <- 1.5
  d$ubc_mwmt[d$ubc_mwmt < (-999)] <- NA
  d$meanmaxaugt5y_c[is_ak] <- d$ubc_mwmt[is_ak]
  d$ubc_map[d$ubc_map < 0] <- NA
  d$meanppt5y_cm[is_ak]    <- d$ubc_map[is_ak]
  inm <- intersect(dimnames(x)[[1]], d$megadbid)
  spe <- spe[dimnames(spe)[[1]] %in% inm, ]
  x   <- x[dimnames(x)[[1]] %in% inm, ]
  d   <- d[d$megadbid %in% inm, ]
  spe <- spe[order(dimnames(spe)[[1]]), ]
  x   <- x[order(dimnames(x)[[1]]),]
  d   <- d[order(d$megadbid),]
  # select the element
  if (which_element == 'n') {
    i   <- which(!is.na(d$cmaq_n_3yroll)) # omit if lacking CMAQ values
  } else if (which_element == 's') {
    i   <- which(!is.na(d$cmaq_s_3yroll)) # omit if lacking CMAQ values
  }
  spe <- spe[i, ] # omit if lacking CMAQ values
  x   <- x[i,]    # omit if lacking CMAQ values
  d   <- d[i,]    # omit if lacking CMAQ values
  # check that rownames and dimensions all equal
  stopifnot({
    identical(dimnames(spe)[[1]], d$megadbid) &
      identical(dimnames(x)[[1]], d$megadbid) &
      all.equal(dim(spe), dim(x)) &
      all.equal(NROW(spe), NROW(d))
  })
  # observed site values: airscore and richness
  d$scr_obs <- apply(x, 1, function(i) mean(i[i>0])) # obs airscore (this way)
  d$sr      <- apply(x, 1, function(i) sum(i>0))     # obs richness of RATED spp
  return(list(spe = spe, x = x, d = d))
}

cat(paste0('start time: ', time_start <- Sys.time()), '\n')

### setup for N ratings
o   <- make_ratings_matrix('n')
o   <- match_site_data(o$x, o$spe, 'n')
x   <- o$x
d   <- o$d
spe <- o$spe
save(d,   file='./res/d_n_scr.rda')   # save descriptor matrix
save(x,   file='./res/x_n_scr.rda')   # save ratings matrix
save(spe, file='./res/spe_n_scr.rda') # save species abundance matrix

### setup for S ratings
o   <- make_ratings_matrix('s')
o   <- match_site_data(o$x, o$spe, 's')
x   <- o$x
d   <- o$d
spe <- o$spe
save(d,   file='./res/d_s_scr.rda')   # save descriptor matrix
save(x,   file='./res/x_s_scr.rda')   # save ratings matrix
save(spe, file='./res/spe_s_scr.rda') # save species abundance matrix

cat(paste0('time elapsed: ', Sys.time()-time_start), '\n') # 12 sec to process

rm(list=ls()) # cleanup

####    END    ####
