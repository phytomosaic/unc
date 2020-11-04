######################################################################
#
#  CLAD WG-2 -- estimating uncertainty -- CMAQ bootstrap resampling
#
#    Rob Smith, phytomosaic@gmail.com, 03 Nov 2020
#
##      GNU General Public License, Version 3.0    ###################



require(rgdal)
require(raster)

### load points
load('./data/d.rda')
str(d)
### load CMAQ values
x <- paste0('precip_adj_bias_adj_withWestCoast_2012_CMAQv5_0_1_',
            'bidi_withNH3MetEPIC_updates_12km_CONUS_N_S_kg_ha')
s <- rgdal::readOGR(dsn = '/home/rob/gis/cmaq_n_dep/', layer = x)
proj4string(s)

### convert points to spatial
`make_xy` <- function(xy, trg_prj, ...){
  coordinates(xy)   <- ~lon+lat
  # proj4string(xy) <- CRS('+init=epsg:4326') # WGS84
  proj4string(xy)   <- CRS('+init=epsg:4269')  # <<--- NAD83 -----
  return(spTransform(xy, projection(trg_prj)))
}
xy <- make_xy(d[c('lon','lat')], trg_prj=proj4string(s))

###
`get_cmaq` <- function(xy, field='TD_N', plotid = dimnames(d)[[1]]) {
  val  <- rep(NA, NROW(xy))
  names(val) <- dimnames(env)[[1]]
  blks <- split(1:NROW(xy), ceiling(seq_along(1:NROW(xy))/20))
  nblk <- length(blks)
  t0   <- Sys.time()
  for(b in 1:nblk){
    cat('\nblock', b, 'of', nblk)
    system.time(val[ blks[[b]] ] <-
                  raster::extract(s, xy[blks[[b]],])[[field]]
    )
  }
  cat('\ntime elapsed:', Sys.time() - t0,'\n')
  val
}


### for each point, generate B=99 random points (Gaussian distributed)
B <- 99
a <- array(NA, dim=c(NROW(d), 2, B))  # array of random coords centered on pts
o <- matrix(NA, nrow=NROW(d), ncol=B) # array of resampled outcomes


d <- d[order(d$lat, d$lon),]
# 0.1 decimal degree is about 10-11 km
for(i in 1:NROW(a)) { # for each site
  a[i,1,] <- d$lon[i] + rnorm(B, 0, 0.1)
  a[i,2,] <- d$lat[i] + rnorm(B, 0, 0.1)
}

plot(d[111:119,c('lon','lat')])

for(i in 1:NROW(a)) { # for each site
  i <- 111
  points(a[i,1,], a[i,2,], cex=0.2, col='red')
}




dim(a)
dim(b)
head(a)
class(a)
array()






dimnames(xy)
rownames(xy@data)


### extract values (kg ha-1 y-1)   ! ! ! TIMEWARN ! ! !
d$n_tot <- get_cmaq(xy, 'TD_N')     # total wet+dry dep of all N
d$s_tot <- get_cmaq(xy, 'TD_S_T')  # total wet+dry dep of all S



rm(s, d, x, xy, make_xy, get_cmaq)
