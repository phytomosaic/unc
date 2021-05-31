######################################################################
#
#  kriging uncertainty for NCLAS tool
#
#    Rob Smith, robert.smith3@usda.gov, 30 May 2021
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(ecole)
require(sp)
require(rgdal)
require(gstat)
require(raster)
rasterOptions(progress='text', timer=TRUE, memfrac=0.9) # memory options

### load data
dn <- read.csv(file='./krig/n_for_nclas.csv')
ds <- read.csv(file='./krig/s_for_nclas.csv')
u  <- viridis::viridis(99)

### remove a missouri outlier
dn <- dn[-which(dn$lon > -95 & dn$lon < -90 & dn$lat < 37),]
ds <- ds[-which(ds$lon > -95 & ds$lon < -90 & ds$lat < 37),]

### soften extreme values for sulfur
ds$ci_rng[which(ds$ci_rng>10)] <- 10 + rnorm(sum(ds$ci_rng>10), 0, 0.1)

# template raster is the CMAQ grid
r             <- raster('~/prj/vuln_gis/gis/ndep_cmaq.tif')
extent(r)     <- extent(c(xmin(r), xmax(r), ymin(r), ymax(r)) / 1000)
projection(r) <- gsub('units=m', 'units=km', projection(r))
trg_prj       <- projection(r) # set target projection for xy points
values(r)     <- NA            # empty values
# r<-raster(extent(xyn)*1.2,nrows=2000,ncols=800,crs=trg_prj) # optional low-rez

# create SpatialPointsDataFrame and reproject (need units=km)
`make_xy` <- function(xy, crs, ...) {    # reproject points
  xy <- as.data.frame(xy)
  xy <- xy[!(is.na(xy$lon) | is.na(xy$lat)),]  # rm NAs
  coordinates(xy) <- ~lon+lat
  proj4string(xy) <- CRS('+init=epsg:4269')    # NAD83 dec deg (FIA)
  return(spTransform(xy, crs))
}
xyn <- make_xy(dn, crs=trg_prj)
xys <- make_xy(ds, crs=trg_prj)


###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
### ! ! ! TIMEWARN ! ! ! ###   ###   ###   ###   ###   ###   ###   ###   ###

### kriging from Spherical model, max dist = 250 km, nmax = 10 neighbors

# N uncertainty [using ordinary kriging]
# model       psill    range
# 1   Nug 0.7753313   0.0000
# 2   Sph 2.2654571 987.9233
gs <- gstat(formula=ci_rng~1, locations=xyn, maxdist=1500)
v  <- variogram(gs, cutoff=1000, width=10)
(fve <- fit.variogram(v, vgm(psill=5, model='Sph', range=1000, nugget=5)))
k  <- gstat(formula=ci_rng~1, locations=xyn, model=fve, nmax=10, maxdist=250)
raster::interpolate(r, k, filename='./krig/n_unc_krig.tif',
                    format='GTiff', options=c('COMPRESS=LZW'), overwrite=T)

# S uncertainty [using ordinary kriging]
#   model     psill    range
# 1   Nug 0.9730974    0.000
# 2   Sph 2.5730151 1078.286
gs <- gstat(formula=ci_rng~1, locations=xys, maxdist=1500)
v  <- variogram(gs, cutoff=1000, width=10)
(fve <- fit.variogram(v, vgm(psill=5, model='Sph', range=1000, nugget=5)))
k  <- gstat(formula=ci_rng~1, locations=xys, model=fve, nmax=10, maxdist=250)
raster::interpolate(r, k, filename='./krig/s_unc_krig.tif',
                    format='GTiff', options=c('COMPRESS=LZW'), overwrite=T)
###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


### reset units back to m
rn <- raster('./krig/n_unc_krig.tif')
rs <- raster('./krig/s_unc_krig.tif')
extent(rn)     <- extent(c(xmin(rn), xmax(rn), ymin(rn), ymax(rn)) * 1000)
projection(rn) <- gsub('units=km', 'units=m', projection(rn))
extent(rs)     <- extent(c(xmin(rs), xmax(rs), ymin(rs), ymax(rs)) * 1000)
projection(rs) <- gsub('units=km', 'units=m', projection(rs))
trg_prj <- projection(rn)



#  ! ! ! TIMEWARN ! ! !  26 min per masking operation  ###   ###   ###   ###
### clip to CONUS borders:
# # download state bounds from gadm website:
# shapefile(getData(name='GADM',download=T,country='USA',level=1),
#    filename='~/gis/boundaries/USA_adm1.shp' )
us <- shapefile('~/gis/boundaries/USA_adm1.shp')
us <- us[!(us$NAME_1 %in% c('Alaska','Hawaii')),]
us <- spTransform(us, trg_prj)
identical(projection(us), projection(rn))
mask(rn, us, filename='./krig/n_unc_krig_clip.tif',
     format='GTiff', options=c('COMPRESS=LZW'), overwrite=T) # ! ! ! TIMEWARN ! ! !
mask(rs, us, filename='./krig/s_unc_krig_clip.tif',
     format='GTiff', options=c('COMPRESS=LZW'), overwrite=T) # ! ! ! TIMEWARN ! ! !
###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


### Fig. 1 --- thumbnail figures - UNMASKED krig
png('./krig/fig_01_unc_krig.png', wid=10, hei=4.75, unit='in',
    res=1080, bg='transparent')
set_par_mercury(2, pty='m')
plot(raster('./krig/n_unc_krig.tif'), col=u, main='N uncertainty', axes=F)
plot(raster('./krig/s_unc_krig.tif'), col=u, main='S uncertainty', axes=F)
dev.off()

### Fig. 2 --- thumbnail figures - MASKED krig
png('./krig/fig_02_unc_krig_clip.png', wid=10, hei=4.75, unit='in',
    res=1080, bg='transparent')
ecole::set_par_mercury(2, pty='m')
plot(raster('./krig/n_unc_krig_clip.tif'), col=u, main='N uncertainty', axes=F)
plot(us,add=T)
points(coordinates(xyn)*1000,col=1,cex=0.3)
plot(raster('./krig/s_unc_krig_clip.tif'), col=u, main='S uncertainty', axes=F)
plot(us,add=T)
points(coordinates(xys)*1000,col=1,cex=0.3)
dev.off()

cat(paste0('\n -----------------\n\n','Done! At: ',
           Sys.time(),'\n\n -----------------\n'))

# ### mask by tree cover... no, handle this at NCLAS instead
# # back-convert to m from km
# kr             <- raster('./krig/n_unc_krig.tif')
# extent(kr)     <- extent(c(xmin(kr), xmax(kr), ymin(kr), ymax(kr)) * 1000)
# projection(kr) <- gsub('units=km', 'units=m', projection(kr))
# msk            <- raster('~/gis/mask_lichen/forest_mask_final.tif') # mask
# msk            <- raster::projectRaster(from = msk, to = r)
# msk            <- crop(msk, extent(kr))
# kr             <- crop(kr, extent(msk))
# krm <- mask(kr, msk, filename='./krig/n_unc_krig_masked.tif',
#             format='GTiff', options=c('COMPRESS=LZW'), overwrite=T)
# plot(raster('./krig/n_unc_krig_masked.tif'), col=u)

####    END    ####
