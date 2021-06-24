######################################################################
#
#  CLAD WG-2 -- estimating uncertainty --
#
#   COMBINED *richness* and *composition* uncertainty
#
#    Rob Smith, phytomosaic@gmail.com, 03 Mar 2021
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(ecole)
require(scales)
require(ggplot2)
require(gridExtra)
require(usmap)
require(viridis)

### setup plot labels
nlabr <- expression(N~uncertainty~(kg~ha^-1~y^-1))
slabr <- expression(S~uncertainty~(kg~ha^-1~y^-1))
nlabc <- expression(N~uncertainty~(kg~ha^-1~y^-1))
slabc <- expression(S~uncertainty~(kg~ha^-1~y^-1))

### *composition* load bootstrap results  -- nitrogen and sulfur
load(file='./res/b_n_scr.rda')   # bootstrap CIs
bn  <- b                         # bootstrap CIs for N
load(file='./res/b_s_scr.rda')   # bootstrap CIs
bs  <- b                         # bootstrap CIs for S
load(file='./res/d_n_scr.rda')   # descriptor matrix
dn <- cbind(d,bn)
load(file='./res/d_s_scr.rda')   # descriptor matrix
ds <- cbind(d,bs)
rm(b,d,bs,bn)
### *richness* load bootstrap results -- nitrogen
fnm <- list.files('./res/rich_boot/', pattern='b_[[:digit:]]')
fnm <- fnm[grep('n_', fnm)] # grab nitrogen
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) })
ci_rng_n <- sapply(lst, '[', , 5, drop=T)  # uncertainty as CI95 range
dimnames(ci_rng_n)[[2]] <- 2:20
dn$ci_rng_n20 <- ci_rng_n[,'20']  # *richness* uncertainties
rm(lst, b, fnm, ci_rng_n)
### *richness* load bootstrap results -- sulfur
fnm <- list.files('./res/rich_boot/', pattern='b_[[:digit:]]')
fnm <- fnm[grep('s_', fnm)] # grab sulfur
fnm <- fnm[order(as.numeric(gsub('[^[:digit:]]','', fnm))-1)]
fnm <- as.list(paste0('./res/rich_boot/', fnm))
lst <- lapply(fnm, function(x) { drop(get(load(x,.GlobalEnv))) })
ci_rng_s <- sapply(lst, '[', , 5, drop=T)  # uncertainty as CI95 range
dimnames(ci_rng_s)[[2]] <- 2:20
ds$ci_rng_s20 <- ci_rng_s[,'20']  # *richness* uncertainties
rm(lst, b, fnm, ci_rng_s)




### --- NOTRUN --- DIAGNOSTICS ----------------------------------------------

###  ignore some zero values
dn$ci_rng[dn$ci_rng == 0]         <- NA   # *composition* uncertainties
ds$ci_rng[dn$ci_rng == 0]         <- NA   # *composition* uncertainties
dn$ci_rng_n20[dn$ci_rng_n20 == 0] <- NA   # *richness* uncertainties
ds$ci_rng_s20[dn$ci_rng_s20 == 0] <- NA   # *richness* uncertainties
# control extreme values, just to aid visualization
dn$ci_rng[dn$ci_rng > 11]         <- NA   # *composition* uncertainties
ds$ci_rng[ds$ci_rng > 11]         <- NA   # *composition* uncertainties
ds$ci_rng_s20[ds$ci_rng_s20 > 7]  <- NA   # *composition* uncertainties
###   a few corrections
dn$ubc_mat[dn$ubc_mat < -99] <- NA
ds$ubc_mat[ds$ubc_mat < -99] <- NA
dn$ubc_map[dn$ubc_map < 20]  <- 20
ds$ubc_map[ds$ubc_map < 20]  <- 20
set.seed(88)
dn$ubc_map <- jitter(dn$ubc_map, factor=2000)
set.seed(88)
ds$ubc_map <- jitter(ds$ubc_map, factor=2000)
dn$ubc_map <- dn$ubc_map + abs(min(dn$ubc_map, na.rm=T)) + 50
ds$ubc_map <- ds$ubc_map + abs(min(ds$ubc_map, na.rm=T)) + 50
dn$sr <- dn$spprich_epimac
ds$sr <- ds$spprich_epimac

### --- Fig. 07 --- diagnostics
png('./fig/fig_07_diagnostics.png',
    wid=5.5, hei=6.5, uni='in', res=1080, bg='transparent')
set_par_mercury(12, CEX=1) ; par(mfrow=c(4,3), oma=c(0,1.5,1,0), cex=0.35)
u <- '#00000060'
plot(jitter(dn$sr, factor=2),  dn$ci_rng_n20, col=u, xlim=c(0,50),
     xlab='Species richness', ylab=nlabr)  ;  add_label('A')
mtext('Bootstrap richness', 2, 4, outer=F, las=3, cex=0.8)
plot(dn$ubc_mat, dn$ci_rng_n20, col=u,
     xlab='Mean ann temp', ylab=nlabr, xlim=c(-3,21))  ;  add_label('B')
plot(dn$ubc_map, dn$ci_rng_n20, col=u,
     xlab='Mean ann precip', ylab=nlabr, log='x', xlim=c(70,7500))  ;  add_label('C')
# N composition
plot(jitter(dn$sr, factor=2),  dn$ci_rng, col=u, xlim=c(0,50),
     xlab='Species richness', ylab=nlabc)  ;  add_label('D')
mtext('Bootstrap identities', 2, 4, outer=F, las=3, cex=0.8)
plot(dn$ubc_mat, dn$ci_rng, col=u,
     xlab='Mean ann temp', ylab=nlabc, xlim=c(-3,21))  ;  add_label('E')
plot(dn$ubc_map, dn$ci_rng, col=u,
     xlab='Mean ann precip', ylab=nlabc, log='x', xlim=c(70,7500))  ;  add_label('F')
plot(jitter(ds$sr, factor=2),  ds$ci_rng_s20, col=u, xlim=c(0,50),
     xlab='Species richness', ylab=slabr)  ;  add_label('G')
mtext('Bootstrap richness', 2, 4, outer=F, las=3, cex=0.8)
plot(ds$ubc_mat, ds$ci_rng_s20, col=u,
     xlab='Mean ann temp', ylab=slabr, xlim=c(-3,21))  ;  add_label('H')
plot(ds$ubc_map, ds$ci_rng_s20, col=u,
     xlab='Mean ann precip', ylab=slabr, log='x', xlim=c(70,7500))  ;  add_label('I')
# S composition
plot(jitter(ds$sr, factor=2),  ds$ci_rng, col=u, xlim=c(0,50),
     xlab='Species richness', ylab=slabc)  ;  add_label('J')
mtext('Bootstrap identities', 2, 4, outer=F, las=3, cex=0.8)
plot(ds$ubc_mat, ds$ci_rng, col=u,
     xlab='Mean ann temp', ylab=slabc, xlim=c(-3,21))  ;  add_label('K')
plot(ds$ubc_map, ds$ci_rng, col=u,
     xlab='Mean ann precip', ylab=slabc, log='x', xlim=c(100,7500))  ;  add_label('L')
dev.off()

####    END    ####
