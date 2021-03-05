######################################################################
#
#  CLAD WG-2 -- parallel bootstrap -- species *richness*
#     bootstrap species that enter the lichen airscore, varying richness
#     answers sampling sufficiency question
#
#    Rob Smith, phytomosaic@gmail.com, 16 Feb 2021
#
##      GNU General Public License, Version 3.0    ###################

rm(list=ls())
require(parallel)

### RICHNESS function
`boot_airscore` <- function(i, B=999, n) {
  xx <-  x[i,]      # subset row of species ratings (traits)
  pp <- pr[i,]      # subset row of species occ probabilities
  a  <- replicate(B, mean(sample(xx, size=n, replace=T, prob=pp))) # airscores
  return(c(quantile(a, c(0.025,0.50,0.975)), sem = sd(a)))
}

### ! ! ! TIMEWARN - Nitrogen - 7.9 min, 11 nodes, B=999, 8875 sites, 346 spp
load(file='./res/x_n_scr.rda')   # ratings matrix
load(file='./res/spe_n_scr.rda') # species abundance matrix
pr  <- as.matrix(sweep(spe, 1, rowSums(spe), '/')) # species site occ probs
# ii <- 1:99      ### TESTING Chop down
# x  <-  x[ii,]
# pr <- pr[ii,]
cat(paste0('start time: ', time_start <- Sys.time()), '\n')
for (N in 2:20) {
  cat(paste0('\n\t\t< < < < Size ',N,' of 20 @ ',Sys.time(),' > > > > \n'))
  b <- t(parallel::mcmapply(FUN = boot_airscore,
                            1:NROW(pr),
                            MoreArgs = list(n=N),
                            mc.cores = 11))
  b <- cbind(b, b[,3] - b[,1])                                # breadth of CIs
  dimnames(b)[[1]] <- dimnames(pr)[[1]]                       # rownames
  dimnames(b)[[2]] <- c('lwr', 'med', 'upr', 'sem', 'ci_rng') # colnames
  save(b, file=paste0('./res/rich_boot/b_',N,'_n_scr.rda'))
  rm(b)
}
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
rm(spe, pr)


### ! ! ! TIMEWARN - Sulfur - 8.1 min, 11 nodes, B=999, 8918 sites, 324 spp
load(file='./res/x_s_scr.rda')   # ratings matrix
load(file='./res/spe_s_scr.rda') # species abundance matrix
pr  <- as.matrix(sweep(spe, 1, rowSums(spe), '/')) # species site occ probs
cat(paste0('start time: ', time_start <- Sys.time()), '\n')
for (N in 2:20) {
  cat(paste0('\n\t\t< < < < Size ',N,' of 20 @ ',Sys.time(),' > > > > \n'))
  b <- t(parallel::mcmapply(FUN = boot_airscore,
                            1:NROW(pr),
                            MoreArgs = list(n=N),
                            mc.cores = 11))
  b <- cbind(b, b[,3] - b[,1])                                # breadth of CIs
  dimnames(b)[[1]] <- dimnames(pr)[[1]]                       # rownames
  dimnames(b)[[2]] <- c('lwr', 'med', 'upr', 'sem', 'ci_rng') # colnames
  save(b, file=paste0('./res/rich_boot/b_',N,'_s_scr.rda'))
  rm(b)
}
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
rm(spe, pr)


####    END    ####
