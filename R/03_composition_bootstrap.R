######################################################################
#
#  CLAD WG-2 -- parallel bootstrap -- species *compositions*
#     bootstrap species that enter the lichen airscore, keep richness constant
#
#    Rob Smith, phytomosaic@gmail.com, 16 Feb 2021
#
##      GNU General Public License, Version 3.0    ###################

rm(list=ls())
require(parallel)

# ### TESTING Chop down
# ii <- 1:999
# x  <- x[ii,]
# pr <- pr[ii,]

### COMPOSITION function
`boot_airscore` <- function(i, B=999) {
  xx <- x[i, ]      # subset row of species ratings (traits)
  pp <- pr[i,]      # subset row of species occ probabilities
  n  <- sum(pp > 0) # sample size = observed richness
  a  <- replicate(B, mean(sample(xx, size=n, replace=T, prob=pp))) # airscores
  return(quantile(a, c(0.025,0.50,0.975))) # bootstrap CI and median of B reps
}

### ! ! ! TIMEWARN - Nitrogen - 23.7 sec, 11 nodes, B=999, 8875 sites, 346 spp
load(file='./res/x_n_scr.rda')   # ratings matrix
load(file='./res/spe_n_scr.rda') # species abundance matrix
pr  <- as.matrix(sweep(spe, 1, rowSums(spe), '/')) # species site occ probs
dim(x)  # ratings matrix to sample, 8875 sites, 346 species
cat(paste0('start time: ', time_start <- Sys.time()), '\n')
b <- t(parallel::mcmapply(FUN = boot_airscore, 1:NROW(pr), mc.cores=11))
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
b <- cbind(b, b[,3] - b[,1])          # breadth of CIs
dimnames(b)[[1]] <- dimnames(pr)[[1]]             # rownames
dimnames(b)[[2]] <- c('lwr', 'med', 'upr', 'ci_rng') # colnames
save(b, file='./res/b_n_scr.rda')
rm(spe, pr, b)


### ! ! ! TIMEWARN - Sulfur - 23.1 sec, 11 nodes, B=999, 8918 sites, 324 spp
load(file='./res/x_s_scr.rda')   # ratings matrix
load(file='./res/spe_s_scr.rda') # species abundance matrix
pr  <- as.matrix(sweep(spe, 1, rowSums(spe), '/')) # species site occ probs
dim(x)  # ratings matrix to sample, 8918 sites, 324 species
cat(paste0('start time: ', time_start <- Sys.time()), '\n')
b <- t(parallel::mcmapply(FUN = boot_airscore, 1:NROW(pr), mc.cores=11))
cat(paste0('time elapsed: ', Sys.time()-time_start), '\n')
b <- cbind(b, b[,3] - b[,1])          # breadth of CIs
dimnames(b)[[1]] <- dimnames(pr)[[1]]             # rownames
dimnames(b)[[2]] <- c('lwr', 'med', 'upr', 'ci_rng') # colnames
save(b, file='./res/b_s_scr.rda')
rm(spe, pr, b)


####    END    ####
