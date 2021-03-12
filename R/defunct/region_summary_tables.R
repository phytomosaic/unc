library(dplyr)
library(plyr)

rm(list=ls())

### load bootstrapped 95% CI for site airscores (each row = one site)
load(file='~/res/b_n_scr.rda', verbose=T)
bn  <- t(b) # bootstrap CIs for N
load(file='~/res/b_s_scr.rda', verbose=T)
bs  <- t(b) # bootstrap CIs for S
rm(b)

### load matching site-descriptor matrices
load(file='~/res/d_n_scr.rda', verbose=T)
dn <- cbind(d,bn)
load(file='~/res/d_s_scr.rda', verbose=T)
ds <- cbind(d,bs)
rm(d,bs,bn)

### remove degenerate solutions
i  <- is.na(dn$upr)
sum(i) #  37 degenerate solutions for N
dn <- dn[!i,]
i  <- is.na(ds$upr)
sum(i) # 480 degenerate solutions for S
ds <- ds[!i,]
rm(i)

### calc CI breadth = range of possible values = UNCERTAINTY
dn$ci_rng <- dn$upr - dn$lwr # calc breadth of CIs
ds$ci_rng <- ds$upr - ds$lwr # calc breadth of CIs

### calc site exceedances: bootstrapped mean minus the fixed CL
dn$exc <- dn$mean - 1.5   # N airscores CL = 1.5
ds$exc <- ds$mean - 2.7   # S airscores CL = 2.7


##Tables of useful information

#First need to sort each state into the regions being used in the manuscript
#These regions can be changed by sorting the state's into different regions

dn$unc.region <- dn$state
dn$unc.region <- dplyr::recode(dn$unc.region,
                               "Utah" = "Southwest",
                               "Arizona" = "Southwest",
                               "Colorado" = "Southwest",
                               "New Mexico" = "Southwest",
                               "Nevada" = "Southwest",
                               
                               "Alaska" = "Alaska",
                               
                               "Oregon" = "Northwest",
                               "Montana" = "Northwest",
                               "Wyoming" = "Northwest",
                               "Idaho" = "Northwest",
                               "Washington" = "Northwest",
                               
                               "California" = "West.Coast",
                               
                               "Michigan" = "Upper.Midwest",
                               "Wisconsin" = "Upper.Midwest",
                               "Minnesota" = "Upper.Midwest",
                               
                               "Indiana" = "Lower.Midwest",
                               "Ohio" = "Lower.Midwest",
                               "Illinois" = "Lower.Midwest",
                               "Missouri" = "Lower.Midwest",
                               
                               "Georgia" = "Southeast",
                               "Alabama" = "Southeast",
                               "South Carolina" = "Southeast",
                               "Tennessee" = "Southeast",
                               "Florida" = "Southeast",
                               
                               "Pennsylvania" = "Midatlantic",
                               "North Carolina" = "Midatlantic",
                               "Virginia" = "Midatlantic",
                               "West Virginia" = "Midatlantic",
                               "Maryland" = "Midatlantic",
                               "Delaware" = "Midatlantic",
                               
                               "New York" = "Northeast",
                               "Massachusetts" = "Northeast",
                               "Rhode Island" = "Northeast",
                               "Connecticut" = "Northeast",
                               "New Hampshire" = "Northeast",
                               "Vermont" = "Northeast",
                               "Maine" = "Northeast",
                               "New Jersey" = "Northeast"
)


ds$unc.region <- ds$state
ds$unc.region <- dplyr::recode(ds$unc.region,
                               "Utah" = "Southwest",
                               "Arizona" = "Southwest",
                               "Colorado" = "Southwest",
                               "New Mexico" = "Southwest",
                               "Nevada" = "Southwest",
                               
                               "Alaska" = "Alaska",
                               
                               "Oregon" = "Northwest",
                               "Montana" = "Northwest",
                               "Wyoming" = "Northwest",
                               "Idaho" = "Northwest",
                               "Washington" = "Northwest",
                               
                               "California" = "West.Coast",
                               
                               "Michigan" = "Upper.Midwest",
                               "Wisconsin" = "Upper.Midwest",
                               "Minnesota" = "Upper.Midwest",
                               
                               "Indiana" = "Lower.Midwest",
                               "Ohio" = "Lower.Midwest",
                               "Illinois" = "Lower.Midwest",
                               "Missouri" = "Lower.Midwest",
                               
                               "Georgia" = "Southeast",
                               "Alabama" = "Southeast",
                               "South Carolina" = "Southeast",
                               "Tennessee" = "Southeast",
                               "Florida" = "Southeast",
                               
                               "Pennsylvania" = "Midatlantic",
                               "North Carolina" = "Midatlantic",
                               "Virginia" = "Midatlantic",
                               "West Virginia" = "Midatlantic",
                               "Maryland" = "Midatlantic",
                               "Delaware" = "Midatlantic",
                               
                               "New York" = "Northeast",
                               "Massachusetts" = "Northeast",
                               "Rhode Island" = "Northeast",
                               "Connecticut" = "Northeast",
                               "New Hampshire" = "Northeast",
                               "Vermont" = "Northeast",
                               "Maine" = "Northeast",
                               "New Jersey" = "Northeast"
)


summ_N <- ddply(dn, .(unc.region),
                function(x)data.frame(
                  avg.exceedence = mean(x$exc),
                  avg.confint = mean(x$ci_rng),
                  min.confint = min(x$ci_rng),
                  max.confint = max(x$ci_rng),
                  n = length(x$exc)
                ))

summ_N$CI_percentofexc <- summ_N$avg.confint / summ_N$avg.exceedence



summ_S <- ddply(ds, .(unc.region),
                function(x)data.frame(
                  avg.exceedence = mean(x$exc),
                  avg.confint = mean(x$ci_rng),
                  min.confint = min(x$ci_rng),
                  max.confint = max(x$ci_rng),
                  n = length(x$exc)
                ))

summ_S$CI_percentofexc <- summ_S$avg.confint / summ_S$avg.exceedence


