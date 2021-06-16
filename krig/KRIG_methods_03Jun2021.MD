### Description of methods for lichen exceedance uncertainty
 
31 May 2021

# Community airscores CL model

Geiser et al. (2021, in press) calculated a community composition "airscore" for each site as the community-weighted mean of species' deposition optima. Optima are the deposition value at which each species reached its peak frequency from 9,000+ nationwide sites.  Therefore, airscores depict air quality as the average central tendency of deposition tolerances among all species in a community, weighted by their abundances.  Using nonparametric multiplicative regression (McCune 2006) of airscores in response to N or S deposition, Geiser et al. (2021, in press) defined this CL as the uppermost deposition value at which the slope of the regression surface did not significantly depart from zero, and reported CLs at 1.5 kg N ha−1 y−1 and 2.7 kg S ha−1 y−1.  Each CL is interpreted as the deposition value below which no detrimental change in lichen community compositions occurs.  Exceedances were calculated for each site as the deviation of observed CMAQ deposition from the CL (exceedance = observed – CL).

# Uncertainty of exceedances

We defined exceedance uncertainty as the breadth of a 95% bootstrap variability band (2.5th to 97.5th percentiles) per site.  Bootstrap resampling evaluates uncertainty owing to the identity of species that could have been obtained with random probability equal to observed frequencies (i.e., the effects of measurement error).  We resampled species with replacement within sites while keeping richness the same as observed, and calculated the resulting airscore and CL exceedance each iteration.  This process was iterated 999 times per site (999 replicates × 8,875 eligible sites = 8,866,125 bootstrap replicates for N or S separately), giving a 95% bootstrap variability band interpreted as uncertainty in the presence of measurement error.  

# Interpolated uncertainty maps

To map interpolated values of uncertainty, we used ordinary kriging.  First, we fit an empirical variogram model (using the observed values at measured plots) to parameterize a final variogram model based on constant spherical covariance describing the spatial decay of similarity among sites (R functions `gstat::variogram`  and `gstat::vgm` ).  From the variogram model, we applied ordinary kriging (R function `gstat::krige` ) to predict interpolated values for unsampled locations between sites (maximum 10 neighboring sites, maximum 250-km neighbor distance).

# Uncertainty, reclassified as 1-5 scores

To be consistent with other ecosystem components, we converted the full-breadth confidence interval to values on a 5-point scale.  Cutpoints were: < 1 kg/ha/y = 5, 
1-2 kg/ha/y = 4, 2-3 kg/ha/y = 3, 3-4 kg/ha/y = 2, >4 kg/ha/y = 1.  The lowest uncertainty (highest scores) tended to occur in the Pacific coastal states, the intermountain Northwest, and the far Northeast (for both nitrogen and sulfur), in addition to the Southeast (for nitrogen) and Southwest (for sulfur).

# References

Geiser, L.H., H.T. Root, R.J. Smith, S. Jovan, L. St. Clair, K.L. Dillman and D. Schwede. 2021. Lichen-based critical loads for atmospheric deposition in US Forests. *Global Change Biology*. In review.

McCune, B. 2006. Non-parametric habitat models with automatic interactions. Journal of Vegetation Science 17:819–830.


####    END    ####
