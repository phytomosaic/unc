# unc
Uncertainty for Critical Loads

### Motivation
Lichen-based critical loads (CLs) are threshold deposition levels at which ecological harm is averted.  CLs (which are absolute values) must be accompanied by uncertainty estimates to avoid a false sense of precision.  Here, we quantify and map uncertainty in lichen-based CLs and exceedances for nitrogen and sulfur based on 9,000+ systematic, whole-community lichen surveys across the United States.

### Contributors
Rob Smith  
Tim Ohlert  
Linda Geiser  

### Comments and feedback
robert.smith3@usda.gov

### Usage

##### Install required packages
```
install.packages('remotes')
install.packages('crs')
install.packages('quantreg')
install.packages('rgdal')
install.packages('raster')
install.packages('labdsv')
install.packages('scales')
remotes::install_github('phytomosaic/ecole') # for helper functions
```

##### Create directories
```
pth <- '/path/to/your/project/'
setwd(pth)
dir.create('./data_raw') # manually place the two raw CSV files here
dir.create('./data')     # for the processed data file used in most analyses
dir.create('./fig')      # for figures as output
```

##### Begin to run scripts....
```
source('./R/00_data_processing.R') # creates file './data/d.rda'
```
##### ....then step thru remaining scripts in `./R`

