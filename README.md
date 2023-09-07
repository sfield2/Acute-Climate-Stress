# Acute-Climate-Stress
All scripts relating to Field et al. 2023 JQS Submission "Climate and Community in Central Mesa Verde".
This directory contains code for: 
1. Creating acute stress regimes for ancestral Pueblo villages in the Central Mesa Verde region (or any region located in the northern US Southwest).
2. Identifying periods of acute stress within acute stress regimes.
3. Applying Granger Causality tests to demographic time series and acute stress regimes.

# Data
Portions of the data used in this research is provided infolder entitled "DATA". Sensitive site location data is not included.

# Analysis 
Script used for analysis is an R script, entitled "ANALYSIS.R", that contains all necessary commands for the creation of acute stress regimes and Grange causality tests. 

See [Glowacki and Field (2023)](https://github.com/sfield2/NSJ-MV-CeramicSeriation) for methods related to ceramic-based occupation assignments.

# Functions
Paleocar function used for paleoclimate reconstructions and maize niche models. Originally built by [Bocinsky (2015)](https://github.com/bocinsky/paleocar), updated by [Reese (2020)](https://github.com/kmreese-io/Reese_2020-JASR/tree/master/FUNCTIONS), and updated again for this analysis. 



If you have questions, email Sean Field (Sean.Field@uwyo.edu).

*NOTE: Several issues related to the upcoming deprecation of the SP, rgdal, rgeos packages have been addressed. Full integration of SF and Terra packages will be completed as [r-spatial](https://r-spatial.org/) updates continue. 
