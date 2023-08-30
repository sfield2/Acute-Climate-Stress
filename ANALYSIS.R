############## DEFINING ACUTE CLIMATE STRESS ################

      ## RELEASE 1.0 NOTES ###

### Current release does not directly rely on 'sp' package (retired Oct 2023). However,  
### current release still uses 'raster' packages, which depend on 'sp'. 'Raster' 
### not be deprecated but should be fully replaced at some point by 'terra' package. 
### As 'terra' replaces raster, the paleocar_update function and the following script 
### will also be updated.

### Step 1-4 include processes needed for identification of acute stress regimes.
### Step 5 includes processes for granger causality tests.

      ## END RELEASE 1.0 NOTES ###


######### STEP 0.0: # ENVIRONMENTS & DATA ######################################

packages <-c('tidyverse','sf','terra','ggplot2','magrittr','plyr','FedData','raster',
             'dplyr','leastcostpath','gdistance', 'ggpubr','zoo','lmtest',
             'tseries','car','stats')

for(p in packages) if(p %in% rownames(installed.packages()) == F) { install.packages(p) }
for(p in packages) suppressPackageStartupMessages(library(p,quietly=T,character.only=T))

setwd("C:/Users/seanp/Documents/PROJECT/AcuteClimateStress")
theme_set(theme_bw())

######### STEP 1.0: IMPORT FUNCTIONS & DATA #####################################

base::source('./FUNCTIONS/paleocar_update.R')

### Import raster data for defining study area
dem <- raster('./DATA/VEPIIN_NED_1.tif',na.rm=T)

### import location of village and boundary of surrounding land form
fv <- st_read("./DATA/farview.shp")%>%
  vect()

mvl <- st_read("./DATA/Mesa_verde_landform.shp")%>%
  vect()

### re-project data to ensure accurate
dem <- raster::projectRaster(dem,crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
fv <- project(fv,"+proj=longlat +datum=WGS84 +ellps=WGS84")
mvl <- project(mvl,"+proj=longlat +datum=WGS84 +ellps=WGS84")

### check data locations
plot(dem)
points(fv,col="black")
lines(mvl,col="blue")

######### STEP 2.0: DEFINE STUDAY AREA WITH LEAST COST ANALYSIS ################
### Define areas within 2 hr round-trip travel from center of village
### First, set limits for computational savings 
neigh <- 16
dem.aggregate <- 5
dem.aggregated <- aggregate(dem, fact=dem.aggregate, fun=mean)

### Build cost surface using Tobler's hiking function
tobler_cs <- leastcostpath::create_slope_cs(dem = dem.aggregated, cost_function = "tobler", neighbours = neigh)

### Derive location data for center of Far View Community
geom_fv<-as.matrix(geom(fv))%>%
  .[1,3:4]

### Build accumulated cost surface (cost accumulated from center of Far View Community)
acc_fv_cs<- accCost(tobler_cs,geom_fv)

### Define catchment by limiting accumulated cost surface to one-hour single-direction travel
### Transform into polygon and determine extent of polygon. 
### NOTE: Catchment or extent can be used as the paleocar template input
limit_fv<- acc_fv_cs<3600

fv_catchment<-rasterToContour(limit_fv)%>%
  .[1,]%>%
  st_as_sf()%>%
  st_cast("POLYGON")%>%
  sf::st_set_crs(.,4326)%>%
  as_Spatial()

fv_extent <- extent(fv_catchment)%>%
  as(.,"SpatialPolygons")

### Check the limits of catchment and extent
plot(fv_catchment,col="white",add=T)
plot(fv_extent,col="grey70",add=T)

### Export catchment areas example
fv_catchment_sf <- st_as_sf(fv_catchment)
st_write(fv_catchment_sf, "./fv_catchment.shp")

fv_extent_sf <- st_as_sf(fv_extent)
st_write(fv_extent_sf, "./fv_extent.shp")

######### STEP 3.0: CLIMATE AND MAIZE NICHE RECON ##############################
### Build reconstruction from local catchment extent
growing_niche_fv <- paleocar(template = fv_extent,
                          label = 'paleoproductivity',
                          raw.dir = paste0('./output/FV_EXTENT/spatial-products/paleocar/raw/'),
                          extraction.dir = paste0('./output/FV_EXTENT/spatial-products/paleocar/extractions/'),
                          prcp_threshold = 350,
                          gdd_threshold = 1800,
                          years = 1:2000,
                          force.redo = F)

### Find average size of niche, PPT, and GDD from reconstruction
niche_fv <- growing_niche_fv$niche%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(sum = rowSums(across(where(is.numeric))))%>%
  mutate(pct_in_niche = (sum/110)*100)%>%
  mutate(year = 600:1300)

precip_fv <- growing_niche_fv$PPT%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_precip = rowMeans(across(where(is.numeric))))%>%
  mutate(average_precip_cm = average_precip/10)%>%
  mutate(year = 600:1300)

gdd_fv <- growing_niche_fv$GDD%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_gdd = rowMeans(across(where(is.numeric))))%>%
  mutate(year = 600:1300)

### Build data frame just for local reconstruction
fv_recon <- as.data.frame(cbind(niche_fv$year,niche_fv$pct_in_niche,precip_fv$average_precip_cm,gdd_fv$average_gdd))
colnames(fv_recon) <- c("Year", "Percent_in_niche", "Average_precipitation_in_cm", "Average_growing_degree_days")

### Build reconstruction from local catchment extent
mvl_extent<-st_as_sf(mvl)%>%
  extent()%>%
  as(.,"SpatialPolygons")

### Check the limits of extent
plot(mvl_extent,add=T,col="green")

### Build reconstruction from non-local extent
growing_niche_mvl <- paleocar(template = mvl_extent,
                             label = 'paleoproductivity',
                             raw.dir = paste0('./output/MESAVERDELANDFORM/spatial-products/paleocar/raw/'),
                             extraction.dir = paste0('./output/MESAVERDELANDFORM/spatial-products/paleocar/extractions/'),
                             prcp_threshold = 350,
                             gdd_threshold = 1800,
                             years = 1:2000,
                             force.redo = F)

niche_mvl <- growing_niche_mvl$niche%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(sum = rowSums(across(where(is.numeric))))%>%
  mutate(pct_in_niche = (sum/980)*100)%>%
  mutate(year = 600:1300)

precip_mvl <- growing_niche_mvl$PPT%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_precip = rowMeans(across(where(is.numeric))))%>%
  mutate(average_precip_cm = average_precip/10)%>%
  mutate(year = 600:1300)

gdd_mvl <- growing_niche_mvl$GDD%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_gdd = rowMeans(across(where(is.numeric))))%>%
  mutate(year = 600:1300)

### Build data frame for non-local reconstruction
mvl_recon <- as.data.frame(cbind(niche_mvl$year,niche_mvl$pct_in_niche,precip_mvl$average_precip_cm,gdd_mvl$average_gdd))
colnames(mvl_recon) <- c("Year", "Percent_in_niche", "Average_precipitation_in_cm", "Average_growing_degree_days")


######### STEP 4.0: ANALYZE CLIMATE RECONSTRUCTIONS AND   ######################
######### IDENTIFY ACUTE STRESS PERIODS    #####################################
### Run following calculations:
### 4.1 calculate generational averages in both conditions
### 4.2 calculate percent difference between generational av and annual conditions in each year
### 4.3 transform into quartiles
### 4.4 incorporate percentage of region in niche

fv_variab<- fv_recon %>%
  # 4.1
  mutate(Gen_Precipitation = rollmean(Average_precipitation_in_cm,20,fill=NA,align="right"))%>%
  mutate(Gen_Temperature = rollmean(Average_growing_degree_days,20,fill=NA,align="right"))%>%
  # 4.2
  mutate(Precipitation_diff = abs(1 - (Average_precipitation_in_cm/Gen_Precipitation)))%>%
  mutate(Temperature_diff = abs(1 - (Average_growing_degree_days/Gen_Temperature)))%>%
  # 4.3
  mutate(Precipitation_diff_quartile = ntile(Precipitation_diff,4))%>%
  mutate(Temperature_diff_quartile = ntile(Temperature_diff,4))%>%
  # 4.4 
  mutate(Region_low_niche = mvl_recon$Percent_in_niche)

### Tally years when any of the four measures are present:
### Measure 1) High variab in precip
### Measure 2) High variab in temp
### Measure 3) Low local niche size
### Measure 4) Low local and non-local niche size 

fv_stress <- list(fv_variab$Year[fv_variab$Precipitation_diff_quartile == "4"],
             fv_variab$Year[fv_variab$Temperature_diff_quartile == "4"],
             fv_variab$Year[fv_variab$Percent_in_niche < 50],
             fv_variab$Year[fv_variab$Region_low_niche < 50 & fv_variab$Percent_in_niche <50])

### Count how many measures in each year
fv_stress <- as.data.frame(unlist(fv_stress))
stress_ct <- count(fv_stress$`unlist(fv_stress)`)

### Build stress count into time series (i.e., stress regime)
fv_stress_regime <- fv_variab$Year%>%
  data.frame(Year=.)%>%
  mutate(stress = stress_ct$freq[match(Year,stress_ct$x)])%>%
  mutate_all(funs(ifelse(is.na(.), 0, .)))%>%
  mutate(gen_av = rollmean(stress, 9, align="right", fill=NA))

### Identify periods with 5 or more consecutive years where more than one 
### stress measure is present 
potential_stress_periods <- subset(fv_stress_regime, gen_av > "1")

sp <- NULL 
sapply(seq(nrow(potential_stress_periods)), function(x)
{
  ifelse((sum(diff(potential_stress_periods[x:(x+4), "Year"], 1)) == 4 &
            sum(diff(potential_stress_periods[x:(x+4), "Year"], 1) == 1) == 4),
         sp <<- rbind(sp, potential_stress_periods[x:(x+4),]),"")
})

fv_stress_periods <- unique(sp)


######### STEP 5.0: APPLY GRANGER CAUSALITY TEST  ##############################
### Import demographic estimates
pop <- read.csv("./DATA/all_community_occupations.csv",header=T)

### Import Far View Stress regime if skipping to this step from Step 0
fv_stress_regime <- read.csv("./DATA/FVC_stress_regime.csv",header=T)

### Smooth demographic curve for Far View by minimum phase lenght 
### (i.e., minimum number of years of occupation phase based on Glowacki 
### and Field 2023)
pop$smoothed <- ceiling(rollapply(pop$fv_pop_av,40,mean,fill=NA,align="center",partial=T))

### Limit Far View demographic curve to study period 
fv_stress_pop <- subset(fv_stress_regime,Year >=700 & Year <= 1280)

### Test if Far View stress regime is stationary
### Interpretation: if p-value is smaller than 0.05, we reject null hypothesis,
### indicating that the time series is stationary
adf.test(fv_stress_pop$gen_av)

### Test if demographic curve is stationary
adf.test(pop$smoothed)

### Because demographic curve is not stationary, apply time lags and re-run analysis 
### Do this for time lags of 1 year to 25 years
granger_tests<- as.data.frame(matrix(NA,25,3))
colnames(granger_tests) <- c("ADF p-value","Granger F-value"," Granger p-value")

for(i in 1:25){
  pop$smoothed_diff <- pop$smoothed-lag(pop$smoothed,i)
  pop$smoothed_diff[is.na(pop$smoothed_diff)]<-0
  adf_t<-adf.test(pop$smoothed_diff)
  granger_tests[i,1]<- adf_t$p.value
  
  gt<-grangertest(pop$smoothed_diff~fv_stress_pop$gen_av, order=1)
  
  granger_tests[i,2]<- gt[2,3]
  granger_tests[i,3]<- gt[2,4]
}


