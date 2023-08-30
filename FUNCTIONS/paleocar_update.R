## PALEOCAR ##
## R. KYLE BOCINSKY

## FIRST MODIFIED BY KELSEY M. REESE (2020) FOR USE, ORIGINAL CODE IS COMMENTED OUT




        ### NOTES FOR KELSEY ### 


## AGAIN MODIFIED BY S. FIELD (2023) TO REMOVE EXPLICIT DEPENDENCIES ON THE 'SP' PACKAGE
## THAT WILL BE DEPRECATED IN OCTOBER 2023. ONLY CHANGE INCLUDES ADDITIONAL COMMENTING OUT 
## OF ALL FOLLOWING CODE, WHICH IS USED TO BUILD EXTENT MASKS (PRESENT ON LINE 106 & 107):
    ##  template <- sp::spTransform(FedData::polygon_from_extent(template,proj4string='+proj=longlat +datum=WGS84 +ellps=WGS84'),
    ##       sp::CRS('+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs'))
## EXTENT MASKS CAN BE BUILT OUTSIDE OF THIS ANALYSIS WITH FOLLOWING CODE:
    ## ex_poly <- st_read("./DATA/EXAMPLE_POLYGON_STUDYAREA.SHP)
    ## ex_poly <- project(ex_poly,"+proj=longlat +datum=WGS84 +ellps=WGS84")
    ## ex_poly_extent <- st_as_sf(ex_poly)%>%
    ##  extent() %>%
    ##  as(.,"SpatialPolygons)
## THEN ex_poly_extent CAN BE USED AS THE @param template IN paleocar_update FUNCTION ONCE CALLED

## NOTE: RASTER ALSO DEPENDS ON SP, ALTHOUGH IT IS NOT BE DEPRECATED, AND SHOULD BE FULLY REPLACED 
## AT SOME POINT BY TERRA PACKAGE. AS TERRA BECOMES MORE PRODUCTIVE, RASTER-BASED FUNCTIONS WILL ALSO
## BE REPLACED

## ALSO REMOVED GDALUTILS LIBRARY AS NOT NECESSARY WITH SF LIBRARY


      ### END OF NOTES FOR KELSEY ###


library(raster)
library(FedData)
library(ncdf4)
library(magrittr)
library(foreach)
library(utils)

###################################################################################
globalVariables(c("tile", "type"))
#' Download, mask, and calculate the maize growing niche from the Bocinsky 2016 data
#' available from the NOAA paleoclimate database:
#' 
#' Bocinsky, R. K., J. Rush, K. W. Kintigh, and T. A. Kohler. 2016. 
#' Exploration and exploitation in the macrohistory of the pre-Hispanic Pueblo Southwest. 
#' Science Advances 2:e1501532.
#'
#' @param template A Spatial* or Raster* object from which to define the area of interest.
#' If omitted, download entire reconstruction.
#' @param label A character string naming the study area.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/NED/'.
#' @param extraction.dir A character string indicating where the extracted and cropped DEM should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/NED/'.
#' @param prcp_threshold The minimum amount of water-year precipitation, in mm, required to be in the farming niche.
#' Defaults to 300 mm. Use 'NA' to suppress niche calculation.
#' @param gdd_threshold The minimum number of Fahrenheit growing degree daysrequired to be in the farming niche.
#' Defaults to 1800 FGDD. Use 'NA' to suppress niche calculation.
#' @param years An integer vector of years, between AD 1 and 2000, you wish to extract.
#' Defaults to 1:2000.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A logical RasterBrick of the niche, through time.
#' @importFrom magrittr %>% %<>%
#' @importFrom foreach %do%
#' @importFrom utils head tail
#' @importFrom purrr map
#' @export
# paleocar <- function(template = NULL,
#                      label = 'PROJECT_NAME',
#                      raw.dir = paste0('./ANALYSIS/PALEOCAR/RAW/'),
#                      extraction.dir = paste0('./ANALYSIS/PALEOCAR/EXTRACTIONS/'),
#                      prcp_threshold = 300,
#                      gdd_threshold = 1800,
#                      years = 1:2000,
#                      force.redo = F)
paleocar <- function(template = NULL,
                     label = 'PROJECT_NAME',
                     raw.dir = paste0('./ANALYSIS/PALEOCAR/RAW/'),
                     extraction.dir = paste0('./ANALYSIS/PALEOCAR/EXTRACTIONS/'),
                     prcp_threshold = 300,
                     gdd_threshold = 1800,
                     years = 1:2000,
                     force.redo = F){
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo & file.exists(paste0(extraction.dir,"/","niche_",head(years,1),"-",tail(years,1),".tif")))
    return(raster::brick(paste0(extraction.dir,"/","niche_",head(years,1),"-",tail(years,1),".tif")))
  
  # Force Raster to load large rasters into memory
  raster::rasterOptions(chunksize=2e+07,maxmemory=2e+08)
  
  url <- 'http://www1.ncdc.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/usa/bocinsky2016/'
  req <- httr::GET(url)
  files <- XML::readHTMLTable(rawToChar(req$content),stringsAsFactors = FALSE)[[1]]$Name # Get the file listing
  files <- files[grep("nc4",files)] # Only download netcdf files
  
  # template <- sp::spTransform(FedData::polygon_from_extent(template),
  #                             sp::CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
  
  if(!is.null(template)){
    # template <- sp::spTransform(FedData::polygon_from_extent(template,proj4string='+proj=longlat +datum=WGS84 +ellps=WGS84'),
    #                             sp::CRS('+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs'))
    extent.latlon <- raster::extent(template)
    
    # Open USGS NED download service.
    # NED tiles are labeled by their northwest corner. Thus, coordinate 36.42N, -105.71W is in grid n37w106
    wests <- seq(ceiling(abs(extent.latlon@xmax)),ceiling(abs(extent.latlon@xmin)))
    norths <- seq(floor(abs(extent.latlon@ymin)),floor(abs(extent.latlon@ymax)))
    
    tilesLocations <- as.matrix(expand.grid(norths,wests,stringsAsFactors = FALSE)) %>%
      apply(1,function(x){paste0(x[2],"W",x[1],"N")})
    
    
    message("Area of interest includes ",length(tilesLocations)," reconstruction tiles.")
    
    files <- sapply(tilesLocations, grep, x = files, value = T) %>% as.vector()
  }
  
  # out <- foreach::foreach(tile = paste0(url,files)) %do% {
  #   FedData::download_data(url = tile,
  #                          destdir = raw.dir)
  # }
  
  download_without_overwrite <- function(url, folder)
  {
    filename <- basename(url)
    base <- tools::file_path_sans_ext(filename)
    ext <- tools::file_ext(filename)
    
    file_exists <- grepl(base, list.files(folder), fixed = TRUE)
    
    if (any(file_exists))
    {
      filename <- paste0(base, " (", sum(file_exists), ")", ".", ext)
    }
    
    download.file(url, file.path(folder, filename), mode = "wb", method = "libcurl")
  }

  out <- foreach::foreach(tile = paste0(url,files)) %do% {
    download_without_overwrite(url = tile,
                               folder = raw.dir )
  }
  
  # mosaick
  out_bricks <- foreach::foreach(type = c("PPT","GDD")) %do% {
    
    if(force.redo | !file.exists(paste0(extraction.dir,"/",type,"_",head(years,1),"-",tail(years,1),".tif"))){
      
      # system(paste0("gdalbuildvrt temp.vrt ", paste0(raw.dir,"/",grep(type, files, value = T), collapse=" ")))
      # system(paste0("gdal_translate -q -ot UInt16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND temp.vrt ",
      #               paste0("'",extraction.dir,"/",type,"_merged.tif","'")))
      # unlink("temp.vrt")
      
      # tile_brick <- suppressWarnings(raster::brick(paste0(extraction.dir,"/",type,"_merged.tif")))
      
      fname <- paste0(raw.dir,'/',grep(type,files,value=T),collapse=' ')
      nc<-nc_open(fname)
      tile_brick <- suppressWarnings(raster::brick(fname,varname=type))
      
      if(!is.null(template)){
        tile_brick %<>%
          raster::crop(template) %>%
          raster::mask(template)
      }
      
      raster::writeRaster(tile_brick, paste0(extraction.dir,"/",type,"_",head(years,1),"-",tail(years,1),".tif"),
                          datatype="INT2U",
                          options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                          overwrite=T,
                          setStatistics=FALSE)
      
      unlink(paste0(extraction.dir,"/",type,"_merged.tif"))
      
    }
    
    return(
      raster::brick(paste0(extraction.dir,"/",type,"_",head(years,1),"-",tail(years,1),".tif"))
    )
  }
  
  names(out_bricks) <- c("PPT","GDD")
  
  out_bricks %<>% purrr::map(function(x){
    projection(x) <- CRS('+proj=longlat +datum=WGS84')
    return(x)
  }) 
  
  if(!is.na(prcp_threshold) & !is.na(gdd_threshold)){
    
    if(force.redo | !file.exists(paste0(extraction.dir,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"))){
      ppt_niche <- out_bricks$PPT >= prcp_threshold
      
      raster::writeRaster(ppt_niche, paste0(extraction.dir,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"),
                          datatype="INT1U",
                          options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                          overwrite=T,
                          setStatistics=FALSE)
    }
    ppt_niche <- raster::brick(paste0(extraction.dir,"/PPT_niche_",head(years,1),"-",tail(years,1),".tif"))
    
    
    if(force.redo | !file.exists(paste0(extraction.dir,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"))){
      gdd_niche <- out_bricks$GDD >= gdd_threshold
      
      raster::writeRaster(gdd_niche, paste0(extraction.dir,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"),
                          datatype="INT1U",
                          options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                          overwrite=T,
                          setStatistics=FALSE)
      
      
    }
    gdd_niche <- raster::brick(paste0(extraction.dir,"/GDD_niche_",head(years,1),"-",tail(years,1),".tif"))
    
    # out <- precip_niche * gdd_niche
    
    out <- ppt_niche * gdd_niche
    
    raster::writeRaster(out, paste0(extraction.dir,"/","niche_",head(years,1),"-",tail(years,1),".tif"),
                        datatype="INT1U",
                        options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                        overwrite=T,
                        setStatistics=FALSE)
    
    out_bricks$niche <- out
  }
  
  return(out_bricks)
}
###################################################################################
# growing.niche <- get_bocinsky2016(template = extent(mv.checkdam.channels.longlat),
#                                   label = 'REESE_2018_JCAA',
#                                   raw.dir = './ANALYSIS/PALEOCAR/RAW',
#                                   extraction.dir = paste0('./ANALYSIS/PALEOCAR/EXTRACTIONS/'),
#                                   prcp_threshold = 300,
#                                   gdd_threshold = 1800,
#                                   years = 1:2000,
#                                   force.redo = T)
###################################################################################
