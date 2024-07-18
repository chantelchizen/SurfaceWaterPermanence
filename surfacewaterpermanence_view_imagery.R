#### Surface Water Permanence - View Imagery
#### View surface water permanence for polygons using Sentinel-2 imagery from GEE

#### Last updated: June 8, 2024
#### Author: Chantel Chizen

# Output:
# 1. plot with pond permanence for a specific area

rm(list = ls())

# Packages ----

library(dplyr)
library(raster)
library(rgee)
library(reticulate)
library(sf)


# Constants ----

## image years ====
# years you want imagery for
years <- c(2017, 2018, 2019, 2020, 2021, 2022)

## image months ====
# months you want imagery for within a given year
start <- c(4, 5, 6, 7, 8, 9, 10)

## water frequency thresholds ====
# create bins for each month in the "start"
waterfreq_th <- c(-0.1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5)

## plot colours ====
ndwiParams <- list(min = 0,
                   max = 7,
                   palette = c("#71F348", "#59C346", "#469C5C", "#347981", "#1F4EA4", "#0C26CB", "#0000CB"))


# Functions for Image Acquisition ----

# get polygon bounding from coordinates
bbox_wrap <- function(x) {
  (st_bbox(x))
}

# check how many months in a year don't have an image
is_integer <- function(x) {
  is.integer(x) || (is.numeric(x) && x == as.integer(x))
}

# calculate NDWI from Sentinel-2 imagery
getNDWI <- function(image) {
  image$normalizedDifference(c("B3", "B8"))$rename("NDWI")
}

# determines the QA value for subsequent cloud masking using s2_clean()
getQABits <- function(image, qa) {
  # Convert decimal (character) to decimal (little endian)
  qa <-
    sum(2 ^ (which(rev(
      unlist(strsplit(as.character(qa), "")) == 1
    )) - 1))
  # Return a single band image of the extracted QA bits, giving the qa value.
  image$bitwiseAnd(qa)$lt(1)
}

# cloud masking
s2_clean <- function(img) {
  # Select only band of interest, for instance, B2,B3,B4,B5,B8,B11,B12
  img_band_selected <- img$select("B[2-5|8]")
  
  # quality band
  qa <- img$select("QA60")
  
  # Select pixels to mask
  quality_mask <- getQABits(qa, "110000000000")
  
  # Mask pixels with value zero.
  img_band_selected$updateMask(quality_mask)
}

# water mask thresholds

## Kaplan & Avdan (2017) (https://doi.org/10.5194/isprs-annals-IV-4-W4-271-2017)
# -0.15 to 1: water
# -1 to -0.15: no water

mask_water <- function(img) {
  mask <-
    img$gte(-0.15)$multiply(ee$Image(1))$add(img$lt(-0.15)$multiply(ee$Image(0)))
  return(mask)
  
}

## EOS thresholds for NDWI (https://eos.com/make-an-analysis/ndwi/)
# 0.2 to 1: Water surface,
# 0.0 to 0.2: Flooding, humidity,
# -0.3 to 0.0: Moderate drought, non-aqueous surfaces,
# -1 to -0.3: Drought, non-aqueous surfaces

# mask_water <- function(img) {
#   mask <-
#     img$gte(-0.20)$multiply(ee$Image(1))$add(img$lt(-0.20)$multiply(ee$Image(0)))
#   return(mask)
#
# }


# Set-Up rgee ----

# for the first rgee use on your computer follow the rgee set-up.R instructions

# ee_check() # Check non-R dependencies

# If rgee() or any non-R dependencies are not working, run ee_clean_pyenv()
# Restart R by closing
# Reopen R and run rgee::ee_install() and follow instructions

#ee_clean_pyenv()
#rgee::ee_install()
#ee_check() # Check non-R dependencies

# Initialize the rgee 
ee_Initialize()


# Load polygon data ----
## important column names to include:
## "wetland" - wetland ID number
## "geometry" - polygon geometry; automatically populates from a .shp file, CRS 4326 for bounding box to work properly

## provide polygon shape file path here ====
site_polygons <- st_read("site_polygons.shp") %>%
  mutate(wetland = as.numeric(wetland)) %>% # convert wetland # to numeric
  arrange(wetland) # reorder list in wetland order

st_crs(site_polygons)


# Get polygon bounding box ----

box_sf <- site_polygons %>%
  dplyr::group_by(wetland) %>%
  tidyr::nest()

box_sf <- box_sf %>%
  dplyr::mutate(bbox = purrr::map(data, bbox_wrap))

box_sf_coords <- subset(box_sf, select = c("wetland", "bbox"))


# Get NDWI images and values for the polygons ----

## for loop prep ====

# get lengths to create for loops
no_areas <- length(box_sf_coords$bbox)
no_years <- length(years)
no_months <- length(start)

# create a [n X p] matrix to add empty dataframes to
ndwi_img = as.list(1:(no_areas * no_years))
dim(ndwi_img) <- c(no_areas, no_years)

# create a [n X p] matrix to add monthly images to
image_number = as.list(1:(no_months))
img_list = as.list(1:(no_months))
dim(img_list) <- c(no_months)
rownames(img_list) <-
  c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

# month image list structure
month_img <- as.list(1:no_months)
dim(month_img) <- c(no_months)
rownames(month_img) <-
  c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

# create pond permanence dataframe structure

columns_wc_annual <- c("wetland", "year", "0", "1", "2", "3", "4", "5", "6", "7")
wetland_wc_annual_summary <- data.frame(matrix(nrow = 0, ncol = 10))
colnames(wetland_wc_annual_summary) = columns_wc_annual

# create wetland summary dataframe structure for consecutive months with water
wetland_consecutive_summary <-  data.frame(site_polygons) %>%
  mutate(
    year_2017 = NA,
    year_2018 = NA,
    year_2019 = NA,
    year_2020 = NA,
    year_2021 = NA,
    year_2022 = NA
  ) %>%
  arrange(wetland)

# create wetland summary dataframe structure for total months with water
wetland_total_summary <-  data.frame(site_polygons) %>%
  mutate(
    year_2017 = NA,
    year_2018 = NA,
    year_2019 = NA,
    year_2020 = NA,
    year_2021 = NA,
    year_2022 = NA
  ) %>%
  arrange(wetland)

## retrieve imagery and annual summary for each wetland ====

for (i in 1:no_areas) {
  print(paste("Working on", i, "of", no_areas))
  
  # select the wetland polygon
  sites <- site_polygons[i,]
  wetland_id <- sites$wetland
  
  wetland_site_consecutive_summary <-
    data.frame(subset(site_polygons, wetland == wetland_id))
  
  wetland_site_total_summary <-
    data.frame(subset(site_polygons, wetland == wetland_id))
  
  # convert wetland polygon to gee vector
  sites_gee <- sites %>%
    st_set_crs(4326) %>%
    sf_as_ee()
  
  # get the bbox coordinates for the wetland polygon
  bbox = box_sf_coords[i, ]$bbox[[1]]
  
  # create a bbox that is 2x larger than the wetland polygon
  # calculate the center
  center_x <- (bbox["xmin"] + bbox["xmax"]) / 2
  center_y <- (bbox["ymin"] + bbox["ymax"]) / 2
  
  # calculate the new width and height of the bbox
  width <- bbox["xmax"] - bbox["xmin"]
  height <- bbox["ymax"] - bbox["ymin"]
  
  # create the new bbox
  bbox <- c(
    xmin = center_x - width,
    xmax = center_x + width,
    ymin = center_y - height,
    ymax = center_y + height
  )
  
  names(bbox) <- c("xmin", "xmax", "ymin", "ymax")
  
  bbox = bbox_wrap(bbox)
  
  ROI <- c(
    bbox$xmin,
    bbox$ymin,
    bbox$xmax,
    bbox$ymin,
    bbox$xmax,
    bbox$ymax,
    bbox$xmin,
    bbox$ymax,
    bbox$xmin,
    bbox$ymin
  )
  
  # convert wetland bounding box to gee vector
  ROI <- matrix(ROI, ncol = 2, byrow = TRUE) %>%
    list() %>%
    st_polygon() %>%
    st_sfc() %>%
    st_set_crs(4326) %>%
    sf_as_ee()
  
  # get imagery for each year
  for (j in 1:no_years) {
    print(paste("Working on", years[j], "-", j, "of", no_years, "years"))
    
    # month image list structure - creates a blank one for each year
    month_img <- as.list(1:no_months)
    dim(month_img) <- c(no_months)
    rownames(month_img) <-
      c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    
    # acquire a median image of each month for a given year
    for (k in 1:no_months) {
      # find images acquired during a given month
      collection <-
        ee$ImageCollection("COPERNICUS/S2")$# retrieve Sentinel-2 imagery
        filterBounds(ROI)$# only retrieve within the wetland polygon
        filter(ee$Filter$calendarRange(years[[j]], years[[j]], "year"))$# filter for the given year
        filter(ee$Filter$calendarRange(start[[k]], start[[k]], "month"))$# filter for the given month
        filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", 20))$map(s2_clean)$sort("DATE_ACQUIRED")
      
      image_number[[k]] <- collection$size()$getInfo()
      print(image_number[[k]])
      #ifelse(image_number[[k]][1] > 0, print(ee_get_date_ic(collection, time_end = FALSE)$time_start), print("NA"))
      
      # take the median image for the month
      collection_median <- collection$median()
      
      # calculate the NDWI
      if (as.numeric(image_number[[k]]) > 0) {
        month_img[[k]] <- getNDWI(collection_median)$clip(sites_gee)
      } else {
        month_img[[k]] <- NA
      }
    }
    
    # create a list of the monthly images for a given year
    month_img_count <-
      month_img[!is.na(month_img) & !is.logical(month_img)]
    
    # check how many months we have images for
    no_month_img <- length(month_img_count)
    no_month_img
    
    missing_date <-
      setdiff(start, match(names(month_img_count), month.abb))
    paste(years[[j]], paste("missing month:", missing_date))
    
    ## interpolate missing months ====
    
    if (length(missing_date) > 1) {
      print("multiple months missing")
      
    } else if (length(missing_date) == 1 && missing_date < 5) {
      
      # identify the position to insert the missing month
      missing_month_pos <- which(start == missing_date)
      
      # create a constant image with the value 1
      previous_raster <- ee$Image(1)$clip(sites_gee)
      
      # create an average raster for the missing month
      next_raster <- month_img[[missing_month_pos + 1]]
      
      # ensure the constant image has the same bands and properties as the original image
      average_raster <-
        previous_raster$rename(next_raster$bandNames())
      average_raster$reproject(crs = next_raster$projection(),
                               scale = next_raster$select(0))
      
      # insert the average raster into the stack at the missing month position
      month_img[[missing_month_pos]] <- average_raster
      
    } else if (length(missing_date) == 1 &&
               missing_date > 4 && missing_date < 10) {
      
      # identify the position to insert the missing month
      missing_month_pos <- which(start == missing_date)
      
      # create an average raster for the missing month
      previous_raster <- month_img[[missing_month_pos - 1]]
      next_raster <- month_img[[missing_month_pos + 1]]
      average_raster <- previous_raster$add(next_raster)$divide(2)
      
      # insert the average raster into the stack at the missing month position
      month_img[[missing_month_pos]] <- average_raster
      
    } else if (length(missing_date) == 1 && missing_date == 10) {
      
      # identify the position to insert the missing month
      missing_month_pos <- which(start == missing_date)
      
      # create an average raster for the missing month
      previous_raster <- month_img[[missing_month_pos - 1]]
      
      # insert the average raster into the stack at the missing month position
      month_img[[missing_month_pos]] <- previous_raster
      
    } else if (length(missing_date) > 1 && missing_date[1] == 4) {
      
      # identify the position to insert the missing month
      missing_month_pos <- which(start == missing_date[1])
      
      # create a constant image with the value 1
      previous_raster <- ee$Image(1)$clip(sites_gee)
      
      # create an average raster for the missing month
      next_raster <- month_img[[missing_month_pos + 1]]
      
      # ensure the constant image has the same bands and properties as the original image
      average_raster <-
        previous_raster$rename(next_raster$bandNames())
      average_raster$reproject(crs = next_raster$projection(),
                               scale = next_raster$select(0))
      
      # insert the average raster into the stack at the missing month position
      month_img[[missing_month_pos]] <- average_raster
      
      
    } else {
      print("all months acquired")
      
    }
    
    
    # put months together in a year collection
    
    month_img <- month_img[!is.na(month_img) & !is.logical(month_img)]
    no_month_img <- length(month_img)
    
    no_month_img <- length(month_img)
    no_month_img
    
    ## create image collection for year ====
    
    collection_year <-
      ee$ImageCollection$fromImages(c(month_img[[1]]))
    
    for (l in 2:no_month_img) {
      temp_collection <- ee$ImageCollection$fromImages(c(month_img[[l]]))
      collection_year <- collection_year$merge(temp_collection)
    }
    
    ## apply the water mask ====
    masked_collection <- collection_year$map(mask_water)
    
    ## visualize the water mask ====
    Map$centerObject(sites_gee)
    Map$addLayer(year_sum, ndwiParams)
    
    ## total months with water for a given year ====
    
    # sum all the months in the collection for the given year
    year_sum <- masked_collection$reduce(ee$Reducer$sum())
    
    ## visualize the total months with surface water permanence for a given year ====
    Map$centerObject(sites_gee)
    Map$addLayer(year_sum, ndwiParams)
  }
}