# This script contains all the functions used in the project

# Debugging and alignment checks for smooth_rast
smooth_rast <- function(raster, extent, ref_rast) {
  extent <- st_transform(extent, crs(raster))
  
  # Transform, crop, and resample raster
  crop_rast <- crop(raster, extent)
  crs(crop_rast) <- "EPSG:4326"
  
  rast_resampled <- resample(crop_rast, ref_rast)  # Match resolution
  return(rast_resampled)
}

# Function to choose field sites 
site_select <- function(cells, dat, ref_layer) {
  
  baseline <- dat[c(cells), ]
  comparison <- dat[-c(cells), ]
  
  # regularized covariance matrix 
  reg.cov <- cov(baseline) + diag(0.01, nrow(cov(baseline)))
  
  # calculate distance matrix 
  new_matrix <- stats::mahalanobis(comparison, apply(baseline, 2, mean, na.rm = F), reg.cov)
  
  tmp.r <- as.matrix(ref_layer)
  # tmp.r[!is.na(tmp.r)] <- NA
  tmp.r[-c(cells)] <- new_matrix
  
  # Find the index of the site with the maximum distance
  max.diff <- which.max(tmp.r)
  
  cells_new <- c(cells, max.diff)
  
  return(cells_new)
}


clim_load<-function(nc, polygon_sf){
  x<-climetric::brick_nc(nc)
  ext.region<-raster::extent(x)
  dim_region<-dim(x)[1:2]
  raster_mean <- climetric::extract_data(x, "mean", "annual", 
                                         writefile = FALSE, layer_type = "month")
  mean_raster <- rast(array(as.numeric(raster_mean$annual$mean), dim = dim_region), 
                      extent = ext.region, crs = "+proj=longlat +datum=WGS84 +no_defs")
  
  # Align and crop raster with sample_extent
  poly <- st_transform(polygon_sf, crs(mean_raster))
  crop_raster <- crop(mean_raster, poly) %>% 
    project(.,"EPSG:4326")  # Ensure CRS consistency
  
  plot(crop_raster)
  return(crop_raster)
}




