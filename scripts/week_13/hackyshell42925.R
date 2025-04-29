# hacky 29 apr 2025

# install packages
# sf - point, line, polygon data 
# terra - raster data 
# leaflet - interactive mapping (like google)

library(tidyverse)
library(sf)
library(terra)
library(leaflet)

# turn off spherical geometry 


# plotting sites 

# https://www.rothamsted.ac.uk/national-capability/the-insect-survey
# Rothamsted Insect Survey 


# write as sf object 


# assign crs to the xy coord


# using rnaturalearth 



# using leaflet 


# monitoring network -----------------------------------------------------------
# using geodata to create a monitoring network that based on geodata
# objective is to evenly consider all environmental data when picking field sites

# write xy coord as sf object 

# set crs 



# define a sampling area
# to draw a circle around our site we need to use a cylindrical projection and 
# then reproject to wgs 84 -mercator projection

sample_extent<-st_buffer(st_transform(center_location, st_crs("+proj=merc")),
                         dist=100000) %>% st_transform(st_crs("EPSG:4326"))


# view site and sampling area in leaflet 








# load in raster data 
# population 




plot(population)
plot(center_location, add=T, col="red", pch=19)


# climate data ----------------------------------------------------------------- 
# https://www.climatologylab.org/datasets.html
tmax<-rast("terraclimate_tmax_2023.nc")
plot(tmax)
plot(tmax)[1]

# remotes::install_github("elizagrames/climetric")
library(climetric)
library(raster) # raster is the same as terra, but climetric was written using raster functions

climetric::generate_wget(extent(sample_extent), c("tmin", "tmax", "ppt", 
                                                  "pdsi", "swe"), 2020, 2020)

source("01_functions.R")

# tmax<-clim_load("terraclimate_tmax_2023.nc", sample_extent) <- this will take 5-10 min

tmax<-rast("tmax_NY.tif")
ppt<-rast("ppt_NY.tif")


# elevation


# distance matrix --------------------------------------------------------------
# resample raster - lowest res 
res(tmax)
res(population)

# resample all rasters to the saem extent, crs, and resolution to align data 
# we can't compare something with 1500x1500 cells to 200x200 - a crs misalignment 
# by a couple km can screw up what data is connected to what geo location


# create a raster stack 

names(r_stack)<-list("max temp", "precipitation", "elevation", "population")

plot(r_stack)

# arrange data in r_stack as a matrix 
# replace nan values with NA
dat<-as.matrix(r_stack)
dat[is.nan(dat)]<-NA

# something I am working on - weighting variables 
# dat[, 'population']<-dat[, 'population'] * 2

cov_matrix <- cov(dat, use="pairwise.complete.obs") + diag(0.01, ncol(dat))  # Regularized covariance
center_mean <- colMeans(dat, na.rm = TRUE)

# Compute Mahalanobis distance



# reference layer to pull spatial location from
ref_layer<-r_stack[[1]]
values(ref_layer)<-NA

# Create an empty raster to store distances
 

plot(mahal_raster)

# select sites by repeatedly choosing the highest dissimilarity from the previously 
# chosen site 
x <- -75.920647
y <- 42.098843

# cell from xy using ref_layer 


# for the first round of site selection cell.dat = distance matrix values


# use which max to find the most dissimilar cell from cell.dat


# run through site_select function for n_iterations
cells <- c(cell1, max.diff)  # set baseline cells (xy1, xy2)

# loop through (x) iterations of field site selection



# Get coordinates xy from cell


# Plot the raster
plot(r_stack[[1]])

# Plot the valid points (those within the sampling extent)


# make leaflet map 


# check sites against data for evenness ----------------------------------------
site_val<-terra::extract(r_stack, xy_values) %>% data.frame()

plot(r_stack[[i]], r_stack[[j]])
points(site_val[,i], site_val[,j], col="red",pch=19)


# plot density of data


# merge rasters and crop to set sampling area?






