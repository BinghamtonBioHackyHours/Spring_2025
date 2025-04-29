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
sf_use_s2(F)

# plotting sites 

# https://www.rothamsted.ac.uk/national-capability/the-insect-survey
# Rothamsted Insect Survey 

dat<-read.csv("Rothamsted_locations.csv")

xy<-sf::st_as_sf(data.frame(lat=dat$x, long=dat$y), coords=c("long", "lat"))
sf::st_crs(xy)<-sf::st_crs("EPSG:4326")

# using rnaturalearth 
uk<-rnaturalearth::ne_countries(country = c("United Kingdom", "Scotland"), returnclass = "sf")[1] %>% 
  st_transform(crs="EPSG:4326")

plot(uk)
xy<-st_transform(xy, st_crs(uk))
plot(uk)
plot(st_geometry(xy), col="red", pch=19, add=T)


# using leaflet 
leaflet() %>% addTiles()

leaflet() %>% addTiles() %>% addMarkers(data = xy)

# monitoring network -----------------------------------------------------------
# using geodata to create a monitoring network that evenly considers all variables 
# or unevenly, depending on what your objectives are. 

center_location<-sf::st_as_sf(data.frame(lat=42.098843, long=-75.920647), coords=c("long", "lat"))
sf::st_crs(center_location)<- sf::st_crs("+proj=latlon")

# define a sampling area
sample_extent<-st_buffer(st_transform(center_location, st_crs("+proj=merc")),
                                      dist=100000) %>% st_transform(st_crs("EPSG:4326"))

leaflet() %>% addTiles() %>% addMarkers(data=center_location) %>% 
  addPolygons(data=sample_extent, fillOpacity = 0.3, col='red')

# load in raster data 
# population 
population<-rast("gpop.tif")
poly<-st_transform(sample_extent, crs(population))
population<-crop(population, poly) %>% terra::project("EPSG:4326")
                 
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
elev<-rast("elevation.tif")
plot(elev)
poly<-st_transform(sample_extent, crs(elev))
elev<-crop(elev, poly) %>% project("EPSG:4326")

# distance matrix --------------------------------------------------------------

# resample raster - lowest res 
res(tmax)
res(population)

r_resampled <- lapply(list(tmax, ppt, elev, population), function(x) {
                             smooth_rast(x, sample_extent, tmax)  # Use x8 as reference raster
                           })

r_stack<-rast(r_resampled)
plot(r_stack)

names(r_stack)<-list("max temp", "precipitation", "elevation", "population")

dat<-as.matrix(r_stack)
dat[is.nan(dat)]<-NA

dat[, 'population']<-dat[, 'population'] * 2

cov_matrix <- cov(dat, use="pairwise.complete.obs") + diag(0.01, ncol(dat))  # Regularized covariance
center_mean <- colMeans(dat, na.rm = TRUE)

# Compute Mahalanobis distance
distance_matrix <- mahalanobis(dat, center_mean, cov_matrix) 
distance_matrix <- as.numeric(distance_matrix)

ref_layer<-r_stack[[1]]
values(ref_layer)<-NA


# Create an empty raster to store distances
mahal_raster <- ref_layer
values(mahal_raster) <- distance_matrix  

plot(mahal_raster)


# select sites by repeatedly choosing the highest dissimilarity from the previously 
# chosen site 
x <- -75.920647
y <- 42.098843
xy1<-cbind(x, y)
cell1<-cellFromXY(ref_layer, xy1)

cell.dat<-distance_matrix
max.diff<-which.max(cell.dat)

# run through site_select function for sites 3-8 
n_iterations <- 10
cells <- c(cell1, max.diff)  # set baseline cells (xy1, xy2)

for (i in 1:n_iterations) {  # loops through (x) iterations of field site selection
  cells <- site_select(cells, dat, ref_layer)
}

# Get coordinates for valid cells
xy_values <- xyFromCell(ref_layer, cells)

# Plot the raster
plot(r_stack[[3]])

# Plot the valid points (those within the sampling extent)
points(xy_values, pch = 13, cex = 1.5, col = "red")

leaflet() %>% addTiles() %>% addMarkers(data=xy_values) %>% addPolygons(data=sample_extent, fillOpacity = 0.2) 


# check sites against data for evenness ----------------------------------------
i<-1
j<-2

site_val<-terra::extract(r_stack, xy_values) %>% data.frame()

# to view all layers, highlight and run lines 271-274
i<-i+1
j<-j+1
plot(r_stack[[i]], r_stack[[j]])
points(site_val[,i], site_val[,j], col="red",pch=19)

plot(density(r_stack[[1]]))

abline(v = site_val[, 1], col="tomato2", lty=2)

# figure showing distribution between raster layers ----------------------------
# install.packages("MASS")
library(MASS)

parcoord(site_val, var.label = F, col = factor(site_val$max.temp))
mtext("Rank order of chosen sites", side = 2, line = 3) # This adds the y-axis label
mtext("Geospatial layers", side=1, line=3)

siteval2 <- apply(site_val, 2, function(x){
  as.numeric(scale(x))
})

for(i in 1:12){
  if(i==1){
    plot(siteval2[i,], type="l", ylim=c(-3,3))
  }else{
    lines(siteval2[i,], col=i)
  }
}





