# Altitude data for Europe

# EU-DEM & EU-DEM Hillshade Data was downloaded from: 
#https://ec.europa.eu/eurostat/de/web/gisco/geodata/reference-data/elevation/eu-dem/eu-dem-laea
#https://ec.europa.eu/eurostat/de/web/gisco/geodata/reference-data/elevation/eu-dem/hillshade

# Set working directory
setwd("C:/Users/Admin/Documents/edc")

# Load packages
library(terra); library(dplyr)
library(sp); library(sf)

load("data/europe.rda")

####################

## Load hillshade data
if(!file.exists("data/hillshade_eur_2km.rda")){
  if(!file.exists("extdata/hillshade_eur.tif")){
    filedir <- "extdata/EU_DEM_Mosaic_hillshade/"
    file <- list.files(filedir, pattern=".tif$", full.names=T)
    hillshade_eur <- rast(file)
    
    # Crop & mask data
    europe_laea <- st_transform(europe, crs=crs(hillshade_eur))
    europe_laea <- vect(europe_laea)
    hillshade_eur <- crop(hillshade_eur,europe_laea); invisible(gc())
    #hillshade_eur <- mask(hillshade_eur, europe_laea); invisible(gc())
    writeRaster(hillshade_eur, filename="extdata/hillshade_eur.tif"); rm(hillshade_eur); invisible(gc())
  } else{
    hillshade_eur <- terra::rast("extdata/hillshade_eur.tif")
    # Transform europe
    europe_laea <- st_transform(europe, crs=crs(hillshade_eur))
  }
  
  # Remove NA values
  #NAvalue(hillshade_eur) <- 0
  
  # Plot data
  #plot(hillshade_eur)
  #plot(europe_laea, add=T)
  
  # Turn into 2 km raster
  hillshade_eur_mean <- aggregate(hillshade_eur, fact=80, fun=mean, na.rm=T)
  hillshade_eur_median <- aggregate(hillshade_eur, fact=80, fun=median, na.rm=T)
  hillshade_eur_min <- aggregate(hillshade_eur, fact=80, fun=min, na.rm=T)
  hillshade_eur_max <- aggregate(hillshade_eur, fact=80, fun=max, na.rm=T)
  hillshade_eur_sd <- aggregate(hillshade_eur, fact=80, fun=sd,na.rm=T); rm(hillshade_eur); gc()
  hillshade_eur_2km <- rast(list(hillshade_eur_mean, hillshade_eur_median, hillshade_eur_min, 
                                 hillshade_eur_max, hillshade_eur_sd)); rm(hillshade_eur_mean, hillshade_eur_median, hillshade_eur_min, 
                                                                           hillshade_eur_max, hillshade_eur_sd); gc()
  hillshade_eur_2km <- mask(hillshade_eur_2km, vect(europe_laea)); gc()
  #plot(hillshade_eur_2km[[1]])
  #plot(europe_laea, add=T)
  hillshade_eur_2km
  
  # Turn into data.frame
  hillshade_eur_2km <- as.data.frame(hillshade_eur_2km, xy=T)
  colnames(hillshade_eur_2km) <- c("x", "y", "hillshade_mean", "hillshade_median", "hillshade_min", 
                                   "hillshade_max", "hillshade_sd")
  
  #library(ggplot2)
  #hillshade_eur_2km %>% ggplot() + geom_tile(aes(x=x, y=y, fill=hillshade_mean))
  
  # Save to file
  save(hillshade_eur_2km, file="data/hillshade_eur_2km.rda", compress="xz")
}

####################

## Load DEM data
if(!file.exists("extdata/eu_dem_eur.tif")){
  filedir <- "extdata/EU_DEM_mosaic_1000K/"
  file <- list.files(filedir, pattern=".tif$", full.names=T)
  eu_dem_eur <- rast(file)
  # Crop & mask data
  europe_laea <- st_transform(europe, crs=crs(eu_dem_eur))
  europe_laea <- vect(europe_laea)
  eu_dem_eur <- crop(eu_dem_eur,europe_laea); invisible(gc())
  writeRaster(eu_dem_eur, filename="extdata/eu_dem_eur.tif"); rm(eu_dem_eur); invisible(gc())
} else{
  eu_dem_eur <- terra::rast("extdata/eu_dem_eur.tif")
  # Transform europe
  europe_laea <- st_transform(europe, crs=crs(eu_dem_eur))
}

# Create 2km eu_dem data
if(!file.exists("data/eu_dem_eur_2km.rda")){
  # Plot data
  #plot(eu_dem_eur)
  #plot(europe_laea, add=T)
  
  # Turn into 2 km raster
  eu_dem_eur_mean <- aggregate(eu_dem_eur, fact=80, fun=mean, expand=T, na.rm=T)
  eu_dem_eur_median <- aggregate(eu_dem_eur, fact=80, fun=median, expand=T, na.rm=T)
  eu_dem_eur_min <- aggregate(eu_dem_eur, fact=80, fun=min, expand=T, na.rm=T)
  eu_dem_eur_max <- aggregate(eu_dem_eur, fact=80, fun=max, expand=T, na.rm=T)
  eu_dem_eur_sd <- aggregate(eu_dem_eur, fact=80, fun=sd, expand=T, na.rm=T); gc()
  
  eu_dem_eur_2km <- rast(list(eu_dem_eur_mean, eu_dem_eur_median, eu_dem_eur_min, 
                              eu_dem_eur_max, eu_dem_eur_sd)); rm(eu_dem_eur_mean, eu_dem_eur_median, eu_dem_eur_min, 
                                                                  eu_dem_eur_max, eu_dem_eur_sd); gc()
  eu_dem_eur_2km <- mask(eu_dem_eur_2km, vect(europe_laea))
  #plot(eu_dem_eur_1km[[1]])
  #plot(europe_laea, add=T)
  eu_dem_eur_2km
  
  # Turn into data.frame
  eu_dem_eur_2km <- as.data.frame(eu_dem_eur_2km, xy=T)
  colnames(eu_dem_eur_2km) <- c("x", "y", "altitude_mean", "altitude_median", "altitude_min", 
                                "altitude_max", "altitude_sd")
  
  #library(ggplot2)
  #eu_dem_eur_2km %>% ggplot() + geom_tile(aes(x=x, y=y, fill=altitude_mean))
  
  # Save to file
  save(eu_dem_eur_2km, file="data/eu_dem_eur_2km.rda", compress="xz")
}

eu_dem_eur <- project(eu_dem_eur, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
eu_dem_eur

if(!file.exists("data/eu_dem_eur_p44deg.rda")){
  # Turn into 0.44 deg raster
  eu_dem_eur_mean <- aggregate(eu_dem_eur, fact=1340, fun=mean, expand=T, na.rm=T)
  eu_dem_eur_median <- aggregate(eu_dem_eur, fact=1340, fun=median, expand=T, na.rm=T)
  eu_dem_eur_min <- aggregate(eu_dem_eur, fact=1340, fun=min, expand=T, na.rm=T)
  eu_dem_eur_max <- aggregate(eu_dem_eur, fact=1340, fun=max, expand=T, na.rm=T)
  eu_dem_eur_sd <- aggregate(eu_dem_eur, fact=1340, fun=sd, expand=T, na.rm=T); gc()
  
  eu_dem_eur_p44deg <- rast(list(eu_dem_eur_mean, eu_dem_eur_median, eu_dem_eur_min, 
                                 eu_dem_eur_max, eu_dem_eur_sd)); rm(eu_dem_eur_mean, eu_dem_eur_median, eu_dem_eur_min, 
                                                                     eu_dem_eur_max, eu_dem_eur_sd); gc()
  
  eu_dem_eur_p44deg <- mask(eu_dem_eur_p44deg, vect(europe))
  eu_dem_eur_p44deg
  
  r <- rast(res=c(0.44, 0.44), xmin=-44.25, xmax=65.31, ymin=22.09, ymax=71.37)
  eu_dem_eur_p44deg <- resample(eu_dem_eur_p44deg, r)
  
  # Turn into data.frame
  eu_dem_eur_p44deg <- as.data.frame(eu_dem_eur_p44deg, xy=T)
  colnames(eu_dem_eur_p44deg) <- c("x", "y", "altitude_mean", "altitude_median", "altitude_min", 
                                   "altitude_max", "altitude_sd")
  
  library(ggplot2)
  eu_dem_eur_p44deg %>% ggplot() + geom_tile(aes(x=x, y=y, fill=altitude_mean))
  
  # Save to file
  save(eu_dem_eur_p44deg, file="data/eu_dem_eur_p44deg.rda", compress="xz")
  
}

if(!file.exists("data/eu_dem_eur_p11deg.rda")){
  # Turn into 0.11 deg raster
  eu_dem_eur_mean <- aggregate(eu_dem_eur, fact=335, fun=mean, expand=T, na.rm=T)
  eu_dem_eur_median <- aggregate(eu_dem_eur, fact=335, fun=median, expand=T, na.rm=T)
  eu_dem_eur_min <- aggregate(eu_dem_eur, fact=335, fun=min, expand=T, na.rm=T)
  eu_dem_eur_max <- aggregate(eu_dem_eur, fact=335, fun=max, expand=T, na.rm=T)
  eu_dem_eur_sd <- aggregate(eu_dem_eur, fact=335, fun=sd, expand=T, na.rm=T); gc()
  
  eu_dem_eur_p11deg <- rast(list(eu_dem_eur_mean, eu_dem_eur_median, eu_dem_eur_min, 
                                 eu_dem_eur_max, eu_dem_eur_sd)); rm(eu_dem_eur_mean, eu_dem_eur_median, eu_dem_eur_min, 
                                                                     eu_dem_eur_max, eu_dem_eur_sd); gc()
  
  eu_dem_eur_p11deg <- mask(eu_dem_eur_p11deg, vect(europe))
  eu_dem_eur_p11deg
  
  r <- terra::rast(res=c(0.11, 0.11), xmin=-10.48, xmax=40.12, ymin=34.96, ymax=71.15)
  eu_dem_eur_p11deg <- resample(eu_dem_eur_p11deg, r)
  
  # Turn into data.frame
  eu_dem_eur_p11deg <- as.data.frame(eu_dem_eur_p11deg, xy=T)
  colnames(eu_dem_eur_p11deg) <- c("x", "y", "altitude_mean", "altitude_median", "altitude_min", 
                                   "altitude_max", "altitude_sd")
  
  library(ggplot2)
  eu_dem_eur_p11deg %>% ggplot() + geom_tile(aes(x=x, y=y, fill=altitude_mean)) +
    scale_fill_gradientn(colours = terrain.colors(255))
  
  # Save to file
  save(eu_dem_eur_p11deg, file="data/eu_dem_eur_p11deg.rda", compress="xz")
}
