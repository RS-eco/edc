# Altitude data for Europe

# SRTM 90m Data was downloaded from: http://srtm.csi.cgiar.org/srtmdata/

# Set file directory
filedir <- "extdata/SRTM/"
#setwd(filedir)

# Load packages
library(raster); library(dplyr)
library(sp); library(sf)

load("data/europe.rda")
europe <- as(europe, "Spatial")

####################

## Process low resolution (90 m) DEM data

file <- list.files(filedir, pattern="cut", full.names=T)
alt_eur <- lapply(file, stack)
NAvalue(alt_eur[[1]]) <- -32768
NAvalue(alt_eur[[2]]) <- -32768

if(!file.exists("data/srtm_csi_eur_900m.rda")){
  alt_eur_mean <- lapply(alt_eur, function(x) raster::aggregate(x, fact=c(10,10), 
                                                                fun=mean, expand=T, na.rm=T))
  alt_eur_mean <- lapply(alt_eur_mean, function(x) raster::crop(x, europe))
  alt_eur_mean <- lapply(alt_eur_mean, function(x) raster::mask(x, europe)); gc()
  
  # Define CRS (WGS84)
  raster::crs(alt_eur_mean[[1]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(alt_eur_mean[[2]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  srtm_csi_eur_900m <- raster::merge(alt_eur_mean[[1]], alt_eur_mean[[2]], 
                                     tolerance=0.2); rm(alt_eur_mean); gc()
  
  # Plot data
  raster::plot(srtm_csi_eur_900m)
  plot(europe, add=T)
  
  # Turn into data.frame
  srtm_csi_eur_900m <- as.data.frame(rasterToPoints(srtm_csi_eur_900m))
  colnames(srtm_csi_eur_900m) <- c("x", "y", "altitude_m")
  
  library(ggplot2)
  srtm_csi_eur_900m %>% ggplot() + geom_tile(aes(x=x, y=y, fill=altitude_m))
  
  # Save to file
  save(srtm_csi_eur_900m, file="data/srtm_csi_eur_900m.rda", compress="xz")
}

if(!file.exists("data/srtm_csi_eur_p11deg.rda")){
  # Turn into 0.44 deg raster
  alt_eur_mean <- lapply(alt_eur, function(x) raster::aggregate(x, fact=c(132,132), 
                                                                fun=mean, expand=T, na.rm=T))
  alt_eur_mean <- lapply(alt_eur_mean, function(x) raster::crop(x, europe))
  alt_eur_mean <- lapply(alt_eur_mean, function(x) raster::mask(x, europe)); gc()
  
  # Define CRS (WGS84)
  raster::crs(alt_eur_mean[[1]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(alt_eur_mean[[2]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  srtm_csi_eur_p11deg <- raster::merge(alt_eur_mean[[1]], alt_eur_mean[[2]], 
                                       tolerance=0.2); rm(alt_eur_mean); gc()
  srtm_csi_eur_p11deg <- raster::crop(srtm_csi_eur_p11deg, europe); gc()
  srtm_csi_eur_p11deg <- raster::mask(srtm_csi_eur_p11deg, europe); gc()
  raster::plot(srtm_csi_eur_p11deg[[1]])
  plot(europe, add=T)
  srtm_csi_eur_p11deg
  
  r <- raster::raster(res=c(0.11, 0.11), xmn=-10.48, xmx=40.12, ymn=34.96, ymx=71.15)
  srtm_csi_eur_p11deg <- raster::resample(srtm_csi_eur_p11deg, r)
  
  # Turn into data.frame
  srtm_csi_eur_p11deg <- as.data.frame(rasterToPoints(srtm_csi_eur_p11deg))
  colnames(srtm_csi_eur_p11deg) <- c("x", "y", "alt_mean")
  
  library(ggplot2)
  srtm_csi_eur_p11deg %>% ggplot() + 
    geom_tile(aes(x=x, y=y, fill=alt_mean))
  
  # Save to file
  save(srtm_csi_eur_p11deg, file="data/srtm_csi_eur_p11deg.rda", compress="xz")
  removeTmpFiles(h=2); gc()
}

if(!file.exists("data/srtm_csi_eur_p44deg.rda")){
  # Turn into 0.44 deg raster
  alt_eur_mean <- lapply(alt_eur, function(x) raster::aggregate(x, fact=c(530,530), 
                                                                fun=mean, expand=T, na.rm=T)); rm(alt_eur)
  alt_eur_mean <- lapply(alt_eur_mean, function(x) raster::crop(x, europe))
  alt_eur_mean <- lapply(alt_eur_mean, function(x) raster::mask(x, europe)); gc()
  
  # Define CRS (WGS84)
  raster::crs(alt_eur_mean[[1]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(alt_eur_mean[[2]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  srtm_csi_eur_p44deg <- raster::merge(alt_eur_mean[[1]], alt_eur_mean[[2]], 
                                       tolerance=0.2); rm(alt_eur_mean); gc()
  srtm_csi_eur_p44deg <- raster::mask(srtm_csi_eur_p44deg, europe); rm(srtm_csi_eur_900m); gc()
  raster::plot(srtm_csi_eur_p44deg[[1]])
  plot(europe, add=T)
  srtm_csi_eur_p44deg
  
  r <- raster::raster(res=c(0.44, 0.44), xmn=-10.37, xmx=40.23, ymn=34.85, ymx=70.93)
  srtm_csi_eur_p44deg <- resample(srtm_csi_eur_p44deg, r)
  
  # Turn into data.frame
  srtm_csi_eur_p44deg <- as.data.frame(rasterToPoints(srtm_csi_eur_p44deg))
  colnames(srtm_csi_eur_p44deg) <- c("x", "y", "alt_mean")
  
  library(ggplot2)
  srtm_csi_eur_p44deg %>% ggplot() + geom_tile(aes(x=x, y=y, fill=alt_mean))
  
  # Save to file
  save(srtm_csi_eur_p44deg, file="data/srtm_csi_eur_p44deg.rda", compress="xz")
  removeTmpFiles(h=2); gc()
}
