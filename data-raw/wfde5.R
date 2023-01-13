#' WFDE5 Data for Europe

#' Data was downloaded from: https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.20d54e34?tab=form

rm(list=ls()); gc()

# Load processNC package
library(processNC)

# Load other packages
library(terra)

# Set file directory
filedir <- "/home/matt/Documents/WFDE5"

# Load europe outline
load("data/europe.rda")

# Specify variables
var <- c("Rainf", "Snowf", "Tair")[3]

# List files
files <- list.files(filedir, pattern=paste0(var, "_"), full.names=T)
head(files)

# Load climate data
#dat <- terra::rast(files[1])
#dat
#time(dat)[1:10]

# Check format

#rm(dat); gc()

if(!file.exists(paste0("inst/extdata/", var, "_all_yearmon.nc"))){
  # Loop through individual files
  all_dat <- lapply(files, function(x){
    # Load data file
    dat <- terra::rast(x)
    
    # Crop and mask by outline of Europe
    dat <- terra::mask(terra::crop(dat, europe), vect(sf::st_as_sf(europe)))
    
    if(var %in% c("Tair")){
      day_mean <- terra::tapp(dat, index="days", fun="mean"); rm(dat)
      terra::time(day_mean) <- as.Date(sub("X", "", names(day_mean)), format="%Y.%m.%d")
      avg <- terra::tapp(day_mean, index="months", fun="mean"); rm(day_mean); gc()
      terra::time(avg) <- zoo::as.yearmon(strsplit(basename(x), split="_")[[1]][4], format="%Y%m")
    } else if(var %in% c("Rainf", "Snowf")){
      day_mean <- terra::tapp(dat, index="days", fun="sum"); rm(dat)
      terra::time(day_mean) <- as.Date(sub("X", "", names(day_mean)), format="%Y.%m.%d")
      avg <- terra::tapp(day_mean, index="months", fun="sum"); rm(day_mean); gc()
      terra::time(avg) <- zoo::as.yearmon(strsplit(basename(x), split="_")[[1]][4], format="%Y%m")
    }
    return(avg)
  })
  all_dat <- terra::rast(all_dat)
  all_dat
  names(all_dat) <- basename(files)
  plot(all_dat[[1]])
  
  writeCDF(all_dat ,filename=paste0("inst/extdata/", var, "_all_yearmon.nc"))
}

# Specify variables
var <- c("Rainf", "Snowf", "Tair")

# Specify dates
dates <- as.Date(paste0(as.numeric(sapply(files, function(x) strsplit(basename(x), "_")[[1]][4])), "15"), "%Y%m%d")

# Load files
lapply(var, function(z){
  dat <- terra::rast(paste0("inst/extdata/", z, "_all_yearmon.nc"))
  time(dat) <- dates
  names(dat) <- as.character(zoo::as.yearmon(dates))
  dat_df <- as.data.frame(dat, xy=T)
  assign(paste0("wfde5_", tolower(z), "_eur"), value=dat_df)
  save(list=paste0("wfde5_", tolower(z), "_eur"), file=paste0("data/wfde5_", tolower(z), "_eur.rda"), compress="xz")
})
