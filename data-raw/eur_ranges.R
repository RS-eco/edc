#' ---
#' title: "European species ranges"
#' author: "RS-eco"
#' ---

rm(list=ls()); gc()

filedir <- "/home/matt/Documents/"

setwd("/home/matt/Documents/edc")

library(sf)
library(fasterize)
library(dplyr)

load("data/europe.rda")

#' ### Rasterize IUCN species shapefiles

#' IUCN data was downloaded from: https://www.iucnredlist.org/resources/spatial-data-download, by running the following code:

#if(!dir.exists("ODONATA")) {
#  download.file("http://spatial-data.s3.amazonaws.com/groups/FW_ODONATA.zip", 
#                destfile = paste0(filedir, "/IUCN/FW_ODONATA.zip")) # <-- 594.3 MB
#  unzip(paste0(filedir, "/IUCN/FW_ODONATA.zip"), exdir = paste0(filedir, "/IUCN"))
#  unlink(paste0(filedir, "/IUCN/FW_ODONATA.zip"))
#}

#' You have one shapefile for a group of animals, consisting of individual polygons for each species with different information (presence, origin, seasonal). You can specify the resolution in degrees, here we use 0.11°.

#+ rasterize_odonata, eval=F
dsn <- paste0(filedir, "/IUCN/FW_ODONATA.shp")
sp_ind_shp <- sf::st_read(dsn=dsn)
#sp_ind_shp <- sf::st_make_valid(sp_ind_shp)

dat_ec <- raster::stack(paste0(filedir, "EURO_CORDEX/EUR/prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r1i1p1_KNMI-RACMO22E_v1-SMHI-DBS45-MESAN-1989-2010_mon_209101-210012.nc"))
dat_ec <- raster::mask(raster::crop(dat_ec, as(europe, "Spatial")), as(europe, "Spatial"))
dat_ec <- dat_ec[[1]]

extent <- raster::extent(c(-10.48, 40.12, 34.85, 71.15))
resolution <- 0.11
crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
r <- raster::raster(ext=extent, resolution=resolution, crs=crs)
r_high <- raster::disaggregate(r, fact=10)
r_high

iucn_odonata_eur_p11deg <- lapply(1:nrow(sp_ind_shp), function(x){
  r_poly <- fasterize::fasterize(sp_ind_shp[x,], r_high, background=NA)
  r_poly <- raster::mask(raster::crop(r_poly, as(europe, "Spatial")), as(europe, "Spatial"))
  r_poly <- raster::aggregate(r_poly, fact=10, fun=sum)
  r_poly <- as.data.frame(raster::rasterToPoints(r_poly))
  colnames(r_poly) <- c("x", "y", "perc_present")
  if(nrow(r_poly)>0){
    r_poly$binomial <- sp_ind_shp$binomial[x]
    r_poly$presence <- sp_ind_shp$presence[x]
    r_poly$origin <- sp_ind_shp$origin[x]
    r_poly$seasonal <- sp_ind_shp$seasonal[x]
    return(r_poly)
  } else{
    return(NULL)
  }
  gc()
}); rm(sp_ind_shp); gc()
iucn_odonata_eur_p11deg <- Filter(Negate(is.null), iucn_odonata_eur_p11deg)
iucn_odonata_eur_p11deg <- dplyr::bind_rows(iucn_odonata_eur_p11deg)
# %>% group_by(binomial, presence, origin, seasonal) %>% 
save(iucn_odonata_eur_p11deg, file="data/iucn_odonata_eur_p11deg.rda", compress="xz")

#' ### Rasterize BirdLife species shapefiles

#+ rasterize_birdlife, eval=FALSE
dsn <- paste0(filedir, "/BirdLife_2018")
files <- list.files(dsn, pattern=".shp", full.names=T)

load("data/birdlife_eur_p11deg.rda")
spec <- unique(birdlife_eur_p11deg$binomial)
spec <- gsub(pattern=" ", replacement="_", spec)
files <- sapply(spec, function(x) files[grep(pattern=paste0(x, "[.]"), files)])

dat_ec <- raster::stack(paste0(filedir, "EURO_CORDEX/EUR/prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r1i1p1_KNMI-RACMO22E_v1-SMHI-DBS45-MESAN-1989-2010_mon_209101-210012.nc"))
dat_ec <- raster::mask(raster::crop(dat_ec, as(europe, "Spatial")), as(europe, "Spatial"))
dat_ec <- dat_ec[[1]]

extent <- raster::extent(c(-10.48, 40.12, 34.85, 71.15))
resolution <- 0.11
crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
r <- raster::raster(ext=extent, resolution=resolution, crs=crs)
r_high <- raster::disaggregate(r, fact=10)
r_high

#options(show.error.messages = FALSE)
options(try.outFile = stdout())
#x <- files[298]
dir.create(paste0(filedir, "Aves"))
lapply(files[386:length(files)], function(x){
  if(!file.exists(paste0(filedir, "Aves/", sub(".shp", ".csv", basename(x))))){
    print(x)
    dat <- sf::st_read(x)
    #try(dat <- sf::st_crop(dat, extent)) # Crop data and skip in case file has holes
    if(nrow(dat) > 0){
      pol_df <- lapply(1:nrow(dat), function(y){
        r_poly <- fasterize::fasterize(dat[y,], r_high, background=NA)
        r_poly <- raster::mask(raster::crop(r_poly, as(europe, "Spatial")), as(europe, "Spatial"))
        try(r_poly <- raster::aggregate(r_poly, fact=10, fun=sum))
        #try(r_poly <- raster::resample(r_poly, dat_ec))
        try(r_poly <- as.data.frame(raster::rasterToPoints(r_poly)))
        try(colnames(r_poly) <- c("x", "y", "perc_present"))
        if(nrow(r_poly)>0){
          r_poly$binomial <- dat$SCINAME[y]
          r_poly$presence <- dat$PRESENC[y]
          r_poly$origin <- dat$ORIGIN[y]
          r_poly$seasonal <- dat$SEASONA[y]
          return(r_poly)
        } else{
          return(NULL)
        }
      }); rm(dat); gc()
      pol_df <- Filter(Negate(is.null), pol_df)
      pol_df <- data.table::rbindlist(pol_df)
      if(nrow(pol_df) > 0){
        readr::write_csv(pol_df, paste0(filedir, "Aves/", 
                                        sub(".shp", ".csv", basename(x))))
      }
    }
  }
}); gc()
files <- list.files(paste0(filedir, "Aves/"), full.names=T)
birdlife_eur_p11deg <- vroom::vroom(files)

save(birdlife_eur_p11deg, file="data/birdlife_eur_p11deg.rda", compress="xz")
#file.remove(files)
#dir.remove(paste0(filedir, "Aves"))
