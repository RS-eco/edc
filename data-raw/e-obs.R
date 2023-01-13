#' ECAD/E-OBS Data for Europe

#' Ensemble version:
#' "We acknowledge the E-OBS dataset from the EU-FP6 project UERRA (http://www.uerra.eu) and 
#' the Copernicus Climate Change Service, and the data providers in the ECA&D project (https://www.ecad.eu)"

#' "Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: 
#' An Ensemble Version of the E-OBS Temperature and Precipitation Datasets, J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200"

#' E-OBS Data was downloaded from: https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=form
#' 
#' But see here for more information on the data: http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php

rm(list=ls()); gc()

##########

### Run in R, not in RStudio!!!

##########

# Load processNC package
#devtools::install_github("RS-eco/processNC")
library(processNC)

# Load other packages
library(terra)

# Set working directory
setwd("/home/matt/Documents/edc/")

# Set file directory
filedir <- "/home/matt/Documents/e-obs"

# Load europe outline
load("data/eur_countries.rda")

# List files
files <- list.files(filedir, pattern=".nc", full.names=T)

# Load, crop & mask elevation data
elev <- terra::rast(files[grep(files, pattern="elev")])
elev <- terra::crop(elev, vect(sf::st_as_sf(eur_countries)), snap="out", mask=T)
plot(elev); gc()
plot(sf::st_geometry(eur_countries), add=T)
writeCDF(elev ,filename=paste0("rawdata/", sub("reg_v23.1e.nc", "", 
                                               basename(files[grep(files, pattern="elev")])), 
                               "eur.nc"), overwrite=T)
rm(elev); gc()

# Specify file
lapply(2:length(files), function(z){
  if(!file.exists(paste0("rawdata/", sub("reg_v23.1e.nc", "", basename(files[z])), 
                         "eur_yearmon.nc"))){
    print(files[z])
    name <- c(NA, "pr", "pr", "tas", "tas", "tasmin", "tasmin", "tasmax", "tasmax")[z]
    #name <- c("tasmin", "tasmin", "tasmax", "tasmax")[z]
    print(name)
    
    # Load and summarise climate data
    #years <- seq(1950,2020, by=5)
    #dat <- lapply(years, function(y){
    #  temp <- tempfile(fileext=".nc")
    #  if(y == 2020){
    #    #summariseNC(files[z], ext=vect(sf::st_as_sf(europe)), cores=5,
    #    #            group_col="month", name=name, startdate=y, enddate=y)
    #  } else{
    #    #summariseNC(files[z], ext=vect(sf::st_as_sf(europe)), cores=5,
    #    #            group_col="month", name=name, startdate=y, enddate=y+4)
    #  }
    #  temp2 <- terra::rast(temp); rm(temp)
    #  terra::mask(temp2, vect(sf::st_as_sf(europe)))
    #}); gc()
    
    years <- seq(1950,2020, by=5)
    dat <- lapply(years, function(y){
      temp <- tempfile(fileext=".nc")
      if(y == 2020){
        aggregateNC(files[z], outfile=temp, group_col="yearmon", var=name, startdate=y, enddate=y)
      } else{
        aggregateNC(files[z], outfile=temp, group_col="yearmon", var=name, startdate=y, enddate=y+4)
      }
      temp2 <- terra::rast(temp); rm(temp); gc()
      temp2 <- terra::crop(temp2, vect(sf::st_as_sf(eur_countries)), snap="out", mask=T); gc()
      return(temp2)
    }); gc()
    dat <- terra::rast(dat); gc()
    writeCDF(x=dat,filename=paste0("rawdata/", sub("reg_v23.1e.nc", "", basename(files[z])), "eur_yearmon.nc"), 
             overwrite=T)
    rm(dat); gc()
  }
})

dat1 <- terra::rast("rawdata/elev_ens_0.1deg_eur.nc")
dat1_df <- as.data.frame(dat1, xy=T)
colnames(dat1_df)
assign(paste0("e-obs_elev_ens_eur_p1deg"), value=dat1_df)
save(list=paste0("e-obs_elev_ens_eur_p1deg"), 
     file=paste0("data/e-obs_elev_ens_eur_p1deg.rda"), compress="xz")

files <- list.files("rawdata/", pattern="mean_0.1deg_eur_yearmon.nc", full.names=T)
lapply(1:length(files), function(z){
  if(!file.exists(paste0("inst/extdata/e-obs_", sub(".nc", ".rda", basename(files[z]))))){
    dat <- terra::rast(files[z])
    if(z > 1){
     dat2 <- terra::rast(files[1])
     dat <- terra::mask(dat, dat2)
    }
    dat <- dat[[terra::time(dat) < "2021-01-15 UTC"]]
    dat_df <- as.data.frame(dat, xy=T)
    colnames(dat_df)[3:ncol(dat_df)] <- as.character(zoo::as.yearmon(terra::time(dat)))
    assign(paste0("e-obs_",  sub(".nc", "", basename(files[z]))), value=dat_df)
    save(list=paste0("e-obs_",  sub(".nc", "", basename(files[z]))), 
         file=paste0("inst/extdata/e-obs_",  sub(".nc", ".rda", basename(files[z]))), compress="xz")
  }
})

## Calculate bioclimatic variables

rm(list=ls()); gc()

library(dplyr); library(tidyr)

# Load data
rr_1967_1997 <- get(load("inst/extdata/e-obs_rr_ens_mean_0.1deg_eur_yearmon.rda")) %>%
  pivot_longer(names_to="yearmon", values_to="val", -c(x,y)) %>% 
  mutate(yr = substr(yearmon,5,8), mon = substr(yearmon,1,3)) %>%
  filter(yr %in% c(1967:1997)) %>% select(-c(yearmon, yr)) %>% 
  group_by(x, y, mon) %>% summarise(mn=mean(val)); gc()
rm("e-obs_rr_ens_mean_0.1deg_eur_yearmon"); gc()
tn_1967_1997 <- get(load("inst/extdata/e-obs_tn_ens_mean_0.1deg_eur_yearmon.rda")) %>% 
  pivot_longer(names_to="yearmon", values_to="val", -c(x,y)) %>% 
  mutate(yr = substr(yearmon,5,8), mon = substr(yearmon,1,3)) %>%
  filter(yr %in% c(1967:1997))  %>% select(-c(yearmon, yr)) %>% 
  group_by(x, y, mon) %>% summarise(mn=mean(val))
rm("e-obs_tn_ens_mean_0.1deg_eur_yearmon"); gc()
tx_1967_1997 <- get(load("inst/extdata/e-obs_tx_ens_mean_0.1deg_eur_yearmon.rda")) %>% 
  pivot_longer(names_to="yearmon", values_to="val", -c(x,y)) %>% 
  mutate(yr = substr(yearmon,5,8), mon = substr(yearmon,1,3)) %>%
  filter(yr %in% c(1967:1997)) %>% select(-c(yearmon, yr)) %>% 
  group_by(x, y, mon) %>% summarise(mn=mean(val)); gc()
rm("e-obs_tx_ens_mean_0.1deg_eur_yearmon"); gc()

prec <- rr_1967_1997 %>% group_by(x,y) %>% tidyr::spread(mon, mn)
prec_xy <- prec %>% ungroup() %>% dplyr::select(c(x,y))
prec <- prec %>% ungroup() %>% dplyr::select(-c(x,y))
tasmin <- tn_1967_1997 %>% group_by(x,y) %>% tidyr::spread(mon, mn) %>% 
  ungroup() %>% dplyr::select(-c(x,y)) 
tasmax <- tx_1967_1997 %>% group_by(x,y) %>% tidyr::spread(mon, mn) %>% 
  ungroup() %>% dplyr::select(-c(x,y))

`e-obs_bioclim_1967_1997_eur` <- dismo::biovars(as.matrix(prec), as.matrix(tasmin), as.matrix(tasmax)) %>%
  as.data.frame()
rm(rr_1967_1997, tn_1967_1997, tx_1967_1997, prec, tasmin, tasmax); invisible(gc())

`e-obs_bioclim_1967_1997_eur`$x <- prec_xy$x
`e-obs_bioclim_1967_1997_eur`$y <- prec_xy$y
head(`e-obs_bioclim_1967_1997_eur`)
summary(`e-obs_bioclim_1967_1997_eur`)

#' Save to file
save(list="e-obs_bioclim_1967_1997_eur", file="data/e-obs_bioclim_1967_1997_eur.rda", compress="xz")
rm(prec_xy, `e-obs_bioclim_1967_1997_eur`); gc()

rr_1987_2017 <- get(load("inst/extdata/e-obs_rr_ens_mean_0.1deg_eur_yearmon.rda")) %>%
  pivot_longer(names_to="yearmon", values_to="val", -c(x,y)) %>% 
  mutate(yr = substr(yearmon,5,8), mon = substr(yearmon,1,3)) %>%
  filter(yr %in% c(1987:2017)) %>% select(-c(yearmon, yr)) %>% 
  group_by(x, y, mon) %>% summarise(mn=mean(val)); gc()
rm("e-obs_rr_ens_mean_0.1deg_eur_yearmon"); gc()
tn_1987_2017 <- get(load("inst/extdata/e-obs_tn_ens_mean_0.1deg_eur_yearmon.rda")) %>%
  pivot_longer(names_to="yearmon", values_to="val", -c(x,y)) %>% 
  mutate(yr = substr(yearmon,5,8), mon = substr(yearmon,1,3)) %>%
  filter(yr %in% c(1987:2017))  %>% select(-c(yearmon, yr)) %>% 
  group_by(x, y, mon) %>% summarise(mn=mean(val))
rm("e-obs_tn_ens_mean_0.1deg_eur_yearmon"); gc()
tx_1987_2017 <- get(load("inst/extdata/e-obs_tx_ens_mean_0.1deg_eur_yearmon.rda")) %>% 
  pivot_longer(names_to="yearmon", values_to="val", -c(x,y)) %>% 
  mutate(yr = substr(yearmon,5,8), mon = substr(yearmon,1,3)) %>%
  filter(yr %in% c(1987:2017)) %>% select(-c(yearmon, yr)) %>% 
  group_by(x, y, mon) %>% summarise(mn=mean(val))
rm("e-obs_tx_ens_mean_0.1deg_eur_yearmon"); gc()

prec <- rr_1987_2017 %>% group_by(x,y) %>% tidyr::spread(mon, mn)
prec_xy <- prec %>% ungroup() %>% dplyr::select(c(x,y))
prec <- prec %>% ungroup() %>% dplyr::select(-c(x,y))
tasmin <- tn_1987_2017 %>% group_by(x,y) %>% tidyr::spread(mon, mn) %>% 
  ungroup() %>% dplyr::select(-c(x,y))
tasmax <- tx_1987_2017 %>% group_by(x,y) %>% tidyr::spread(mon, mn) %>% 
  ungroup() %>% dplyr::select(-c(x,y))
`e-obs_bioclim_1987_2017_eur` <- dismo::biovars(as.matrix(prec), as.matrix(tasmin), 
                                                as.matrix(tasmax)) %>% as.data.frame()
rm(rr_1987_2017, tn_1987_2017, tx_1987_2017, prec, tasmin, tasmax); invisible(gc())
`e-obs_bioclim_1987_2017_eur`$x <- prec_xy$x
`e-obs_bioclim_1987_2017_eur`$y <- prec_xy$y
head(`e-obs_bioclim_1987_2017_eur`)
summary(`e-obs_bioclim_1987_2017_eur`)

#' Save to file
save(list="e-obs_bioclim_1987_2017_eur", file="data/e-obs_bioclim_1987_2017_eur.rda", compress="xz")
rm(prec_xy, `e-obs_bioclim_1987_2017_eur`); gc()
