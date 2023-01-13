#' ---
#' title: "Calculate Euro-Cordex bioclim data for Europe"
#' author: "RS-eco"
#' ---

rm(list=ls()); gc()
#remotes::install_github("RS-eco/processNC")
library(raster); library(processNC); library(dplyr); library(lubridate); 
library(dismo); library(tidyr); library(rnaturalearthhires)

# Load high resolution country shapefile
data(countries10)
countries10 <- sf::st_as_sf(countries10)

# Please also make sure the raster and dismo package are installed

# To run this code, you will also require the processNC package (available from Github.com/RS-eco/processNC) 
# and a Linux machine with cdo installed!

# Specify file directory
filedir <- "/home/matt/Documents/"
filedir <- paste0(filedir, "EURO_CORDEX")

####################

#' ### Region (center of boundaries):

#' ~ 27N - 72N, ~ 22W - 45E

#' ### For rotated polar RCMs (in rotated coordinates):

#' RotPole (198.0; 39.25)
#' TLC (331.79; 21.67)
#' Nx=106 Ny=103

#' ### Periods

#' Hindcast (ERA Interim): 1989 -2008
#' Control: 1951 - 2005 (1951- 1980, 1981 - 2010)
#' Scenario: 2006 - 2100 (2011-2040, 2041-71, 2071-2100)

# Domain
domain <- "EUR-11" # 0.11 degree

# Experiment
rcps <- c("rcp26", "rcp45", "rcp85")
# historical is not available for prAdjust and tasAdjust

# Driving model
#gcms <- c("CNRM-CERFACS-CNRM-CM5", "ICHEC-EC-EARTH", "IPSL-IPSL-CM5A-MR", "MOHC-HadGEM2-ES", "MPI-M-MPI-ESM-LR") 
gcms <- c("ICHEC-EC-EARTH", "IPSL-IPSL-CM5A-MR", "MPI-M-MPI-ESM-LR")
#' Models for which no files have been downloaded: "CCCma-CanESM2", "ECMWF-ERAINT", "IPSL-IPSL-CM5A-LR", 
#' "MIROC-MIROC5", "NCC-NorESM1-M", "NOAA-GFDL-GFDL-ESM2G"

# Ensemble
#ensembles <- c("r1i1p1", "r2i1p1", "r3i1p1", "r12i1p1")
ensembles <- "r1i1p1"

# RCM Model
#rcm_models <- c("ARPEGE51", "CCLM4-8-17", "HIRHAM5", "RACMO22E", "RCA4", "REMO2009", "WRF331F")
rcm_models <- c("RACMO22E", "RCA4", "REMO2009")

# Downscaling realisation
rs <- "v1"
#rs <- c("v1", "v1a", "v2")

# Time Frequency
tm_freq <- "mon"

# Variable
#vars <- c("prAdjust", "tasAdjust", "tasmaxAdjust", "tasminAdjust")
vars <- c("prAdjust", "tasmaxAdjust", "tasminAdjust")

# Time period (of observed data???)
year_period <- c("1989-2010") # Is constant!!!

#' ## Identify available file combinations
rcms_long <- c("CLMcom-CCLM4-8-17", "KNMI-RACMO22E", "MPI-CSC-REMO2009", "SMHI-RCA4", "DMI-HIRHAM5")
rcms_long <- c("KNMI-RACMO22E", "MPI-CSC-REMO2009", "SMHI-RCA4")

#' ## Subset month period according to 30-yr periods
month_period <- c(paste0(c(1991,2001,2011,2041,2051,2061,2071,2081,2091), "01-", 
                         c(2000,2010,2020,2050,2060,2070,2080,2090,2100, by=10), "12"),
                  "209101-209911", "209101-209912")
combinations <- expand.grid(variable=vars, gcm=gcms, ensemble=ensembles, rcp=rcps, rcm=rcms_long, rs=rs, time_period=month_period)
unique(combinations$time_period)

#' ## Add time_period to combinations
combinations$time_frame <- ifelse(combinations$time_period %in% c("199101-200012", "200101-201012", "201101-202012"), "1991-2020", 
                                  ifelse(combinations$time_period %in% c("204101-205012", "205101-206012", "206101-207012"), "2041-2070", 
                                         "2071-2100"))
unique(combinations$time_frame)

filenames <- combinations %>% rowwise() %>%
  do(filename = paste0(.$variable, "_", domain, "_", .$gcm, "_", .$rcp, "_", .$ensemble, "_", 
                       .$rcm, "_", .$rs, "-SMHI-DBS45-MESAN-", year_period, "_", 
                       tm_freq, "_", .$time_period, ".nc"))
combinations$filename <- unlist(filenames$filename)

avail_combinations <- combinations[sapply(filenames, function(x){x %in% list.files(filedir, pattern=".nc")}) == TRUE,]; rm(filenames)
avail_combinations$filename <- paste0(filedir, "/", avail_combinations$filename)
avail_combinations <- avail_combinations %>% arrange(filename)
head(avail_combinations$filename)

#' ## Remap data and save to file
lapply(1:nrow(avail_combinations), function(x){
  if(!file.exists(sub("EURO_CORDEX/", "EURO_CORDEX/EUR/", avail_combinations$filename[x]))){
    if(!file.exists(sub(".nc", "_remap.nc", avail_combinations$filename[x]))){
      remapNC(gridfile=list.files(paste0(system.file(package="processNC"), "/extdata"), 
                                  pattern="euro-cordex_grid_eur.txt", full.names=T), infile=avail_combinations$filename[x],
              outfile=sub("[.]nc", "_remap.nc", avail_combinations$filename[x]))
    }
    dat <- raster::stack(sub("[.]nc", "_remap.nc", avail_combinations$filename[x]), 
                         varname = paste(avail_combinations$variable[x])) %>% raster::mask(countries10)
    writeRaster(dat, filename=sub("EURO_CORDEX/", "EURO_CORDEX/EUR/", avail_combinations$filename[x]), 
                varname=paste(avail_combinations$variable[x]), overwrite=T)
    if(x != 1){
      file.remove(sub("[.]nc", "_remap.nc", avail_combinations$filename[x]))   
    }
  }
})
raster::plot(raster::stack(sub("EURO_CORDEX/", "EURO_CORDEX/EUR/", avail_combinations$filename[1]), 
                           varname=paste(avail_combinations$variable[1]))[[1]])
#europe <- countries10 %>% filter(REGION_UN == "Europe")
plot(sf::st_geometry(countries10), add=T)
#sp::plot(bavaria, add=T)

#' ## Merge data files of subsequent years
avail_combinations$outfile <- sub("EURO_CORDEX/", "EURO_CORDEX/EUR/", avail_combinations$filename)
comb_files <- avail_combinations %>% ungroup() %>% dplyr::select(-c(time_period,filename)) %>% 
  group_by(variable, gcm, ensemble, rcp, rcm, rs, time_frame) %>% tidyr::nest()
comb_files$infiles <- sapply(1:nrow(comb_files), function(x){list(comb_files$data[[x]]$outfile)})
comb_files$outfile <- paste0(filedir, "/EUR/", comb_files$variable, "_", comb_files$gcm, "_", comb_files$rcp, "_",
                             comb_files$ensemble, "_", comb_files$rcm, "_", comb_files$rs, "-SMHI-DBS45-MESAN-", year_period, "_", 
                             tm_freq, "_", comb_files$time_frame, ".nc")

# Read in all files of subsequent time periods and write to one file
library(lubridate)
lapply(1:nrow(comb_files), function(x){
  print(x)
  if(!file.exists(comb_files$outfile[x])){
    dat <- lapply(comb_files$infiles[[x]], stack)
    dat <- stack(dat)
    dat <- setZ(dat, as.Date(gsub("[.]", "-", sub("X", "", names(dat)))))
    date <- getZ(dat)
    
    # Subset date
    if(comb_files$time_frame[x] == "1991-2020"){
      dat <- dat[[which(lubridate::year(date) %in% c(1991:2015))]]
      comb_files$time_frame[x] <- "1991-2015"
      date <- date[which(lubridate::year(date) %in% c(1991:2015))]
    }
    
    if(comb_files$variable[x]== "prAdjust"){
      ndays <- lubridate::days_in_month(getZ(dat))
      dat <- calc(dat, fun=function(x){x*ndays*86400}) # Convert precipitation from kg m-2 s-1 to mm month-1
      dat <- setZ(dat, date)
      dat <- zApply(dat, by=month, fun=mean)
      #plot(dat)
      writeRaster(dat, filename=comb_files$outfile[x], overwrite=T, varname="prAdjust")
    } else if(comb_files$variable[x]== "tasAdjust"){
      dat <- calc(dat, fun=function(x)x-273.15) # Convert temperature to degree C
      dat <- setZ(dat, date)
      dat <- zApply(dat, by=month, fun=mean)
      writeRaster(dat, filename=comb_files$outfile[x], overwrite=T, varname="tasAdjust")
    } else if(comb_files$variable[x]== "tasminAdjust"){
      dat <- calc(dat, fun=function(x)x-273.15) # Convert temperature to degree C
      dat <- setZ(dat, date)
      dat <- zApply(dat, by=month, fun=mean)
      writeRaster(dat, filename=comb_files$outfile[x], overwrite=T, varname="tasminAdjust")
    } else if(comb_files$variable[x]== "tasmaxAdjust"){
      dat <- calc(dat, fun=function(x)x-273.15) # Convert temperature to degree C
      dat <- setZ(dat, date)
      dat <- zApply(dat, by=month, fun=mean)
      writeRaster(dat, filename=comb_files$outfile[x], overwrite=T, varname="tasmaxAdjust")
    }
    rm(dat); gc()
  }
})
raster::plot(raster::stack(comb_files$outfile[11], varname="prAdjust")[[1]])
plot(sf::st_geometry(countries10), add=T)

########################################

## Create bioclim files

comb_files_wide <- comb_files %>% dplyr::select(-c(data,infiles)) %>% spread(variable, outfile)

prec <- lapply(comb_files_wide$prAdjust, stack)
tasmin <- lapply(comb_files_wide$tasminAdjust, stack)
tasmax <- lapply(comb_files_wide$tasmaxAdjust, stack)

bioclim <- lapply(1:length(prec), function(x){dismo::biovars(prec[[x]], tasmin[[x]], tasmax[[x]])})
rm(prec, tasmin, tasmax); invisible(gc())

# Mask by europe outline
load("data/europe.rda")
bioclim <- lapply(1:length(bioclim), function(x){raster::mask(raster::crop(bioclim[[x]], as(europe, "Spatial")), as(europe, "Spatial"))})

# Turn into data.frame
cordex_bioclim_eur <- lapply(1:length(bioclim), function(x){as.data.frame(raster::rasterToPoints(bioclim[[x]]))}); invisible(gc())
cordex_bioclim_eur <- bind_rows(cordex_bioclim_eur, .id="id")
colnames(cordex_bioclim_eur) <- c("id", "x", "y", paste0("bio", 1:19))
comb_files_wide$id <- as.character(1:length(bioclim))
colnames(comb_files_wide)
cordex_bioclim_p11deg_eur <- cordex_bioclim_eur %>% left_join(comb_files_wide, by="id") %>% 
  dplyr::select(x, y, gcm, ensemble, rcm, rs, rcp, time_frame, bio1, bio2, bio3, bio4, bio5, bio6, bio7,
                bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19) %>%
  dplyr::mutate_at(paste0("bio", 1:19), round, 2) # Round numbers for better storage performance
head(cordex_bioclim_p11deg_eur)
summary(cordex_bioclim_p11deg_eur)

#' Save to file
save(cordex_bioclim_eur, file="data/cordex_bioclim_eur.rda", compress="xz")
#readr::write_csv(cordex_bioclim_eur, "/home/matt/Desktop/cordex_bioclim_eur.csv.xz")
rm(cordex_bioclim_eur); gc()

#####

#' Create bioclim data at 0.5° resolution

prec <- lapply(comb_files_wide$prAdjust, stack)
tasmin <- lapply(comb_files_wide$tasminAdjust, stack)
tasmax <- lapply(comb_files_wide$tasmaxAdjust, stack)

prec <- lapply(prec, function(x) aggregate(x, fact=4, fun=mean, expand=T))
tasmin <- lapply(tasmin, function(x) aggregate(x, fact=4, fun=min, expand=T))
tasmax <- lapply(tasmax, function(x) aggregate(x, fact=4, fun=max, expand=T))

bioclim <- lapply(1:length(prec), function(x){dismo::biovars(prec[[x]], tasmin[[x]], tasmax[[x]])})
rm(prec, tasmin, tasmax); invisible(gc())

# Mask by europe outline
load("data/europe.rda")
bioclim <- lapply(1:length(bioclim), function(x){raster::mask(raster::crop(bioclim[[x]], as(europe, "Spatial")), as(europe, "Spatial"))})

cordex_bioclim_eur_p44deg <- lapply(1:length(bioclim), function(x){as.data.frame(raster::rasterToPoints(bioclim[[x]]))}); invisible(gc())
cordex_bioclim_eur_p44deg <- bind_rows(cordex_bioclim_eur_p5deg, .id="id")
colnames(cordex_bioclim_eur_p44deg) <- c("id", "x", "y", paste0("bio", 1:19))
comb_files_wide$id <- as.character(1:length(bioclim))
colnames(comb_files_wide)
cordex_bioclim_eur_p44deg <- cordex_bioclim_eur_p44deg %>% left_join(comb_files_wide, by="id") %>% 
  dplyr::select(x, y, gcm, ensemble, rcm, rs, rcp, time_frame, bio1, bio2, bio3, bio4, bio5, bio6, bio7,
                bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19) %>%
  dplyr::mutate_at(paste0("bio", 1:19), round, 2) # Round numbers for better storage performance
head(cordex_bioclim_eur_p44deg)
summary(cordex_bioclim_eur_p44deg)

#' Save to file
save(cordex_bioclim_eur_p44deg, file="data/cordex_bioclim_eur_p44deg.rda", compress="xz")
rm(cordex_bioclim_eur_p44deg); gc()

#dat <- readr::read_csv(comb_files$outfile[1]) 
#dat %>% ggplot2::ggplot() + ggplot2::geom_raster(ggplot2::aes(x=x, y=y, fill=Jan)) + 
#  ggplot2::coord_sf()