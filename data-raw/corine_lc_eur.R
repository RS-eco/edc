#' ---
#' title: "Create Corine land cover data for Europe"
#' author: "RS-eco"
#' ---

# Downloaded from https://land.copernicus.eu/pan-european/corine-land-cover

rm(list=ls()); invisible(gc())

# Specify file directory
filedir <- "/home/matt/Documents/Corine/DATA"
#filedir <- "C:/Users/Admin/Documents/Corine/DATA"

# Set working directory
#setwd("C:/Users/Admin/Documents/edc")
setwd("/home/matt/Documents/edc")

# Load packages
library(terra); library(rworldmap)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(sf)

# Load world map
#devtools::install_github("ropensci/rnaturalearthhires")
data(countries10, package="rnaturalearthhires")
world <- sf::st_as_sf(countries10) %>% filter(CONTINENT == "Europe")
'%!in%' <- function(x,y)!('%in%'(x,y))
world <- world %>% filter(NAME %!in% c("Russia", "Iceland"))
world$NAME

box <- rgeos::bbox2SP(n = 72, s = 30, w = -12, e = 42,
                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
world <- st_crop(world, box)

# Load files
clc_files <- list.files(filedir, pattern=".*_CLC.*\\.tif$", full.names=T)

# Re-project outline/world data
world_laea <- sf::st_transform(world, "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

if(!file.exists("rawdata/corine_lc_eur.tif")){
  clc_dat <- terra::rast(clc_files)
  world_laea <- vect(world_laea)
  clc_sub <- crop(clc_dat, world_laea); rm(clc_dat); invisible(gc())
  clc_sub <- mask(clc_sub, world_laea); invisible(gc())
  clc_sub <- terra::project(clc_sub, "+proj=longlat +datum=WGS84")
  writeRaster(clc_sub, filename="rawdata/corine_lc_eur.tif"); rm(clc_sub); invisible(gc())
  gc(); terra::tmpFiles(remove=T)
}

# Layerize data and aggregate
years <- c(1990, 2000, 2006, 2012, 2018)

lapply(1:length(years), function(x){
  if(!file.exists(paste0("rawdata/corine_lc_eur_", years[x], "_p11deg.tif")) | 
     !file.exists(paste0("rawdata/corine_lc_eur_", years[x], "_p44deg.tif"))){
    dat <- raster::raster("rawdata/corine_lc_eur.tif", band=x)
    dat <- raster::layerize(dat); gc()
    
    # Layerize data
    res <- "p44deg"
    fact <- c(335,335)
    r <- raster::raster(res=c(0.44, 0.44), xmn=-10.37, xmx=40.23, ymn=34.85, ymx=70.93)
    dat2 <- raster::aggregate(dat, fact=fact, fun=sum, expand=T, na.rm=T); raster::removeTmpFiles(h=4); gc()
    dat2 <- dat2/(fact[1]*fact[2]/100)
    dat2 <- raster::resample(dat2, r, method="ngb"); gc(); raster::removeTmpFiles(h=4)
    raster::writeRaster(dat2, filename=paste0("rawdata/corine_lc_eur_", years[x], "_", res, ".tif"), 
                        overwrite=TRUE); rm(dat2); gc(); raster::removeTmpFiles(h=4)
    
    res <- "p11deg"
    fact <- c(84,84)
    r <- raster::raster(res=c(0.11, 0.11), xmn=-10.48, xmx=40.12, ymn=34.96, ymx=71.15)
    dat2 <- raster::aggregate(dat, fact=fact, fun=sum, expand=T, na.rm=T); raster::removeTmpFiles(h=4); gc()
    dat2 <- dat2/(fact[1]*fact[2]/100)
    dat2 <- raster::resample(dat2, r, method="ngb"); gc(); raster::removeTmpFiles(h=4)
    raster::writeRaster(dat2, filename=paste0("rawdata/corine_lc_eur_", years[x], "_", res, ".tif"), 
                        overwrite=TRUE); rm(dat, dat2); gc(); raster::removeTmpFiles(h=0.01)
  }
})

# Turn into data.frame
res <- "p11deg"
years <- c(1990, 2000, 2006, 2012, 2018)
dat2 <- lapply(1:length(years), function(x){
  dat <-  raster::stack(paste0("rawdata/corine_lc_eur_", years[x], "_", res, ".tif"))
  dat <- crop(dat, world)
  dat <- mask(dat, world); invisible(gc())
  plot(dat)
  
  clc_sub <- dat %>% raster::rasterToPoints() %>% as.data.frame()
  colnames(clc_sub) <- sub(pattern=paste0("corine_lc_eur_", years[x], "_", 
                                          res, "."), "", colnames(clc_sub))
  colnames(clc_sub)
  clc_sub$year <- years[x]
  
  # Re-structure data.frame & define categories
  clc_sub <- clc_sub %>% as.data.frame() %>% group_by(x,y,year) %>% 
    tidyr::pivot_longer(names_to="var", values_to="perc_cover", -group_cols()) %>% 
    tidyr::drop_na() %>% mutate(var = as.numeric(var))
  head(clc_sub)
  clc_sub <- clc_sub %>% mutate(var = factor(var, levels = c(1:45), 
                                             labels=c("Continuous urban fabric","Discontinuous urban fabric","Industrial or commercial units",
                                                      "Road and rail networks and associated land","Port areas","Airports","Mineral extraction sites","Dump sites",
                                                      "Construction sites","Green urban areas","Sport and leisure facilities","Non-irrigated arable land",
                                                      "Permanently irrigated land", "Rice fields", "Vineyards", "Fruit trees and berry plantations", "Olive groves",
                                                      "Pastures", "Annual crops associated with permanent crops", "Complex cultivation patterns",
                                                      "Land principally occupied by agriculture with significant areas of natural vegetation",
                                                      "Agro-forestry areas","Broad-leaved forest", "Coniferous forest", "Mixed forest", "Natural grasslands",
                                                      "Moors and heathland", "Sclerophyllous vegetation", "Transitional woodland-shrub", "Beaches dunes sands",
                                                      "Bare rocks", "Sparsely vegetated areas", "Burnt areas", "Glaciers and perpetual snow", "Inland marshes",
                                                      "Peat bogs", "Salt marshes", "Salines", "Intertidal flats", "Water courses", "Water bodies",
                                                      "Coastal lagoons", "Estuaries", "Sea and ocean", "NODATA")))
  return(clc_sub)
}); gc()
dat2 <- bind_rows(dat2) %>% ungroup()
dat3 <- dat2 %>% tidyr::pivot_wider(names_from = year, values_from = perc_cover, values_fill=0)
dat3 <- dat3[rowSums(dat3[,-c(1,2,3)], na.rm=T) != 0,]
head(dat3)

# Check range of values
summary(dat3)
if(res)
  dat3 <- dat3 %>% group_by(x,y,var) %>% mutate_at(vars(-group_cols()), ~./(fact[1]*fact[2]/100)) %>%
  ungroup()
summary(dat3)

# Plot data
library(ggplot2)
dat3 %>% filter(var == "Mixed forest") %>% 
  ggplot() + geom_tile(aes(x=x, y=y, fill=`2018`)) + coord_sf()
#ggsave("corine_lc_eur_1km.png", p, width=10, height=8, dpi=600)

# Save data
filename <- paste0("corine_lc_eur_", res)
assign(filename, dat3)
save(list=filename, file=paste0("data/corine_lc_eur_", res, ".rda"), compress="xz")
