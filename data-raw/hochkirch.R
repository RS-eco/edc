rm(list=ls()); gc()

library("dtplyr")
library("dplyr")
library("tidyr")
library("sf")
library("magrittr")
library("stringr")

# Specify file directory
indir <- "/home/matt/t6p/group_hof/@BayKlif/data/"
#indir <- "Z://group_hof/@BayKlif/data/"

# Set working directory
setwd("/home/matt/Documents/edc")

# Check filelist
filelist <- list.files(paste0(indir, "Orthoptera_Europe_Hochkirch/"), pattern=".shp$")
filelist %<>% str_replace_all(".shp","") 
filelist %<>% str_replace_all("data/Orthoptera_Europe_Hochkirch/","") 

load("data/europe.rda")

# "Gomphocerippus_rufus"

# load species data
shapefile_list <- lapply(filelist, function(x){
  print(x)
  dat <- sf::read_sf(paste0(indir, "Orthoptera_Europe_Hochkirch/", x, ".shp", sep=""), type=6) %>%
    st_make_valid() %>% 
    select(any_of(c("BINOMIAL", "PRESENCE", "ORIGIN", "SEASONAL", "geometry")))
  #outline <- rgeos::gPolygonize(rgeos::gNode(as(as(world, "Spatial"), "SpatialLines")))
  #outline <- rgeos::gUnaryUnion(outline)
  
  # This drops invalid geometries, mostly range parts of Denmark!
  dat <- dat[st_is_valid(dat)==TRUE,]
  if(is.null(dat$BINOMIAL)){
    dat$BINOMIAL <- sub("_", " ", x)
  }
  if(is.null(dat$PRESENCE)){
    dat$PRESENCE <- NA
  }
  if(is.null(dat$ORIGIN)){
    dat$ORIGIN <- NA
  }
  if(is.null(dat$SEASONAL)){
    dat$SEASONAL <- NA
  }
  dat <- dat %>% filter(!is.na(st_dimension(geometry))) %>% st_cast(to="POLYGON") %>%
    group_by(BINOMIAL, PRESENCE, ORIGIN, SEASONAL) %>% summarise(); gc()
  dat2 <- dat %>% sf::st_make_valid() %>% sf::st_intersection(europe); gc()
  if(nrow(dat) == 0){
    dat <- NULL
  }
  return(dat)
}); gc()
shapefile_list <- Filter(Negate(is.null), shapefile_list); gc()

hochkirch_orthoptera_eur <- sf::st_as_sf(dplyr::bind_rows(shapefile_list, fill=T))
#hochkirch_orthoptera_eur <- hochkirch_orthoptera_eur %>%
#  group_by(BINOMIAL, PRESENCE, ORIGIN, SEASONAL) %>% summarise()
head(hochkirch_orthoptera_eur)
#hochkirch_orthoptera_eur %>% select(BINOMIAL, PRESENCE, ORIGIN, SEASONAL, geometry)
save(hochkirch_orthoptera_eur, file="data/hochkirch_orthoptera_eur.rda", compress="xz")
