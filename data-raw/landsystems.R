#' ---
#' title: "Create landsystems data for Europe"
#' author: "RS-eco"
#' ---

library(tidyverse); library(sf); library(raster)

# Downloaded from https://dataverse.nl/dataset.xhtml?persistentId=doi:10.34894/XNC5KA

landsystem_dat <- raster::raster("extdata/EU_landSystem.tif")
landsystem_dat
NAvalue(landsystem_dat) <- 0
plot(landsystem_dat)

# Load outline of Europe
load("data/europe.rda")

# Crop by extent of europe
europe_laea <- sf::st_transform(europe, raster::projection(landsystem_dat))
landsystem_eur_1km <- raster::mask(raster::crop(landsystem_dat, as(europe_laea, "Spatial"), snap="out"), 
                               as(europe_laea, "Spatial"))
gc()
plot(landsystem_eur_1km)
plot(europe_laea, add=T)

landsystem_eur_1km <- as.data.frame(raster::rasterToPoints(landsystem_eur_1km))
colnames(landsystem_eur_1km) <- c("x", "y", "landsystem")
library(magrittr)
landsystem_eur_1km %<>% 
  mutate_at(vars(-c(x,y)), function(x) 
    factor(x, levels = c(21, 22, 23, 41, 42, 43, 51, 52, 53, 61, 62,
                         63, 31, 32, 71, 72, 74, 75, 731, 732, 733, 11,
                         12, 13, 80, 90), 
           labels=c("low-intensity settlement", "medium intensity settlement", "high intensity settlement",
                    "low-intensity forest", "medium-intensity forest", "high-intensity forest", "low-intensity grassland",
                    "medium-intensity grassland", "high-intensity grassland", "low-intensity cropland", 
                    "medium-intensity cropland", "high-intensity cropland", "extensive perm-crops", "intensive perm-crops",
                    "forest/shrubs and cropland mosaics", "forest/shrubs and grassland mosaics", 
                    "forest/shrubs and bare mosaics", "forest/shrubs and mixed agriculture mosaics",
                    "low-intensity agricultural mosaics", "medium-intensity agricultural mosaics", 
                    "high-intensity agricultural mosaics", "water body", "wetland", "glacier", "shrub",
                    "bare and rocks")))
head(landsystem_eur_1km)
save(landsystem_eur_1km, file="data/landsystem_eur_1km.rda", compress="xz")
rm(landsystem_eur_1km); gc()

# Layerize data
landsystem_w <- raster::layerize(landsystem_dat)
system_id <- sub("X","", names(landsystem_w))

# Aggregate data
landsystem_w <- lapply(1:raster::nlayers(landsystem_w), function(x){
  dat <- raster::aggregate(landsystem_w[[x]], fact=c(10, 10), fun=sum, expand=T, na.rm=T)
  raster::mask(dat, as(europe_laea, "Spatial"))
})
raster::plot(landsystem_w[[10]])

# Turn into data.frame & calculate percentage cover
landsystem_l <- lapply(1:length(landsystem_w), function(x){
  dat <- raster::rasterToPoints(landsystem_w[[x]]) %>% as.data.frame()
  colnames(dat) <- c("x", "y", "perc_cover")
  dat$landsystem <- system_id[x]
  return(dat)
})
landsystem_l <- bind_rows(landsystem_l)
colnames(landsystem_l)
head(landsystem_l)

# Re-structure data.frame & define categories
landsystem_perc_eur_10km <- landsystem_l %>% 
  mutate(landsystem = as.numeric(sub("X", "", landsystem))) %>% 
  mutate(landsystem = factor(landsystem, levels = c(21, 22, 23, 41, 42, 43, 51, 52, 53, 61, 62,
                                                    63, 31, 32, 71, 72, 74, 75, 731, 732, 733, 11,
                                                    12, 13, 80, 90), 
                             labels=c("low-intensity settlement", "medium intensity settlement", "high intensity settlement",
                                      "low-intensity forest", "medium-intensity forest", "high-intensity forest", "low-intensity grassland",
                                      "medium-intensity grassland", "high-intensity grassland", "low-intensity cropland", 
                                      "medium-intensity cropland", "high-intensity cropland", "extensive perm-crops", "intensive perm-crops",
                                      "forest/shrubs and cropland mosaics", "forest/shrubs and grassland mosaics", 
                                      "forest/shrubs and bare mosaics", "forest/shrubs and mixed agriculture mosaics",
                                      "low-intensity agricultural mosaics", "medium-intensity agricultural mosaics", 
                                      "high-intensity agricultural mosaics", "water body", "wetland", "glacier", "shrub",
                                      "bare and rocks"))) %>%
  tidyr::spread(landsystem, perc_cover) %>% ungroup() %>% 
  mutate_at(vars(-c(x,y)), list(~ tidyr::replace_na(., 0)))
head(landsystem_perc_eur_10km)
# Save to file
save(landsystem_perc_eur_10km, file="data/landsystem_perc_eur_10km.rda", compress="xz")
