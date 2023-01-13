#' ---
#' title: "Create UNESCO biosphere reserves data for Europe"
#' author: "RS-eco"
#' ---

library(tidyverse); library(sf); library(raster)

# Downloaded from https://zenodo.org/record/4905532

# Biosphere reserves details
BR_list_details <- vroom::vroom("https://zenodo.org/record/4905532/files/BR_list_details.csv?download=1")

# Biosphere reserves shapes
#download.file("https://zenodo.org/record/4905532/files/Vector_Data.zip?download=1", "extdata)
#unzip("extdata/Vector_Data.zip", exdir="extdata")
unesco_br_zones_eur <- sf::st_read("extdata/UNESCO_BR_Europe_EPSG-3035.shp")
unesco_br_zones_eur <- left_join(unesco_br_zones_eur, BR_list_details, by=c("BR_code"="BR_Code"))

# Load outline of Europe
load("data/europe.rda")
europe_laea <- sf::st_transform(europe, sf::st_crs(unesco_br_zones))

# Crop by extent of europe
library(sf)
unesco_br_zones_eur <- unesco_br_zones_eur %>% st_make_valid() %>% 
  st_crop(st_bbox(europe_laea))
gc()

# Plot shape
plot(sf::st_geometry(unesco_br_zones_eur))
plot(europe_laea, add=T)

# Simplify data to reduce file size
unesco_br_zones_eur <- st_simplify(unesco_br_zones_eur, dTolerance=5)

# Save to file
save(unesco_br_zones_eur, file="data/unesco_br_zones_eur.rda", compress="xz")
