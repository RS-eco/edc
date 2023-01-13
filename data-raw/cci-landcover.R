#' ---
#' title: "Create CCI data for Europe"
#' author: "RS-eco"
#' ---

rm(list=ls()); gc()

## CCI Land-cover data

# Information on the data can be found here:
# http://maps.elie.ucl.ac.be/CCI/viewer/download.php

# Data here is from Aidin Niamir: https://zenodo.org/record/3730469

# Specify file directory
filedir <- "C:/Users/Admin/Documents/CCI/"

# Load files
files <- list.files(filedir, pattern=".nc", full.names=T)
files

# Load data
library(terra)
dat <- terra::rast(files)
dat

# Crop by extent of europe
load("data/europe.rda")
r_cci_eur <- terra::mask(terra::crop(dat, terra::vect(sf::st_as_sf(europe))), terra::vect(sf::st_as_sf(europe)))
r_cci_eur
plot(r_cci_eur[[1]])
rm(dat); gc()

# Save to file
library(dplyr); library(tidyr)
cci_eur <- as.data.frame(r_cci_eur, xy=T) %>% 
  pivot_longer(names_to="var", values_to="fracCover", -c(x,y)) %>% 
  rowwise() %>% mutate(year = strsplit(var, split="[.]")[[1]][2],
                       class = strsplit(var, split="_")[[1]][2]) %>% 
  as.data.frame() %>% select(-var) %>%
  mutate(class = factor(class, levels = seq(0,220, by=10), 
            labels=c("No Data", "Cropland, rainfed", "Cropland, irrigated or post-flooding",
                     "Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)",
                     "Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)",
                     "Tree cover, broadleaved, evergreen, closed to open (>15%)",
                     "Tree cover, broadleaved, deciduous, closed to open (>15%)",
                     "Tree cover, needleleaved, evergreen, closed to open (>15%)",
                     "Tree cover, needleleaved, deciduous, closed to open (>15%)",
                     "Tree cover, mixed leaf type (broadleaved and needleleaved)",
                     "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)",
                     "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)", 
                     "Shrubland", "Grassland", "Lichens and mosses", 
                     "Sparse vegetation (tree, shrub, herbaceous cover) (<15%)", 
                     "Tree cover, flooded, fresh or brakish water", 
                     "Tree cover, flooded, saline water", 
                     "Shrub or herbaceous cover, flooded, fresh/saline/brakish water",
                     "Urban areas", "Bare areas", "Water bodies", "Permanent snow and ice")))
head(cci_eur)
cci_eur <- cci_eur %>% pivot_wider(values_from="fracCover", names_from="year")
head(cci_eur)
save(cci_eur, file="data/cci_eur.rda", compress="xz")
