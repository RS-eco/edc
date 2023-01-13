## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", results="asis",
  warning=F, message=F, fig.width=10, fig.height=8
)

## -----------------------------------------------------------------------------
# Load packages
library(dplyr); library(magrittr); library(tidyr)
library(ggplot2); library(patchwork)
library(sf); library(sp)

# Load edc package
library(edc)

# Load shapefile of Europe
data("europe", package="edc")

## -----------------------------------------------------------------------------
data("hochkirch_orthoptera_eur")
#head(hochkirch_orthoptera_eur)

## -----------------------------------------------------------------------------
data("kalkman_odonata_eur")
#head(kalkman_odonata_eur)

## -----------------------------------------------------------------------------
# Read Euro-Cordex climate data
data("cordex_bioclim_eur_p44deg")

# Split data by scenario and timeframe to be able to use rasterFromXYZ
groupkeys <- cordex_bioclim_eur_p44deg %>% 
  tidyr::unite(col="gcm_rcm_rcp", c("gcm", "ensemble", "rcm", "rs", "rcp"), sep="_") %>% 
  group_by(gcm_rcm_rcp, time_frame) %>% group_keys()
cordex_bioclim_eur_p44deg <- cordex_bioclim_eur_p44deg %>% 
  tidyr::unite(col="gcm_rcm_rcp", c("gcm", "ensemble", "rcm", "rs", "rcp"), sep="_") %>% 
  group_split(gcm_rcm_rcp, time_frame)

#' The rasterFromXYZ function requires that your data.frame is structured in a very specific way,
#' i.e. first column is x, second column is y and the remaining columns can only be numeric variables.

# Turn data into raster
r_bioclim <- lapply(cordex_bioclim_eur_p44deg, raster::rasterFromXYZ)
#raster::plot(r_bioclim[[1]][[1]])
#plot(sf::st_geometry(europe), add=T)

#r_bioclim[[1]]

#' Data has no CRS

# Define CRS (WGS84)
r_bioclim <- lapply(r_bioclim, function(x){
  raster::crs(x) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  return(x)
})
#r_bioclim[[1]]

#' Now data has correct CRS

## -----------------------------------------------------------------------------
# Extract climate data for Acheta hispanicus
ortho_bioclim <- lapply(1:length(r_bioclim), function(x){
  raster::extract(x=r_bioclim[[x]], y=as(hochkirch_orthoptera_eur %>% 
                    filter(BINOMIAL == "Acheta hispanicus"), "Spatial")) %>% as.data.frame() %>% 
    mutate(gcm_rcm_rcp=groupkeys$gcm_rcm_rcp[x], time_frame=groupkeys$time_frame[x])
})
ortho_bioclim <- bind_rows(ortho_bioclim) %>% drop_na()
#head(ortho_bioclim)
#summary(ortho_bioclim)

ortho_bioclim %>% filter(gcm_rcm_rcp == "ICHEC-EC-EARTH_r1i1p1_KNMI-RACMO22E_v1_rcp45") %>% 
  select(time_frame, bio1:bio19) %>%
  ggplot() + geom_violin(aes(x=time_frame, y=bio1)) + 
  labs(x="Presence", y="bio1", title="Acheta hispanicus")

## -----------------------------------------------------------------------------
# Turn Odonata data into Spatial Points
sf_odonata <- sf::st_as_sf(kalkman_odonata_eur, coords = c("Longitude", "Latitude"), 
                           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Plot data of Ophiogomphus cecilia
plot(sf::st_geometry(sf_odonata %>% select("Ophiogomphus cecilia")))
plot(sf::st_geometry(europe), add=T)

# Extract climate data for Odonata data
odo_bioclim <- lapply(1:length(r_bioclim), function(x){
  raster::extract(r_bioclim[[x]], sf_odonata, sp=T) %>% 
    as.data.frame() %>% select(-c(Country, `MGRS.WGS84`)) %>%
    mutate(gcm_rcm_rcp=groupkeys$gcm_rcm_rcp[x], time_frame=groupkeys$time_frame[x])
})
odo_bioclim <- bind_rows(odo_bioclim)
#head(odo_bioclim)
#summary(odo_bioclim)

p1 <- odo_bioclim %>% filter(gcm_rcm_rcp == "ICHEC-EC-EARTH_r1i1p1_KNMI-RACMO22E_v1_rcp45", 
                             time_frame== "1991-2020", Period ==">1990") %>% 
  select(coords.x1, coords.x2, bio1:bio19, Aeshna.cyanea) %>%
  ggplot() + geom_violin(aes(x=as.factor(Aeshna.cyanea), y=bio1)) + 
  labs(x="Presence", y="bio1", title="Aeshna cyanea") + theme_bw()
p2 <- odo_bioclim %>% filter(gcm_rcm_rcp == "ICHEC-EC-EARTH_r1i1p1_KNMI-RACMO22E_v1_rcp45", 
                             time_frame== "1991-2020", Period ==">1990") %>% 
  select(coords.x1, coords.x2, bio1:bio19, Aeshna.cyanea) %>%
  ggplot() + geom_violin(aes(x=as.factor(Aeshna.cyanea), y=bio1)) + 
  labs(x="Presence", y="bio12", title="Aeshna cyanea") + theme_bw()
p1 + p2
rm(list=ls()); gc()

## ---- eval=F------------------------------------------------------------------
#  
#  odonata2 <- odonata
#  odonata2$Longitude <- odonata2$Longitude + 0.44
#  sf_odonata2 <- sf::st_as_sf(odonata2, coords = c("Longitude", "Latitude"),
#                              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#  
#  odonata3 <- odonata
#  odonata3$Longitude <- odonata3$Longitude - 0.44
#  sf_odonata3 <- sf::st_as_sf(odonata3, coords = c("Longitude", "Latitude"),
#                              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#  
#  odonata4 <- odonata
#  odonata4$Latitude <- odonata4$Latitude - 0.44
#  sf_odonata4 <- sf::st_as_sf(odonata4, coords = c("Longitude", "Latitude"),
#                              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#  
#  odonata5 <- odonata
#  odonata5$Latitude <- odonata5$Latitude + 0.44
#  sf_odonata5 <- sf::st_as_sf(odonata5, coords = c("Longitude", "Latitude"),
#                              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#  
#  # Extract climate data for Odonata data
#  cell_no <- raster::extract(r_bioclim, sf_odonata, cellnumbers=T, df=T)
#  odonata$x <- raster::xFromCell(r_bioclim, cell_no$cells)
#  odonata$y <- raster::yFromCell(r_bioclim, cell_no$cells)
#  cell_no2 <- raster::extract(r_bioclim, sf_odonata2, cellnumbers=T, df=T)
#  odonata2$x <- raster::xFromCell(r_bioclim, cell_no2$cells)
#  odonata2$y <- raster::yFromCell(r_bioclim, cell_no2$cells)
#  cell_no3 <- raster::extract(r_bioclim, sf_odonata3, cellnumbers=T, df=T)
#  odonata3$x <- raster::xFromCell(r_bioclim, cell_no3$cells)
#  odonata3$y <- raster::yFromCell(r_bioclim, cell_no3$cells)
#  cell_no4 <- raster::extract(r_bioclim, sf_odonata4, cellnumbers=T, df=T)
#  odonata4$x <- raster::xFromCell(r_bioclim, cell_no4$cells)
#  odonata4$y <- raster::yFromCell(r_bioclim, cell_no4$cells)
#  cell_no5 <- raster::extract(r_bioclim, sf_odonata5, cellnumbers=T, df=T)
#  odonata5$x <- raster::xFromCell(r_bioclim, cell_no5$cells)
#  odonata5$y <- raster::yFromCell(r_bioclim, cell_no5$cells)
#  
#  odonata %<>% dplyr::select(-c(Country, Latitude, Longitude, `MGRS WGS84`, Period)) %>% right_join(r_bioclim_df[,c("x","y")]) %>%
#    mutate_at(vars(-c(x,y)), list(~ replace_na(., 0)))
#  odonata2 %<>% dplyr::select(-c(Country, Latitude, Longitude, `MGRS WGS84`, Period)) %>% right_join(r_bioclim_df[,c("x","y")]) %>%
#    mutate_at(vars(-c(x,y)), list(~ replace_na(., 0)))
#  odonata3 %<>% dplyr::select(-c(Country, Latitude, Longitude, `MGRS WGS84`, Period)) %>% right_join(r_bioclim_df[,c("x","y")]) %>%
#    mutate_at(vars(-c(x,y)), list(~ replace_na(., 0)))
#  odonata4 %<>% dplyr::select(-c(Country, Latitude, Longitude, `MGRS WGS84`, Period)) %>% right_join(r_bioclim_df[,c("x","y")]) %>%
#    mutate_at(vars(-c(x,y)), list(~ replace_na(., 0)))
#  odonata5 %<>% dplyr::select(-c(Country, Latitude, Longitude, `MGRS WGS84`, Period)) %>% right_join(r_bioclim_df[,c("x","y")]) %>%
#    mutate_at(vars(-c(x,y)), list(~ replace_na(., 0)))
#  
#  odonata <- odonata %>% gather(species, presence, -c(x,y)) %>% bind_rows(odonata2 %>% gather(species, presence, -c(x,y))) %>%
#    bind_rows(odonata3 %>% gather(species, presence, -c(x,y))) %>%
#    bind_rows(odonata4 %>% gather(species, presence, -c(x,y))) %>%
#    bind_rows(odonata5 %>% gather(species, presence, -c(x,y))) %>%
#    group_by(x,y, species) %>% summarise(presence=mean(as.numeric(presence), na.rm=T))
#  odonata$presence[odonata$presence > 0] <- 1
#  odonata <- odonata %>% spread(species, presence)
#  head(odonata)
#  odonata %>% ggplot() + geom_tile(aes(x=x, y=y, fill=`Aeshna cyanea`))

