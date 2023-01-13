## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  comment = "#>", echo=T, warning=F, message=F,
  fig.width=9, fig.height=6
)

## ----load_pkgs----------------------------------------------------------------
library(dplyr); library(ggplot2); library(sf)

## -----------------------------------------------------------------------------
# Load shapefile of Europe
data("europe", package="edc")

## ----eur_fine, fig.width=8, fig.height=9--------------------------------------
data("srtm_csi_eur_p11deg", package="edc")

# Plot individual variable from long format
srtm_csi_eur_p11deg %>% ggplot() + 
  geom_tile(aes (x=x, y=y, fill=alt_mean)) +
  geom_sf(data=europe, fill="NA") + 
  scale_fill_gradientn(name="Elevation (m)", colors=terrain.colors(255)) + 
  coord_sf(expand=F) + theme_bw() + labs(x="", y="")
rm(srtm_csi_eur_p11deg); invisible(gc())

## ----eur_rough, fig.width=8, fig.height=9-------------------------------------
data("srtm_csi_eur_p44deg", package="edc")

# Plot individual variable from long format
srtm_csi_eur_p44deg %>% ggplot() + 
  geom_tile(aes (x=x, y=y, fill=alt_mean)) +
  geom_sf(data=europe, fill="NA") + 
  scale_fill_gradientn(name="Elevation (m)", colors=terrain.colors(255)) + 
  coord_sf(expand=F) + theme_bw() + labs(x="", y="")
rm(srtm_csi_eur_p44deg); invisible(gc())

