## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", results="asis",
  warning=F, message=F, fig.width=9, fig.height=4
)

## -----------------------------------------------------------------------------
# Load packages
library(dplyr); library(magrittr); library(tidyr)
library(ggplot2); library(sf); library(patchwork)
library(scico)

# Load edc package
library(edc)

# Load shapefile of Europe
data("europe", package="edc")

## -----------------------------------------------------------------------------
# Load data of all Odonata species
data("iucn_odonata_eur_p11deg")

# Select data for Aeshna viridis and plot presence, origin and season information
p1 <- iucn_odonata_eur_p11deg %>% filter(binomial == "Aeshna viridis") %>%
 ggplot() + geom_tile(aes(x=x, y=y, fill=as.factor(presence))) +
  geom_sf(data=europe, fill=NA) + coord_sf() +
  scale_fill_discrete(name="Presence") + labs(x="", y="") + 
  theme_bw() + theme(legend.position = "bottom")
p2 <- iucn_odonata_eur_p11deg %>% filter(binomial == "Aeshna viridis") %>%
 ggplot() + geom_tile(aes(x=x, y=y, fill=as.factor(origin))) +
  geom_sf(data=europe, fill=NA) + coord_sf() +
  scale_fill_discrete(name="Origin") + labs(x="", y="") + 
  theme_bw() + theme(legend.position = "bottom")
p3 <- iucn_odonata_eur_p11deg %>% filter(binomial == "Aeshna viridis") %>%
 ggplot() + geom_tile(aes(x=x, y=y, fill=as.factor(presence))) +
  geom_sf(data=europe, fill=NA) + coord_sf() +
  scale_fill_discrete(name="Seasonal") + labs(x="", y="") + 
  theme_bw() + theme(legend.position = "bottom")
p1 + p2 + p3

## -----------------------------------------------------------------------------
# Load data for all bird species
data("birdlife_eur_p11deg")

# Select data for all Vanellus species and plot percentage cover data
birdlife_eur_p11deg %>% 
  filter(binomial %in% c("Vanellus gregarius", "Vanellus spinosus", "Vanellus vanellus")) %>%
  ggplot() + geom_tile(aes(x=x, y=y, fill=perc_present)) +
  facet_wrap(.~binomial) + 
  geom_sf(data=europe, fill=NA) + coord_sf() + 
   scale_fill_scico(palette = 'roma') + 
  labs(x="", y="") + theme_bw()

# Plot data of Vanellus spinosus
birdlife_eur_p11deg %>% 
  filter(binomial == "Vanellus spinosus") %>%
  ggplot() + geom_tile(aes(x=x, y=y, fill=perc_present)) +
  scale_fill_scico(palette = 'roma') + 
  labs(x="", y="") + theme_bw()
rm(list=ls()); gc()

