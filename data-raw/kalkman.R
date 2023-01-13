#' ---
#' title: "Read Kalkman et al. data"
#' author: "RS-eco"
#' ---

rm(list=ls()); gc()

library("tidyverse")
library("sf")
library("magrittr")

# Read Odonata data
odonata <- readxl::read_xlsx(paste(indir, "Kalkman et al._Eur_Dragonfly_AtlasData.xlsx",sep=""))
colnames(odonata)[2] <- "Species.name"
odonata$Presence <- 1
odonata$Latitude <- as.numeric(gsub(",", ".", odonata$Latitude))
odonata$Longitude <- as.numeric(gsub(",", ".", odonata$Longitude))
kalkman_odonata_eur <- odonata %>% distinct_all() %>% spread(Species.name, Presence) %>%
  replace(., is.na(.), 0)
save(kalkman_odonata_eur, file="data/kalkman_odonata_eur.rda", compress="xz")
