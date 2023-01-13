#' ---
#' title: "Extract species-specific EBCC data"
#' ---

#' Data was downloaded from GBIF.org on the 15th of January 2021
#' The EBCC data file is located at t6p/group_hof/@BayKlif/data/EBCC-EBBA/
#' and the file is named: 0155958-200613084148143.csv

# Set working directory
filedir <- "/home/matt/t6p/group_hof/@BayKlif/data/EBCC-EBBA/"
#filedir <- "Z:///group_hof/@BayKlif/data/EBCC-EBBA/"

# Read data
library(vroom)
dat <- vroom::vroom(paste0(filedir, "/0155958-200613084148143.csv"))
head(dat)

library(dplyr); library(tidyr)
library(magrittr)

# Filter data by species name & select needed columns
ebcc_Vanellus_vanellus <- dat %>% filter(species == "Vanellus vanellus") %>% 
  dplyr::select(species, countryCode, occurrenceStatus, decimalLatitude, decimalLongitude) %>%
  mutate(occurrenceStatus = factor(occurrenceStatus, levels=c("ABSENT", "PRESENT"), labels=c(0,1))) %>%
  as.data.frame()

# Save data to file
save(ebcc_Vanellus_vanellus, file="data/ebcc_Vanellus_vanellus.rda", compress="xz")
rm(ebcc_Vanellus_vanellus); gc()

# Filter data by species name & select needed columns
ebcc_Phylloscopus_bonelli <- dat %>% filter(species == "Phylloscopus bonelli") %>% 
  dplyr::select(species, countryCode, occurrenceStatus, decimalLatitude, decimalLongitude) %>%
  mutate(occurrenceStatus = factor(occurrenceStatus, levels=c("ABSENT", "PRESENT"), labels=c(0,1))) %>%
  as.data.frame()

# Save data to file
save(ebcc_Phylloscopus_bonelli, file="data/ebcc_Phylloscopus_bonelli.rda", compress="xz")
rm(ebcc_Phylloscopus_bonelli); gc()

# Filter data by species name & select needed columns
ebcc_Sylvia_curruca <- dat %>% filter(species == "Sylvia curruca") %>% 
  dplyr::select(species, countryCode, occurrenceStatus, decimalLatitude, decimalLongitude) %>% 
  mutate(occurrenceStatus = factor(occurrenceStatus, levels=c("ABSENT", "PRESENT"), labels=c(0,1))) %>%
  as.data.frame()

# Save data to file
save(ebcc_Sylvia_curruca, file="data/ebcc_Sylvia_curruca.rda", compress="xz")
rm(ebcc_Sylvia_curruca); gc()

# Filter data by species name & select needed columns
ebcc_Turdus_torquatus <- dat %>% filter(species == "Turdus torquatus") %>% 
  dplyr::select(species, countryCode, occurrenceStatus, decimalLatitude, decimalLongitude) %>%
  mutate(occurrenceStatus = factor(occurrenceStatus, levels=c("ABSENT", "PRESENT"), labels=c(0,1))) %>%
  as.data.frame()

# Save data to file
save(ebcc_Turdus_torquatus, file="data/ebcc_Turdus_torquatus.rda", compress="xz")
rm(ebcc_Turdus_torquatus); gc()

# Select needed columns
dat %<>% dplyr::select(species, occurrenceStatus, lat=decimalLatitude, lon=decimalLongitude)
dat$occurrenceStatus <- factor(dat$occurrenceStatus, levels=c("ABSENT", "PRESENT"), 
                               labels=c(0,1))
dat$occurrenceStatus <- as.numeric(dat$occurrenceStatus)
ebcc_eur <- dat %>% as.data.frame() %>% tidyr::drop_na() %>% 
  tidyr::pivot_wider(names_from=species, values_from=occurrenceStatus, values_fn = sum) %>%
  mutate(across(-c(lon,lat), ~if_else(. > 1, 1, .))) %>% 
  mutate_at(vars(-c(lon,lat)), replace_na, 0)
head(ebcc_eur)

# Save data to file
save(ebcc_eur, file="data/ebcc_eur.rda", compress="xz")

###########################

