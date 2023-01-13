#' ---
#' title: "Create biogeographic regions data for Europe"
#' author: "RS-eco"
#' ---

library(dplyr)

# Downloaded from https://zenodo.org/record/3766175#.XqbP5GgzZPZ
biogeographic_regions_eur <- vroom::vroom("https://bdj.pensoft.net/article/download/suppl/5759521/") %>%
  select(X,Y, CellCode, biogeographical_region)
head(biogeographic_regions_eur)

# Save to file
save(biogeographic_regions_eur, file="data/biogeographic_regions_eur.rda", compress="xz")