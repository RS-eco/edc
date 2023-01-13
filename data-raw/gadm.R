#' Additional package to get shapefile of country outlines

rm(list=ls()); gc()
library(dplyr); library(sf)

# Load countries shapefile
#install.packages("rnaturalearthhires",
#                 repos = "http://packages.ropensci.org",
#                 type = "source")
data(countries10, package="rnaturalearthhires")
europe<- sf::st_as_sf(countries10) %>% filter(CONTINENT == "Europe")
head(europe)
europe$NAME
'%!in%' <- function(x,y)!('%in%'(x,y))
eur <- europe %>% 
  filter(NAME %!in% c("Russia"))
plot(sf::st_geometry(eur))

box <- rgeos::bbox2SP(n = 72, s = 30, w = -12, e = 42,
                      proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(st_geometry(st_as_sf(box)), add=T, border="red")

eur <- st_crop(eur, box)
plot(sf::st_geometry(eur))
unique(eur$NAME)

europe <- st_union(eur)

#' Save to file
save(europe, file="data/europe.rda", compress="xz")

# Select certain countries
eur_countries <- eur %>% 
  filter(NAME_EN %in% c("Austria", "Belgium", "Bulgaria", "Croatia", "Czech Republic", "Denmark", 
                     "Estonia", "France", "Germany", "Greece", "Hungary", "Ireland", "Italy", 
                     "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", "Malta", "Montenegro", 
                     "Netherlands", "Poland", "Portugal", "Romania", "Serbia", "Slovakia", 
                     "Slovenia", "Spain", "Switzerland", "United Kingdom"))
plot(sf::st_geometry(eur_countries))
save(eur_countries, file="data/eur_countries.rda", compress="xz")

eur_27 <- eur %>% 
  filter(NAME_EN %in% c("Austria", "Belgium", "Bulgaria", "Croatia", "Czech Republic", "Denmark", 
                        "Estonia", "France", "Germany", "Greece", "Hungary", "Ireland", "Italy", 
                        "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", "Malta", "Montenegro", 
                        "Netherlands", "Poland", "Portugal", "Romania", "Serbia", "Slovakia", 
                        "Slovenia", "Spain"))
plot(sf::st_geometry(eur_27))
save(eur_27, file="data/eur_27.rda", compress="xz")