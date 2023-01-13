## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(comment = "#>", echo=T, warning=F, message=F)

## ----load_pkgs----------------------------------------------------------------
library(dplyr); library(sf); library(ggplot2)
library(tidyr); library(scico)

## -----------------------------------------------------------------------------
# Load shapefile of Europe
data("europe", package="edc")

## -----------------------------------------------------------------------------
# Named vector of clc colours
clc_cols <- c("Continuous urban fabric" = "#e6004d", "Discontinuous urban fabric" = "#ff0000", "Industrial or commercial units" = "#cc4df2", "Road and rail networks and associated land" = "#cc0000", "Port areas" = "#e6cccc", "Airports" = "#e6cce6", "Mineral extraction sites" = "#a600cc", "Dump sites" = "#a64d00", "Construction sites" = "#ff4dff", "Green urban areas" = "#ffa6ff", "Sport and leisure facilities" = "#ffe6ff", "Non-irrigated arable land" = "#ffffa8", "Permanently irrigated land" = "#ffff00", "Rice fields" = "#e6e600", "Vineyards" = "#e68000", "Fruit trees and berry plantations" = "#f2a64d", "Olive groves" = "#e6a600", "Pastures" = "#e6e64d", "Annual crops associated with permanent crops" = "#ffe6a6", "Complex cultivation patterns" = "#ffe64d", "Land principally occupied by agriculture with \n significant areas of natural vegetation" = "#e6cc4d", "Agro-forestry areas" = "#f2cca6", "Broad-leaved forest" = "#80ff00", "Coniferous forest" = "#00a600", "Mixed forest" = "#4dff00", "Natural grasslands" = "#ccf24d", "Moors and heathland" = "#a6ff80", "Sclerophyllous vegetation" = "#a6e64d", "Transitional woodland-shrub" = "#a6f200", "Beaches - dunes - sands" = "#e6e6e6", "Bare rocks" = "#cccccc", "Sparsely vegetated areas" = "#ccffcc", "Burnt areas" = "#000000", "Glaciers and perpetual snow" = "#a6e6cc", "Inland marshes" = "#a6a6ff", "Peat bogs" = "#4d4dff", "Salt marshes" = "#ccccff", "Salines" = "#e6e6ff", "Intertidal flats" = "#a6a6e6", "Water courses" = "#00ccf2", "Water bodies" = "#80f2e6", "Coastal lagoons" = "#00ffa6", "Estuaries" = "#a6ffe6", "Sea and ocean" = "#e6f2ff", "NODATA" = "#ffffff")

## -----------------------------------------------------------------------------
#data("corine_lc_eur_p11deg", package="edc")
data("corine_lc_eur_p44deg", package="edc")

# Plot data from long format
#corine_lc_eur_p11deg %>% 
#  filter(var %in% c("Pastures", "Broad-leaved forest", "Coniferous forest",
#                    "Mixed forest", "Natural grasslands")) %>% 
#  pivot_longer(names_to="year", values_to="value", -c(x,y,var)) %>%
#  ggplot() + geom_tile(aes(x=x, y=y, fill=value)) +
#  facet_grid(year~var) + scale_fill_scico(name="% Cover", palette="roma") + 
#  geom_sf(data=europe, fill=NA) + coord_sf() + theme_bw() + 
#  theme(strip.background = element_blank(), axis.title = element_blank())
#rm(corine_lc_eur_p11deg); invisible(gc())

# Plot data from long format
corine_lc_eur_p44deg %>% 
  filter(var %in% c("Pastures", "Broad-leaved forest", "Coniferous forest", 
                    "Mixed forest", "Natural grasslands")) %>% 
  pivot_longer(names_to="year", values_to="value", -c(x,y,var)) %>%
  ggplot() + geom_tile(aes(x=x, y=y, fill=value)) +
  facet_grid(year~var) + scale_fill_scico(name="% Cover", palette="roma") + 
  geom_sf(data=europe, fill=NA) + coord_sf() + theme_bw() + 
  theme(strip.background = element_blank(), axis.title = element_blank())

## -----------------------------------------------------------------------------
lc_eur <- corine_lc_eur_p44deg %>% dplyr::select(-c(`1990`, `2000`, `2006`,`2012`)) %>%
  pivot_wider(names_from="var", values_from=`2018`) %>%
  mutate_all(list(~replace_na(., 0)))
head(lc_eur)

lc_eur$forest <- lc_eur$`Agro-forestry areas` + lc_eur$`Coniferous forest` + 
  lc_eur$`Mixed forest` + lc_eur$`Broad-leaved forest`

# Plot individual land-cover percentage cover
lc_eur %>% ggplot() + geom_tile(aes(x=x, y=y, fill=forest)) + 
  scale_fill_scico(name="", palette="roma") + 
  coord_sf() + theme_bw() + labs(x="", y="") + 
  ggtitle("Forest cover (%)") + 
  theme(legend.text = element_text(size=10))
rm(lc_eur); invisible(gc())

## -----------------------------------------------------------------------------
# Calculate dominant land-cover class
dominant_clc <- corine_lc_eur_p44deg %>% group_by(x,y) %>% slice(which.max(`2018`)) %>% 
  dplyr::select(x,y,var, `2018`) %>% mutate(var = as.factor(var)) %>% drop_na()

# Select available land-cover classes
clc_cols_eur <- clc_cols[names(clc_cols) %in% dominant_clc$var]
dominant_clc %>% ggplot() +
  geom_tile(aes(x=x, y=y, fill=var)) + 
  scale_fill_manual(name="CLC2018", values=clc_cols_eur) + 
  coord_sf() + theme_bw() + theme(legend.position="bottom") + 
  labs(x="", y="")
rm(corine_lc_eur_p44deg, dominant_clc); invisible(gc())

