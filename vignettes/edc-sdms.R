## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, warning=F, message=F,
                      fig.width=8, fig.height=6)

## -----------------------------------------------------------------------------
rm(list=ls()); gc()

# Load packages
library(dplyr); library(tidyr)
library(ggplot2); library(scico)

# Load edc package
library(edc)

# Load shapefile of Europe
data("europe", package="edc")

## -----------------------------------------------------------------------------
# Install mecofun package if not available
if(!"mecofun" %in% installed.packages()[,"Package"]){
  remotes::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")
}

# Load mecofun package
library(mecofun)
#' This package includes the following functions:
#' predictSDM, crossvalSDM, evalSDM, TSS, expl_deviance, inflated_response
#' eo_mask, partial_response, range_size, range_centre, select07, select07_cv

## -----------------------------------------------------------------------------
# Load bird species range data for all grid cells
data("ebcc_eur")

# Load climate data
data("cordex_bioclim_eur_p44deg")

# Select for current conditions & calculate ensemble mean
curclim <- cordex_bioclim_eur_p44deg %>% filter(time_frame == "1991-2020") %>%
  group_by(x,y) %>% summarise_at(vars(bio1:bio19), ~mean(.,na.rm=T))

# Add landcover data to curclim
data("corine_lc_eur_p44deg")
landcover_eur <- corine_lc_eur_p44deg %>% dplyr::select(x,y,var,`2018`) %>% 
  pivot_wider(names_from="var", values_from=`2018`)

curclim <- raster::rasterFromXYZ(curclim, crs="+init=EPSG:4326")
landcover_eur <- raster::rasterFromXYZ(landcover_eur, crs="+init=EPSG:4326")
curenv <- raster::stack(curclim, landcover_eur) %>% raster::rasterToPoints() %>%
  as.data.frame() %>% mutate_at(vars(-c(x,y)), replace_na, 0); rm(curclim, landcover_eur)

# Select for future conditions and calculate ensemble mean across GCMs
futclim <- cordex_bioclim_eur_p44deg %>% filter(time_frame != "1991-2020") %>%
  group_by(x,y,time_frame, rcp) %>% summarise_at(vars(bio1:bio19), ~mean(.,na.rm=T))

## -----------------------------------------------------------------------------
# plot species using ggplot
ebcc_eur %>% 
  mutate(`Nucifraga caryocatactes`= as.factor(`Nucifraga caryocatactes`)) %>% 
  ggplot() + geom_point(aes(x=lon, y=lat, color=`Nucifraga caryocatactes`), 
                        shape=15, size=1) + scale_color_manual(values=c("grey80", "blue")) + 
  coord_equal() + theme_bw() + labs(x="", y="") + 
  ggtitle("Nucifraga caryocatactes current range") +
  theme(legend.position = "none")

# plot environmental variables

#define two colour scales for annual mean temperature and annual precipitation
tempcol <- scale_fill_scico(name = "°C", palette="roma", direction=-1,
                            limits = c(min(curenv$bio1),max(cordex_bioclim_eur_p44deg$bio1)))
preccol <- scale_fill_scico(name = "mm", palette="roma",
                            limits = c(min(curenv$bio12),max(cordex_bioclim_eur_p44deg$bio12)))

curenv %>% # which dataset to use
  ggplot() + # start plotting
  # type of plot, define x and y axis, define fill the environmental variable
  geom_tile(aes(x=x, y=y, fill=bio1)) + 
  scale_fill_scico(name="°C", palette="roma", direction=-1) +
  coord_sf() + theme_bw() + labs(x="", y="") + # set x/y to equal, use custom map-theme
  ggtitle("Annual Mean Temperature") + # set plot title
  theme(legend.text = element_text(size=10)) # set font size for the legend

#different colour scales (see above):

#bio1: annual mean temperature

curenv %>% ggplot() + geom_tile(aes(x=x, y=y, fill=bio1)) + 
  tempcol + coord_sf() + theme_bw() + labs(x="", y="") + 
  ggtitle("Annual Mean Temperature") + 
  theme(legend.position = "right",
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent"))

#bio12: annual precipitation

curenv %>% ggplot() + geom_tile(aes(x=x, y=y, fill=bio12)) + 
  preccol + coord_sf() + theme_bw() + labs(x="", y="") + 
  ggtitle("Annual Precipitation") + 
  theme(legend.position = "right",
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent"))

#bio1

futclim %>% filter(time_frame == "2041-2070", rcp == "rcp26") %>% 
  ggplot() + geom_tile(aes(x=x, y=y, fill=bio1)) + 
  tempcol + coord_sf() + theme_bw() + labs(x="", y="") +  
  ggtitle("Annual Mean Temperature RCP2.6, 2055") + 
  theme(legend.position = "none",
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent"))

#bio12

futclim %>% filter(time_frame == "2041-2070", rcp == "rcp26") %>% 
  ggplot() + geom_tile(aes(x=x, y=y, fill=bio12)) + 
  preccol + coord_sf() + theme_bw() + labs(x="", y="") + 
  ggtitle("Annual Precipitation RCP2.6, 2055") + 
  theme(legend.position = "none",
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill="transparent"))

## -----------------------------------------------------------------------------
# Turn bird data into long format
birds_eur <- ebcc_eur %>% pivot_longer(names_to="species", values_to="occurrenceStatus", -c(lon,lat)); rm(ebcc_eur); gc()

# Extract data for certain species
birds_eur <- birds_eur
#transfer bird records to 0.44 degrees grid
birds_eur %<>% filter(occurrenceStatus == 1) %>% as.data.frame()
birds <- letsR::lets.presab.points(xy=birds_eur[,c("lon", "lat")], species=birds_eur$species, 
                                   xmn=min(birds_eur$lon), xmx=max(birds_eur$lon), 
                                   ymn=min(birds_eur$lat), ymx=max(birds_eur$lat), 
                                   resol=1, show.matrix = T)
birds <- as.data.frame(birds) %>% dplyr::select(c("Longitude(x)", "Latitude(y)", "Alauda arvensis", "Perdix perdix", "Saxicola rubetra", "Nucifraga caryocatactes", "Vanellus vanellus","Phylloscopus bonelli", "Sylvia curruca", "Numenius arquata", "Ficedula albicollis"))
colnames(birds) <- c("lon", "lat", "Alauda_arvensis", "Perdix_perdix", "Saxicola_rubetra", 
                     "Nucifraga_caryocatactes", "Vanellus_vanellus", "Phylloscopus_bonelli", 
                     "Sylvia_curruca", "Numenius_arquata", "Ficedula_albicollis")
birds <- birds %>% replace_na(list(Alauda_arvensis = 0, Perdix_perdix = 0, Saxicola_rubetra = 0, Vanellus_vanellus = 0, Phylloscopus_bonelli = 0, Sylvia_curruca = 0, Numenius_arquata = 0, Ficedula_albicollis = 0))

curenv_r <- raster::rasterFromXYZ(curenv) #turn climate data into raster
birds_r <- raster::rasterFromXYZ(birds) #turn species data (here: odonata) into raster
birds_r <- raster::mask(raster::resample(birds_r, curenv_r), curenv_r[[1]])
birds_r[birds_r > 0] <- 1
#raster::plot(birds_r)

spec_env <- raster::stack(birds_r, curenv_r) #combine Numenius.arquata and climate raster
spec_env <- as.data.frame(raster::rasterToPoints(spec_env)) %>% drop_na()

# 2 explanatory variables:
myExpl2 <- c("bio1", "bio12")

# Automatically create model formula from variables
(formExpl2 <- as.formula(paste("Nucifraga_caryocatactes ~ ", 
                               paste(myExpl2, collapse="+"),sep = "")))

# Fit Generalized Linear Model (GLM) in the simplest form
glm_NCaryo_2 <- glm(formExpl2, family='binomial', data=spec_env)

# You would need to specify polynomials and interactions manually in the formula:

# Fit Generalized Linear Model (GLM) in a bit more complex form (with polynomials)
glm_NCaryo_2_Pol <- glm(Nucifraga_caryocatactes ~ bio1 + I(bio1^2) + bio12 + I(bio12^2), 
                        family='binomial', data=spec_env)

# Fit Generalized Linear Model (GLM) in an even more complex form (with polynomials and interactions)
glm_NCaryo_2_PolInt <- glm(Nucifraga_caryocatactes ~ bio1 + I(bio1^2) + bio12 + I(bio12^2) + 
                             bio1*bio12, family='binomial', data=spec_env)

# check model summary ----
(sum_glm_NCaryo_2 <- summary(glm_NCaryo_2))
(sum_glm_NCaryo_2_PolInt <- summary(glm_NCaryo_2_PolInt))
(sum_glm_NCaryo_2_Pol <- summary(glm_NCaryo_2_Pol))

#par(mfrow=c(1,1))
#plot(glm_NCaryo_2)

## -----------------------------------------------------------------------------
# basic plots of occurrence vs. explanatory variable
par(mfrow=c(2,1))
plot(spec_env$bio1, spec_env$Nucifraga_caryocatactes, ylab="", xlab="bio1")
plot(spec_env$bio12, spec_env$Nucifraga_caryocatactes, ylab="", xlab="bio12")

# check partial response curves
par(mfrow=c(3,2))
partial_response(glm_NCaryo_2, predictors = spec_env[,myExpl2])
partial_response(glm_NCaryo_2_Pol, predictors = spec_env[,myExpl2]) 
partial_response(glm_NCaryo_2_PolInt, predictors = spec_env[,myExpl2])

#' This is needed for getting TSS, AUC and Kappa values

# Make cross-validated predictions for GLM:
crosspred_glm_NCaryo_2 <- crossvalSDM(glm_NCaryo_2, kfold=5, 
                                      traindat= spec_env, colname_pred=myExpl2, 
                                      colname_species = "Nucifraga_caryocatactes")
crosspred_glm_NCaryo_2_Pol <- crossvalSDM(glm_NCaryo_2_Pol, kfold=5, 
                                          traindat= spec_env, colname_pred=myExpl2, 
                                          colname_species = "Nucifraga_caryocatactes")

# Assess cross-validated model performance
(eval_glm_NCaryo_2 <- evalSDM(observation = spec_env$Nucifraga_caryocatactes, 
                              predictions = crosspred_glm_NCaryo_2))
(eval_glm_NCaryo_2_Pol <- evalSDM(observation = spec_env$Nucifraga_caryocatactes, 
                                  predictions = crosspred_glm_NCaryo_2_Pol))

# check variable importances
(glm_imp <- caret::varImp(glm_NCaryo_2, scale=T))
par(mfrow=c(1,1))
barplot((glm_imp$Overall/sum(glm_imp$Overall)*100)[2:1], 
        names.arg=rownames(glm_imp)[2:1], horiz=T,
        main="GLM", xlab="Relative influence")

## -----------------------------------------------------------------------------
# Make predictions to current climate:
spec_env$pred_glm_NCaryo_2 <- predictSDM(glm_NCaryo_2, spec_env)

# Make binary/threshholded predictions:
spec_env$bin_glm_NCaryo_2 <- ifelse(spec_env$pred_glm_NCaryo_2 > eval_glm_NCaryo_2$thresh, 1, 0)

par(mfrow=c(1,1), mar=c(5,5,4,1))
boxplot(spec_env$pred_glm_NCaryo_2 ~ spec_env$Nucifraga_caryocatactes, las=1, 
        xlab="Aktuelle Verbreitung", ylab="Vorkommenswahrscheinlichkeit", cex.lab=2,
        col="dodgerblue4", cex.axis=1.5, main="GLM-Modell mit 2 Klimavariablen", cex.main=2)

#---------------------------------------------------
#' ## plot species using ggplot 
#---------------------------------------------------

# plot histogramm

spec_env %>% ggplot() + geom_histogram(aes(x=pred_glm_NCaryo_2), col="grey0", alpha=0.2) +
  labs(y="Number of grid cells", x="Occurrence probability",
       title="Distribution of N. caryocatactes occurrence probability under current climate") + 
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey") + theme_classic()

# plot maps

# probability map
(plot_curclim_GLM <- spec_env %>% ggplot() + 
    geom_tile(aes(x=x, y=y, fill=pred_glm_NCaryo_2)) + 
    scale_fill_scico(name="Probability of\noccurrence", palette="roma", direction=-1)  +  
    coord_sf() + theme_bw() + labs(x="", y="") + 
    ggtitle("GLM current climate") + 
    theme(legend.text = element_text(size=10)))

# binary map
(plot_curclim_GLM_bin <- 
    spec_env %>%
    ggplot() + geom_tile(aes(x=x, y=y, fill=as.factor(bin_glm_NCaryo_2))) + 
    scale_fill_manual(name="Occurrence", values=c("grey80", "blue")) + 
    coord_sf() + theme_bw() + labs(x="", y="") + 
    ggtitle("GLM current climate") + 
    theme(legend.text = element_text(size=10)))

## -----------------------------------------------------------------------------
# Assess novel environments in future climate layer:

# Values of 1 in the eo.mask will indicate novel environmental conditions
futclim$eo_mask <- eo_mask(curenv[,myExpl2],futclim[,myExpl2])
futclim %>% # which dataset to use
  ggplot() +              # start plotting
  # type of plot, define x axis, y axis is automatically a count, define color and transparency
  geom_tile(aes(x=x, y=y, fill=as.factor(eo_mask))) +
  scale_fill_manual(name="", values=c("grey80", "blue")) + 
  coord_sf() + theme_bw() + labs(x="", y="") + # set x/y to equal, use custom map-theme
  ggtitle("Environmental novelty") +       # set plot title
  theme(legend.text = element_text(size=10)) # set font size for the legend

# Make predictions to futclim
futclim$pred_glm_NCaryo_2 <- predictSDM(glm_NCaryo_2, futclim)

# Make binary/threshholded predictions:
futclim$bin_glm_NCaryo_2 <- ifelse(futclim$pred_glm_NCaryo_2 > eval_glm_NCaryo_2$thresh, 1, 0)

# => using this framework you can proceed with other future climate data

#---------------------------------------------------
# plot species using ggplot 
#---------------------------------------------------

# plot histogramm
futclim %>% # which dataset to use
  ggplot() +              # start plotting
  # type of plot, define x axis, y axis is automatically a count, define color and transparency
  geom_histogram(aes(x=pred_glm_NCaryo_2), col="grey0", alpha=0.2) +
  labs(y="Number of grid cells", x="Occurrence probability",        # add axis labels
       title="Distribution of N. caryocatactes occurrence probability under CC2670 GLM") +  
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey") + # add a horizontal line
  theme_classic()

# plot maps

# probability map
(plot_CC6070_GLM <- 
    futclim %>% # which dataset to use
    ggplot() + geom_tile(aes(x=x, y=y, fill=pred_glm_NCaryo_2)) + 
    scale_fill_scico(name="Probability of\noccurrence", palette = "roma", direction=-1) +
    coord_sf() + theme_bw() + labs(x="", y="") + 
    ggtitle("GLM CC6070") +       # set plot title
    theme(legend.text = element_text(size=10))) # set font size for the legend

# binary map
(plot_CC6070_GLM_bin <- 
    futclim %>% # which dataset to use
    ggplot() +  geom_tile(aes(x=x, y=y, fill=as.factor(bin_glm_NCaryo_2))) + 
    scale_fill_manual(name="Occurrence", values=c("grey80", "blue")) + 
    coord_sf() + theme_bw() + labs(x="", y="") + 
    ggtitle("GLM CC6070") +       # set plot title
    theme(legend.text = element_text(size=10))) # set font size for the legend
rm(list=ls()); gc()

