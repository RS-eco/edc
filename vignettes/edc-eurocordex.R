## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  results="asis", cache=F,
  echo=F, warning=F, message=F,
  fig.width=10, fig.height=8
)

## ----load_pkgs----------------------------------------------------------------
rm(list=ls()); invisible(gc())
library(dplyr); library(magrittr); library(tidyr); library(ggplot2)
library(ggthemes); library(patchwork); library(zoo)

# Define standard colour scheme
whitered <- colorRampPalette(c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'))(255)
#pie(rep(1,9),col=c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'))

whiteblue <- colorRampPalette(c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081'))(255)
#pie(rep(1,9),col=c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081'))

bluewhitered <- colorRampPalette(c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e"))(255)
#pie(rep(1,7),col=c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e"))

## -----------------------------------------------------------------------------
# Load shapefile of Europe
data("europe", package="edc")

## -----------------------------------------------------------------------------
#data("cordex_bioclim_eur", package="edc")
data("cordex_bioclim_eur_p44deg", package="edc")

## ---- fig.width=10, fig.height=10---------------------------------------------
tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>% mutate(var = factor(var, levels=paste0("bio", 1:19))) %>% 
  ggplot() + geom_violin(aes(x=time_frame, y=value, fill=rcp)) + 
  facet_wrap(.~var, scales="free_y", strip.position = "left", ncol=4) + 
  theme_few() + theme(strip.background = element_blank(),
                       legend.position=c(0.9,0.05), strip.placement = "outside",
                       axis.text.x = element_text(angle=45, vjust=0.5)) + 
  scale_fill_manual(name="", values=c("#0099B4FF", "#ADB6B6FF", "#AD002AFF")) + labs(x="", y="")

## ---- fig.width=10, fig.height=10---------------------------------------------
tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  filter(var == "bio1") %>% ggplot() + geom_violin(aes(x=time_frame, y=value, fill=rcp)) + 
  facet_wrap(.~gcm+rcm+ensemble, ncol=5) + 
  theme_few() + theme(legend.position=c(0.9,0.2), axis.text.x = element_text(angle=45, vjust=0.5)) + 
  scale_fill_manual(name="", values=c("#0099B4FF", "#ADB6B6FF", "#AD002AFF")) + labs(x="", y="")

## ---- fig.width=10, fig.height=10---------------------------------------------
tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  filter(var == "bio12") %>% ggplot() + geom_violin(aes(x=time_frame, y=value, fill=rcp)) + 
  facet_wrap(.~gcm+rcm+ensemble, ncol=5) + 
  theme_few() + theme(legend.position=c(0.9,0.2), axis.text.x = element_text(angle=45, vjust=0.5)) + 
  scale_fill_manual(name="", values=c("#0099B4FF", "#ADB6B6FF", "#AD002AFF")) + labs(x="", y="")

## -----------------------------------------------------------------------------
dat <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  mutate(var = factor(var, levels=paste0("bio", 1:19))) %>% 
  filter(var == "bio1") %>%
  group_by(x, y, time_frame, rcp, var) %>% summarise(value=mean(value))
col_val <- scales::rescale(quantile(dat$value, probs=seq(0,1,0.12)))
lim <- c(min(dat$value, na.rm=T), max(dat$value, na.rm=T))
dat %>% ggplot() + geom_tile(aes(x=x, y=y, fill=value)) + 
  facet_grid(rcp~time_frame, switch = "y") + geom_sf(data=europe, fill="NA") + 
  scale_fill_gradientn(name="bio1", colours=whitered, limits=lim, values=col_val) + 
  coord_sf() + theme_few() + theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + 
  labs(x="", y="")

## -----------------------------------------------------------------------------
dat <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  mutate(var = factor(var, levels=paste0("bio", 1:19))) %>% 
  filter(var == "bio4") %>%
  group_by(x, y, time_frame, rcp, var) %>% summarise(value=mean(value))
col_val <- scales::rescale(quantile(dat$value, probs=seq(0,1,0.12)))
lim <- c(min(dat$value, na.rm=T), max(dat$value, na.rm=T))
dat %>% ggplot() + geom_tile(aes(x=x, y=y, fill=value)) + 
  facet_grid(rcp~time_frame, switch = "y") + geom_sf(data=europe, fill="NA") + 
  scale_fill_gradientn(name="bio4", colours=whitered, limits=lim, values=col_val) + 
  coord_sf() + theme_few() + theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + 
  labs(x="", y="")

## -----------------------------------------------------------------------------
dat <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  mutate(var = factor(var, levels=paste0("bio", 1:19))) %>% 
  filter(var == "bio5") %>%
  group_by(x, y, time_frame, rcp, var) %>% summarise(value=mean(value))
col_val <- scales::rescale(quantile(dat$value, probs=seq(0,1,0.12)))
lim <- c(min(dat$value, na.rm=T), max(dat$value, na.rm=T))

dat %>% ggplot() + geom_tile(aes(x=x, y=y, fill=value)) + geom_sf(data=europe, fill="NA") + 
  facet_grid(rcp~time_frame, switch = "y") + scale_fill_gradientn(name="bio5", colours=whitered, limits=lim, values=col_val) + 
  coord_sf() + theme_few() + theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + 
  labs(x="", y="")

## -----------------------------------------------------------------------------
dat <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  mutate(var = factor(var, levels=paste0("bio", 1:19))) %>% 
  filter(var == "bio12") %>%
  group_by(x, y, time_frame, rcp, var) %>% summarise(value=mean(value))
col_val <- scales::rescale(quantile(dat$value, probs=seq(0,1,0.12)))
lim <- c(min(dat$value, na.rm=T), max(dat$value, na.rm=T))
dat %>% ggplot() + geom_tile(aes(x=x, y=y, fill=value)) + geom_sf(data=europe, fill="NA") + 
  facet_grid(rcp~time_frame, switch = "y") + scale_fill_gradientn(name="bio12", colours=whiteblue, limits=lim, values=col_val) + 
  coord_sf() + theme_few() + theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + 
  labs(x="", y="")

## -----------------------------------------------------------------------------
dat <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>%
  mutate(var = factor(var, levels=paste0("bio", 1:19))) %>% 
  filter(var == "bio15") %>%
  group_by(x, y, time_frame, rcp, var) %>% summarise(value=mean(value))
col_val <- scales::rescale(quantile(dat$value, probs=seq(0,1,0.12)))
lim <- c(min(dat$value, na.rm=T), max(dat$value, na.rm=T))
dat %>% ggplot() + geom_tile(aes(x=x, y=y, fill=value)) + geom_sf(data=europe, fill="NA") + 
  facet_grid(rcp~time_frame, switch = "y") + scale_fill_gradientn(name="bio15", colours=whiteblue, limits=lim, values=col_val) + 
  coord_sf() + theme_few() + theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + 
  labs(x="", y="")
rm(dat); invisible(gc())

## ---- fig.width=10, fig.height=7----------------------------------------------
# Plot bioclimatic changes between 1991-2020 and 2021 - 2050
bioclim_30yr <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>% 
  group_by(x, y, time_frame, var, rcp) %>% summarise(value=mean(value)) %>%
  tidyr::spread(time_frame,value) %>% filter(var == "bio1")
bioclim_30yr$change <- (bioclim_30yr$`2041-2070`-bioclim_30yr$`1991-2020`)
bioclim_30yr$rel_change <- (bioclim_30yr$`2041-2070`-bioclim_30yr$`1991-2020`)/bioclim_30yr$`1991-2020`*100

col_val <- scales::rescale(quantile(bioclim_30yr$change, probs=seq(0,1,0.12)))
lim <- c(min(bioclim_30yr$change, na.rm=T), max(bioclim_30yr$change, na.rm=T))
p1 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=change)) + facet_wrap(.~rcp, switch="y") + 
  scale_fill_gradientn(name="change", colours=whitered, limits=lim, values=col_val) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio1")

#col_val <- scales::rescale(quantile(bioclim_30yr$rel_change, probs=seq(0,1,0.12)))
#lim <- c(min(bioclim_30yr$rel_change, na.rm=T), max(bioclim_30yr$rel_change, na.rm=T))
p2 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=rel_change)) + facet_wrap(.~rcp, switch="y") + 
  scale_fill_gradientn(name="% change", colours=whitered) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio1")
p1 / p2
rm(p1, p2)

## ---- fig.width=10, fig.height=7----------------------------------------------
# Plot bioclimatic changes between 1991-2020 and 2021 - 2050
bioclim_30yr <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>% 
  group_by(x, y, time_frame, var, rcp) %>% summarise(value=mean(value)) %>%
  tidyr::spread(time_frame,value) %>% filter(var == "bio12")
bioclim_30yr$change <- (bioclim_30yr$`2041-2070`-bioclim_30yr$`1991-2020`)
bioclim_30yr$rel_change <- (bioclim_30yr$`2041-2070`-bioclim_30yr$`1991-2020`)/bioclim_30yr$`1991-2020`*100

col_val <- scales::rescale(quantile(bioclim_30yr$change, probs=seq(0,1,0.12)))
lim <- c(min(bioclim_30yr$change, na.rm=T), max(bioclim_30yr$change, na.rm=T))
p1 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=change)) + facet_wrap(.~rcp, switch="y") + 
  scale_fill_gradientn(name="change", colours=whiteblue, limits=lim, values=col_val) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio12")

col_val <- scales::rescale(quantile(bioclim_30yr$rel_change, probs=seq(0,1,0.12)))
lim <- c(min(bioclim_30yr$rel_change, na.rm=T), max(bioclim_30yr$rel_change, na.rm=T))
p2 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=rel_change)) + facet_wrap(.~rcp, switch="y") + 
  scale_fill_gradientn(name="% change", colours=whiteblue, limits=lim, values=col_val) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio12")
p1 / p2

## ---- fig.width=10, fig.height=7----------------------------------------------
# Plot bioclimatic changes between 1991-2020 and 2021 - 2050
bioclim_30yr <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>% 
  group_by(x, y, time_frame, var, rcp) %>% summarise(value=mean(value)) %>%
  tidyr::spread(time_frame,value) %>% filter(var == "bio1")
bioclim_30yr$change <- (bioclim_30yr$`2071-2100`-bioclim_30yr$`1991-2020`)
bioclim_30yr$rel_change <- (bioclim_30yr$`2071-2100`-bioclim_30yr$`1991-2020`)/bioclim_30yr$`1991-2020`*100

col_val <- scales::rescale(quantile(bioclim_30yr$change, probs=seq(0,1,0.12)))
lim <- c(min(bioclim_30yr$change, na.rm=T), max(bioclim_30yr$change, na.rm=T))
p1 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=change)) + facet_wrap(.~rcp) + 
  scale_fill_gradientn(name="change", colours=whitered, limits=lim, values=col_val) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio1")

#col_val <- scales::rescale(quantile(bioclim_30yr$rel_change, probs=seq(0,1,0.12)))
#lim <- c(min(bioclim_30yr$rel_change, na.rm=T), max(bioclim_30yr$rel_change, na.rm=T))
p2 <- bioclim_30yr %>% 
  ggplot() + geom_tile(aes(x=x, y=y, fill=rel_change)) + facet_wrap(.~rcp) + 
  scale_fill_gradientn(name="% change", colours=whitered) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio1")
p1 / p2

## ---- fig.width=10, fig.height=7----------------------------------------------
# Plot bioclimatic changes between 1991-2020 and 2021 - 2050
bioclim_30yr <- tidyr::gather(cordex_bioclim_eur_p44deg, var, value, -c(x,y,time_frame,gcm,rcp,rcm,ensemble,rs)) %>% 
  group_by(x, y, time_frame, var, rcp) %>% summarise(value=mean(value)) %>%
  tidyr::spread(time_frame,value) %>% filter(var == "bio12")
bioclim_30yr$change <- (bioclim_30yr$`2071-2100`-bioclim_30yr$`1991-2020`)
bioclim_30yr$rel_change <- (bioclim_30yr$`2071-2100`-bioclim_30yr$`1991-2020`)/bioclim_30yr$`1991-2020`*100

col_val <- scales::rescale(quantile(bioclim_30yr$change, probs=seq(0,1,0.12)))
lim <- c(min(bioclim_30yr$change, na.rm=T), max(bioclim_30yr$change, na.rm=T))
p1 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=change)) + facet_wrap(.~rcp) + 
  scale_fill_gradientn(name="change", colours=whiteblue, limits=lim, values=col_val) +
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio12")

col_val <- scales::rescale(quantile(bioclim_30yr$rel_change, probs=seq(0,1,0.12)))
lim <- c(min(bioclim_30yr$rel_change, na.rm=T), max(bioclim_30yr$rel_change, na.rm=T))
p2 <- bioclim_30yr %>% ggplot() + geom_tile(aes(x=x, y=y, fill=rel_change)) + facet_wrap(.~rcp) + 
  scale_fill_gradientn(name="% change", colours=whiteblue, limits=lim, values=col_val) + 
  geom_sf(data=europe, fill="NA") + coord_sf() + theme_few() + 
  theme(strip.placement = "outside", axis.text.x = element_text(angle=45, vjust=0.5)) + labs(x="", y="bio12")
p1 / p2
rm(list=ls()); gc()

