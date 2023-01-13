edc: European Data Cube
================

## Installation

You can install edc from github with:

``` r
# Install remotes if not available
if (!"remotes" %in% installed.packages()[, "Package"]) install.packages("remotes")

# Install edc package from Github
remotes::install_github("RS-eco/edc", build_vignettes = T)
```

If package installation with `install_github()` fails, try this as an
alternative:

``` r
tmp_zip <- tempfile(fileext = ".zip")
source_url <- "https://api.github.com/repos/RS-eco/edc/zipball/main"
utils::download.file(source_url, destfile = tmp_zip, method = "wget")
file.exists(tmp_zip)

remotes::install_local(tmp_zip)
```

After installation, simply load the edc package:

``` r
library(edc)
```

**If you encounter a bug or if you have any problems, please file an
[issue](https://github.com/RS-eco/edc/issues) on Github.**

## Datasets

### GADM data

Shapefile of Europe can be accessed by:

``` r
data("europe")
```

### Euro-Cordex Data

The bioclimatic Euro-Cordex data can be accessed by:

``` r
# Euro-Cordex data for Europe
load(system.file("extdata", "cordex_bioclim_eur.rda", package = "edc"))
# Data at 0.11 degree resolution

# There is also a version in 0.44 degree resolution
data("cordex_bioclim_eur_p44deg")

# Euro-Cordex data for Bavaria can be found in the bavDC package under:
# https://github.com/RS-eco/bavDC
```

**Note:** Please also have a look at the corresponding
[vignette](https://github.com/RS-eco/edc/blob/main/vignettes/edc-eurocordex.Rmd)
for detailed information about the dataset and the related variables.

### E-OBS climate data

E-OBS in situ gridded climate data of Europe at 1 degree resolution:

``` r
# Elevation
data("e-obs_elev_ens_eur_p1deg")

# Rainfall
load(system.file("extdata", "e-obs_rr_ens_mean_0.1deg_eur_yearmon.rda", package = "edc"))

# Mean temperature
load(system.file("extdata", "e-obs_tg_ens_mean_0.1deg_eur_yearmon.rda", package = "edc"))

# Minimum temperature
load(system.file("extdata", "e-obs_tn_ens_mean_0.1deg_eur_yearmon.rda", package = "edc"))

# Maximum temperature
load(system.file("extdata", "e-obs_tx_ens_mean_0.1deg_eur_yearmon.rda", package = "edc"))
```

### WFDE5 data

WFDE5 re-analysis climate (rainflux, snowflux, air temperature) data:

``` r
# Rainflux
data("wfde5_rainf_eur")

# Snowflux
data("wfde5_snowf_eur")

# Air temperature
data("wfde5_tair_eur")
```

### CCI land-cover data

Fraction of ESA CCI land cover time series (1992 - 2018) for Europe:

``` r
data("cci_eur")
```

### Corine Land-cover data

``` r
# Corine data for whole of Europe aggregated for three different resolutions
load(system.file("extdata", "corine_lc_eur_1km.rda", package = "edc"))
data("corine_lc_eur_p11deg")
data("corine_lc_eur_p44deg")
```

**Note:** Please also have a look at the corresponding
[vignette](https://github.com/RS-eco/edc/blob/main/vignettes/edc-landcover.Rmd)
for detailed information about the dataset and the related variables.

### EU-DEM data

``` r
# EU-DEM Elevation data aggregated to 2 km
load(system.file("extdata", "eu_dem_eur_2km.rda", package = "edc"))

# EU-DEM hillshade data aggregated to 2km
load(system.file("extdata", "hillshade_eur_2km.rda", package = "edc"))

# EU-DEM Elevation data aggregated to 0.11 & 0.44 deg
data("eu_dem_eur_p11deg")
data("eu_dem_eur_p44deg")
```

### SRTM data

The SRTM elevation data can be accessed by:

``` r
# SRTM data for Europe
load(system.file("extdata", "srtm_csi_eur_900m.rda", package = "edc"))
data("srtm_csi_eur_p11deg")
data("srtm_csi_eur_p44deg")
```

**Note:** Please also have a look at the corresponding
[vignette](https://github.com/RS-eco/edc/blob/main/vignettes/edc-elevation.Rmd)
for detailed information about the dataset and the related variables.

### IUCN data

The IUCN range map data can be accessed by:

``` r
# Odonata IUCN ranges for Europe
data("iucn_odonata_eur_p11deg")

# BirdLife ranges for Europe
data("birdlife_eur_p11deg")
```

**Note:** Please also have a look at the corresponding
[vignette](https://github.com/RS-eco/edc/blob/main/vignettes/edc-iucn.Rmd)
for detailed information about the dataset and the related variables.

### EBCC data

The EBCC Atlas data can be accessed by:

``` r
# Bird atlas data for Europe
data("ebcc_eur")

# EBCC data for selected bird species
data("ebcc_Sylvia_curruca")
data("ebcc_Turdus_torquatus")
data("ebcc_Phylloscopus_bonelli")
data("ebcc_Vanellus_vanellus")
```

**Note:** Please also have a look at the corresponding
[vignette](https://github.com/RS-eco/edc/blob/main/vignettes/edc-ebcc.Rmd)
for detailed information about the dataset and the related variables.

### Kalkman data

The Kalkman Odonata Atlas data for Europe can be accessed by:

``` r
data("kalkman_odonata_eur")
```

### Hochkirch data

The Hochkirch Orthoptera Atlas data for Europe can be accessed by:

``` r
data("hochkirch_orthoptera_eur")
```

**Note:** Please also have a look at the corresponding
[vignette](https://github.com/RS-eco/edc/blob/main/vignettes/edc-atlas.Rmd)
for detailed information about the datasets and the related variables.

### Landsystem data

The landsystem data for Europe can be accessed by:

``` r
# Original 1km data
data("landsystem_eur_1km")

# Percentage cover of each landsystem class on a 10x10 km
data("landsystem_perc_eur_10km")
```

### Biogeographic regions of Europe

A gridded dataset of the biogeographic regions of Europe can be accessed
by:

``` r
data("biogeographic_regions_eur")
```

### Zoning of UNSECO Biosphere Reserves

The zoning of UNESCO Biosphere Reserves data for Europe can be accessed
by:

``` r
load(system.file("extdata", "unesco_br_zones_eur.rda", package = "edc"))  # simplified to reduce file size
```
