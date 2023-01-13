#' R Package to access European species, climate and land-use data
#' 
#' @name edc package
#' @aliases edcpackage
#' @docType package
#' @title R Package to access European species, climate and land-use data
#' @author RS-eco
#'
#' @importFrom dplyr left_join %>%
#' @references https://www.iucnredlist.org/resources/spatial-data-download; 
#' http://datazone.birdlife.org/species/spcdistPOS; 
#' https://datadryad.org/resource/doi:10.5061/dryad.83s7k
#' @keywords package
#'
NULL
#'
#' @docType data
#' @name biogeographic_regions_eur
#' @title Biogeographic regions of Europe
#' @description Biogeographical regions of Europe
#' @usage data(biogeographic_regions_eur)
#' @details A gird-based map for the biogeographical regions of Europe
#' @format A \code{data.frame} with 251138 observations and 4 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://zenodo.org/record/3766175#.XqbP5GgzZPZ}}
#' @references Cervellini M, Zannini P, Di Musciano M, Fattorini S, Jimenez-Alfaro B, Rocchini D, Field R, R. Vetaas O,
#' Irl SD.H, Beierkuhnlein C, Hoffmann S, Fischer J-C, Casella L, Angelini P, Genovesi P, Nascimbene J, Chiarucci A (2020) 
#' A grid-based map for the Biogeographical Regions of Europe. Biodiversity Data Journal 8: e53720. 
#' \url{https://doi.org/10.3897/BDJ.8.e53720}
NULL
#'
#' @docType data
#' @name birdlife_eur_p11deg
#' @title Distribution of terrestrial birds at 0.11 degree spatial resolution for Europe
#' @description Data with x,y coordinates for each terrestrial bird species (present in 0.11 degree grid cells)
#' @details This dataset contains the distribution of all terrestrial bird species (on a 0.11 degree grid) for Europe
#' according to the IUCN range maps.
#' @format A \code{data.frame} with 10728270 observations and 7 variables.
NULL
#'
#' @docType data
#' @name cci_eur
#' @title Fraction of ESA CCI land cover time series (1992 - 2018) of Europe
#' @description Fraction of ESA CCI land cover time series (1992 - 2018) of Europe
#' @usage data(cci_eur)
#' @details Fraction of ESA CCI land cover time series (1992 - 2018) of Europe.
#' @format A \code{data.frame} with 249634 observations and 30 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://zenodo.org/record/3730469}}
NULL
#'
#' @docType data
#' @name cordex_bioclim_eur_p44deg
#' @title Euro-cordex bioclimatic data of Europe
#' @description Euro-cordex bioclimatic data of Europe gridded to 0.44 degree
#' @usage data(cordex_bioclim_eur_p44deg)
#' @details Euro-cordex climate simulation data of Europe.
#' The data.frame contains bioclimatic data for 3 time periods (1971-2000, 2021-2050, 2071-2100) 
#' under 3 representative concentration pathways (rcps; RCP2.6, RCP4.5, RCP8.5), 
#' 5 global circulation models (gcms; MPI-M-MPI-ESM-LR, CNRM-CERFACS-CNRM-CM5, 
#' MOHC-HadGEM2-ES, ICHEC-EC-EARTH, IPSL-IPSL-CM5A-MR), 5 regional climate models (rcms; CLMcom-CCLM4-8-17, 
#' KNMI-RACMO22E, MPI-CSC-REMO2009, SMHI-RCA4, DMI-HIRHAM5), 
#' 4 ensemble (r1i1p1, r2i1p1, r3i1p1, r12i1p1) and 3 rs (v1, v1a, v2).
#' @format A \code{data.frame} with 82782 observations and 27 variables.
NULL
#'
#' @docType data
#' @name corine_lc_eur_p11deg
#' @title Corine land cover data of Europe
#' @description Corine land cover data of Europe
#' @usage data(corine_lc_eur_p11deg)
#' @details Corine land cover data of Europe resampled onto a 0.11 degree grid.
#' Corine Land cover are land cover maps for Europe mapped at a resolution of 100 x 100 m.
#' @format A \code{data.frame} with 834380 observations and 8 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://land.copernicus.eu/pan-european/corine-land-cover}}
NULL
#'
#' @docType data
#' @name corine_lc_eur_p44deg
#' @title Corine land cover data of Europe
#' @description Corine land cover data of Europe
#' @usage data(corine_lc_eur_p44deg)
#' @details Corine land cover data of Europe resampled onto a 0.44 degree grid.
#' Corine Land cover are land cover maps for Europe mapped at a resolution of 100 x 100 m.
#' @format A \code{data.frame} with 87563 observations and 8 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://land.copernicus.eu/pan-european/corine-land-cover}}
NULL
#'
#' @docType data
#' @name ebba2_shp
#' @title shapefile of EBBA2 grid
#' @description outline of EBBA2 grid
#' @usage data(ebba2_shp)
#' @details Boundaries of EBBA2 grid.
#' @format A \code{sfc_MULTIPPOLYGON} object.
NULL
#' 
#' @docType data
#' @name ebcc_eur
#' @title European Bird Census Council data of Europe
#' @description EBCC bird data of Europe
#' @usage data(ebcc_eur)
#' @details European Breeding Bird Atlas data
#' @format A \code{data.frame} with 3920 observations and 497 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://www.gbif.org}}. 
#' More information can be found here: \itemize{\item \url{https://www.ebcc.info/what-we-do/ebba2/ebcc-atlas-of-european-breeding-birds/}}.
NULL
#'
#' @docType data
#' @name ebcc_Phylloscopus_bonelli
#' @title European Bird Census Council data of Phylloscopus bonelli
#' @description EBCC bird data of Phylloscopus bonelli
#' @usage data(ebcc_Phylloscopus_bonelli)
#' @details European Breeding Bird Atlas data of  Phylloscopus bonelli
#' @format A \code{data.frame} with 2622 observations and 5 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://www.gbif.org}}. 
#' More information can be found here: \itemize{\item \url{https://www.ebcc.info/what-we-do/ebba2/ebcc-atlas-of-european-breeding-birds/}}.
NULL
#'
#' @docType data
#' @name ebcc_Sylvia_curruca
#' @title European Bird Census Council data of Sylvia curruca
#' @description EBCC bird data of Sylvia curruca
#' @usage data(ebcc_Sylvia_curruca)
#' @details European Breeding Bird Atlas data of Sylvia curruca
#' @format A \code{data.frame} with 2862 observations and 5 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://www.gbif.org}}. 
#' More information can be found here: \itemize{\item \url{https://www.ebcc.info/what-we-do/ebba2/ebcc-atlas-of-european-breeding-birds/}}.
NULL
#'
#' @docType data
#' @name ebcc_Turdus_torquatus
#' @title European Bird Census Council data of Turdus torquatus
#' @description EBCC bird data of Turdus torquatus
#' @usage data(ebcc_Turdus_torquatus)
#' @details European Breeding Bird Atlas data of Turdus torquatus
#' @format A \code{data.frame} with 2642 observations and 5 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://www.gbif.org}}. 
#' More information can be found here: \itemize{\item \url{https://www.ebcc.info/what-we-do/ebba2/ebcc-atlas-of-european-breeding-birds/}}.
NULL
#'
#' @docType data
#' @name ebcc_Vanellus_vanellus
#' @title European Bird Census Council data of Vanellus vanellus
#' @description EBCC bird data of Vanellus vanellus
#' @usage data(ebcc_Vanellus_vanellus)
#' @details European Breeding Bird Atlas data of Vanellus vanellus
#' @format A \code{data.frame} with 2856 observations and 5 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://www.gbif.org}}. 
#' More information can be found here: \itemize{\item \url{https://www.ebcc.info/what-we-do/ebba2/ebcc-atlas-of-european-breeding-birds/}}.
NULL
#'
#' @docType data
#' @name e-obs_bioclim_1967_1997_eur
#' @title E-OBS gridded bioclimatic data for Europe at 0.1 degree resolution from 1967 to 1997
#' @description E-OBS gridded bioclimatic data for Europe at 0.1 degree resolution from 1967 to 1997
#' @usage data(`e-obs_bioclim_1967_1997_eur`)
#' @details E-OBS gridded elevation data for Europe at 0.1 degree resolution from 1967 to 1997.
#' @format A \code{data.frame} with 43008 observations and 21 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=form}}
#' @references Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: 
#' An Ensemble Version of the E-OBS Temperature and Precipitation Datasets, 
#' J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200
NULL
#'
#' @docType data
#' @name e-obs_bioclim_1987_2017_eur
#' @title E-OBS gridded bioclimatic data for Europe at 0.1 degree resolution from 1987 to 2017
#' @description E-OBS gridded bioclimatic data for Europe at 0.1 degree resolution from 1987 to 2017
#' @usage data(`e-obs_bioclim_1987_2017_eur`)
#' @details E-OBS gridded bioclimatic data for Europe at 0.1 degree resolution from 1987 to 2017.
#' @format A \code{data.frame} with 43008 observations and 21 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=form}}
#' @references Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: 
#' An Ensemble Version of the E-OBS Temperature and Precipitation Datasets, 
#' J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200
NULL
#'
#' @docType data
#' @name e-obs_elev_ens_eur_p1deg
#' @title E-OBS gridded elevation data for Europe at 0.1 degree resolution
#' @description E-OBS gridded elevation data for Europe at 0.1 degree resolution
#' @usage data(`e-obs_elev_ens_eur_p1deg`)
#' @details E-OBS gridded elevation data for Europe at 0.1 degree resolution.
#' @format A \code{data.frame} with 44200 observations and 3 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=form}}
#' @references Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: 
#' An Ensemble Version of the E-OBS Temperature and Precipitation Datasets, 
#' J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200
NULL
#'
#' @docType data
#' @name eu_dem_eur_p11deg
#' @title Eurostat data of elevation across Europe
#' @description European elevation data at a resolution of 25 m resampled to 0.11 degree
#' @usage data(eu_dem_eur_p11deg)
#' @details European elevation data derived from Eurostat at a resolution of 25 m resampled to 0.11 degree.
#' @format A \code{data.frame} with 63258 observations and 7 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://ec.europa.eu/eurostat/de/web/gisco/geodata/reference-data/elevation/eu-dem/eu-dem-laea}}.
NULL
#'
#' @docType data
#' @name eu_dem_eur_p44deg
#' @title Eurostat data of elevation across Europe
#' @description European elevation data at a resolution of 25 m resampled to 0.44 degree
#' @usage data(eu_dem_eur_p44deg)
#' @details European elevation data derived from Eurostat at a resolution of 25 m resampled to 0.44 degree.
#' @format A \code{data.frame} with 4868 observations and 7 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://ec.europa.eu/eurostat/de/web/gisco/geodata/reference-data/elevation/eu-dem/eu-dem-laea}}.
NULL
#'
#' @docType data
#' @name eur_27
#' @title administrative boundary of EU27 countries
#' @description GADM outline of EU27 countries
#' @usage data(eur_27)
#' @details Administrative boundary of EU27 countries derived from GADM data.
#' @format A \code{sfc_MULTIPPOLYGON} object.
#' @source This data has been obtained from: \itemize{\item \url{https://www.gadm.org/data.html}}
NULL
#'
#' @docType data
#' @name eur_countries
#' @title administrative boundary of EU27 countries plus Switzerland and UK
#' @description GADM outline of EU27 countries plus Switzerland and UK
#' @usage data(eur_countries)
#' @details Administrative boundary of EU27 countries plus Switzerland and UK 
#' derived from GADM data.
#' @format A \code{sfc_MULTIPPOLYGON} object.
#' @source This data has been obtained from: \itemize{\item \url{https://www.gadm.org/data.html}}
NULL
#'
#' @docType data
#' @name europe
#' @title administrative boundary of Europe
#' @description GADM outline of Europe
#' @usage data(europe)
#' @details Administrative boundary of Europe derived from GADM data.
#' @format A \code{sfc_MULTIPPOLYGON} object.
#' @source This data has been obtained from: \itemize{\item \url{https://www.gadm.org/data.html}}
NULL
#'
#' @docType data
#' @name hochkirch_orthoptera_eur
#' @title Distribution of Orthoptera species for Europe
#' @description Data with range maps for each Orthoptera species in Europe
#' @details This dataset contains the distribution of all Orthoptera species in Europe
#' according to the range polygons by Hochkirch et al.
#' @format A \code{data.frame} with 1355 observations and 5 variables.
NULL
#'
#' @docType data
#' @name iucn_odonata_eur_p11deg
#' @title Distribution of Odonata species at 0.11 degree spatial resolution for Europe
#' @description Data with x,y coordinates for each Odonata species (present in 0.11 degree grid cells) in Europe
#' @details This dataset contains the distribution of all Odonata species (on a 0.11 degree grid) in Europe
#' according to the IUCN ranges.
#' @format A \code{data.frame} with 662221 observations and 7 variables.
NULL
#'
#' @docType data
#' @name kalkman_odonata_eur
#' @title Distribution of Odonata species in Europe
#' @description Data with x,y coordinates for each Odonata species in Europe
#' @details This dataset contains the distribution of Odonata species in Europe according to Kalkman et al. 
#' @format A \code{data.frame} with 6552 observations and 156 variables.
NULL
#'
#' @docType data
#' @name landsystem_eur_1km
#' @title Landsystem data of Europe
#' @description Landsystem data of Europe
#' @usage data(landsystem_eur_1km)
#' @details Landsystem data of Europe
#' @format A \code{data.frame} with 4894729 observations and 3 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://dataverse.nl/dataset.xhtml?persistentId=doi:10.34894/XNC5KA}}
#' @references Dou, Y., Cosentino, F., Malek, Z. et al. (2021) 
#' A new European land systems representation accounting for landscape characteristics. Landscape Ecology. 
#' \url{https://doi.org/10.1007/s10980-021-01227-5}.
NULL
#'
#' @docType data
#' @name landsystem_perc_eur_10km
#' @title Percentage cover landsystem data of Europe gridded to 10x10 km
#' @description Percentage cover landsystem data of Europe gridded to 10x10 km
#' @usage data(landsystem_perc_eur_10km)
#' @details Percentage cover landsystem data of Europe gridded to 10x10 km
#' @format A \code{data.frame} with 57335 observations and 28 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://dataverse.nl/dataset.xhtml?persistentId=doi:10.34894/XNC5KA}}
#' @references Dou, Y., Cosentino, F., Malek, Z. et al. (2021) 
#' A new European land systems representation accounting for landscape characteristics. Landscape Ecology. 
#' \url{https://doi.org/10.1007/s10980-021-01227-5}.
NULL
#'
#' @docType data
#' @name srtm_csi_eur_p11deg
#' @title srtm data of Europe
#' @description NASA SRTM v3.0 elevation data of Europe
#' @usage data(srtm_csi_eur_p11deg)
#' @details Elevation data of Europe derived from NASA SRTM 3 Arc Sec v3.0 data and aggregated to 0.11 degree resolution.
#' @format A \code{data.frame} with 46181 observations and 7 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://srtm.csi.cgiar.org/srtmdata/}}
NULL
#'
#' @docType data
#' @name srtm_csi_eur_p44deg
#' @title srtm data of Europe
#' @description NASA SRTM v3.0 elevation data of Europe
#' @usage data(srtm_csi_eur_p44deg)
#' @details Elevation data of Europe derived from NASA SRTM 3 Arc Sec v3.0 data and aggregated to 0.44 degree resolution.
#' @format A \code{data.frame} with 2996 observations and 7 variables.
#' @source This data has been obtained from: \itemize{\item \url{http://srtm.csi.cgiar.org/srtmdata/}}
NULL
#'
#' @docType data
#' @name wfde5_rainf_eur
#' @title WFDE5 reanalysis yearmon rainfall flux data of Europe
#' @description WFDE5 reanalysis yearmon rainfall flux data of Europe
#' @usage data(wfde5_rainf_eur)
#' @details WFDE5 reanalysis yearmon rainfall flux data of Europe.
#' The data.frame contains rainfall flux in kg m-2 s-1 from January 1979 - December 2019 at a 0.5° grid
#' @format A \code{data.frame} with 3507 observations and 494 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.20d54e34?tab=form}}
NULL
#'
#' @docType data
#' @name wfde5_snowf_eur
#' @title WFDE5 reanalysis yearmon snowfall flux data of Europe
#' @description WFDE5 reanalysis yearmon snowfall flux data of Europe
#' @usage data(wfde5_rainf_eur)
#' @details WFDE5 reanalysis yearmon snowfall flux data of Europe.
#' The data.frame contains snowfall flux in kg m-2 s-1 from January 1979 - December 2019 at a 0.5° grid
#' @format A \code{data.frame} with 3507 observations and 494 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.20d54e34?tab=form}}
NULL
#'
#' @docType data
#' @name wfde5_tair_eur
#' @title WFDE5 reanalysis yearmon air temperature data of Europe
#' @description WFDE5 reanalysis yearmon air temperature data of Europe
#' @usage data(wfde5_rainf_eur)
#' @details WFDE5 reanalysis yearmon near-surface air temperature in K, 
#' the temperature of air at 2 metres above the surface of land, sea or inland waters, data of Europe.
#' The data.frame contains air temperature data from January 1979 - December 2019 at a 0.5° grid
#' @format A \code{data.frame} with 3507 observations and 494 variables.
#' @source This data has been obtained from: \itemize{\item \url{https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.20d54e34?tab=form}}
NULL