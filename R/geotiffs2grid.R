#!/usr/bin/env Rscript
#
## ---------------------------
##
## Script name: geotiffs2grid.R
##
## Purpose of script: save geotiffs as grd files.
##
## Author: Leonidas Liakos
##
## Date Created: 11/02/2020
##
##
## ---------------------------
##
## Run order: 1
##
##
## ---------------------------
library(magrittr)
library(rasterVis)
library(mapview)
library(mapedit)
library(purrr)
library(ggplot2)
library(lunar)
library(here)

source(here::here("R", 'myfunctions.R'))

####====== Read settings ==================== ####
cfg <- config::get(file = here::here("R", "config.yml"))
####========================================= ####


# run the script for each directory of geotiffs

#DNB_geotiffs_DIR <- "BlackMarble_VNP46A1_BIG_BRDF_vnp43ma4v001_h19v05"
DNB_geotiffs_DIR <- "BlackMarble_VNP46A1_BIG_BRDF_vnp43ma4v001_h20v05"



# Set directories variables -----------------------------------------------
geotiffs_DIR <- here::here("data","geotiffs", DNB_geotiffs_DIR)



# Dates with dummy data -----------------------------------------------
dummy_files <- file.path(geotiffs_DIR, c("VNP46A1.A2017204.h20v05.001.2019184180811.DNB_201707230108_#2017-204#_2100.tif",
                                                 "VNP46A1.A2012092.h20v05.001.2019088062532.DNB_201204010024_#2012-092#_2100.tif"
                                                 ))

dnb_files <- list.files(geotiffs_DIR , pattern = sprintf("^VNP46A1.A%s.*\\.DNB*.*tif$",cfg$YEAR), full.names = T) #lunar corrected files (roman method)
original_dnb_files <- list.files(geotiffs_DIR, pattern = sprintf("^VNP46A1.A%s.*original*.*tif$",cfg$YEAR), full.names = T) #original dnb files.no lunar correction


# remove dummy files 
indexes_to_keep <- which(!dnb_files %in%  dummy_files)
dnb_files <- dnb_files[indexes_to_keep]
original_dnb_files <- original_dnb_files[indexes_to_keep]


# Generare daily composites -----------------------------------------------
# 
# The original data:
original_daily_DNB <- generate_daily_dnb(original_dnb_files) # original data

### Write in disk
writeRaster(
    original_daily_DNB,
    here::here("data","grd", DNB_geotiffs_DIR, cfg$original_DNB_grd),
    format = "raster",
    overwrite = TRUE,
    datatype = 'INT2S' # χάνουμε λίγη πληροφορία με τους ακέραιους
) #read it with stack()



#### The lunar corrected data:
daily_DNB <- generate_daily_dnb(dnb_files) # lunar corrected (roman method) data

writeRaster(
    daily_DNB,
    here::here("data","grd", DNB_geotiffs_DIR, cfg$DNB_grd),
    format = "raster",
    overwrite = TRUE,
    datatype = 'INT2S' # χάνουμε λίγη πληροφορία με τους ακέραιους
) #read it with stack()

