#!/usr/bin/env Rscript
#
## ---------------------------
##
## Script name: miller_plor.R
##
## Purpose of script: to generate Lunar irradiance plot, based on mt2009 script
##
## Author: Leonidas Liakos
##
## Date Created: 05/04/2020
##
##
## ---------------------------
##
## Run order: -
##
##
## ---------------------------
#####  

library(reticulate)
library(lunar)
library(ggplot2)

####====== Read settings ==================== ####
Sys.setenv(R_CONFIG_ACTIVE = name)
cfg <- config::get(file = here::here("R", "config.yml"))
####========================================= ####

source(here::here("R", 'myfunctions.R'))


use_condaenv(condaenv = 'earth', required = TRUE) # Name of Conda environment to activate, conda_list() to list available envs
os <- import("os")
os$chdir(here::here("Python"))
source_python("lunar_irrad_DNB.py")

#mt2009("202001010101") #Lunar Irradiance (mW/m^2-micron)
d <- seq(as.Date(sprintf("%s/1/1", cfg$YEAR)), as.Date(sprintf("%s/12/31", cfg$YEAR)), "days")
d_str <- paste(format(d, "%Y%m%d")  , "0101", sep = "")
miller <- sapply(d_str, function(x)
    mt2009(x)$lun_irrad_scl)
my_df <- data.frame(mydates = d, lunar_ir = miller)



m <- moon_phases(d)
dateRanges =data.frame(from=d[m$fullmoon_starts], to=d[m$fullmoon_ends])

(g <- ggplot2::ggplot(data = my_df) +
    scale_fill_manual('',
                      values = 'grey',
                      guide = guide_legend(override.aes = list(alpha = 1))) +
        geom_rect(
        data = dateRanges,
        aes(
            xmin = from,
            xmax = to,
            ymin = -Inf,
            ymax = Inf,
            fill = "Full moon"
        ),
        alpha = 0.5
    ) +
    ggtitle(sprintf("Lunar Irradiance (mW/m^2-micron), 2018, at 01:00")) +
    geom_line( data = my_df,aes(x = as.Date(d , "X%Y.%m.%d"), y =lunar_ir), size=0.3) +
    xlab("Date") +
    ylab("Lunar Irradiance (mW/m^2-micron)") +
    scale_y_continuous(limit = c(0, max(my_df$lunar_ir))) +theme_gray(base_size = 10))



# save as png
ggsave(
    sprintf('fig.1_lunar_irradiance_mt2009_%s.tif', "2018"),
    plot = g,
    device = "tiff",
    path = here::here('output'),
    scale = 1,
    width = 9,
    height = 5,
    units = c("cm"),
    dpi = 600,
    limitsize = TRUE
)
