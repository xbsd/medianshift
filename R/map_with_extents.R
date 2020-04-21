#!/usr/bin/env Rscript
#
## ---------------------------
##
## Script name: map_with_extents.R
##
## Purpose of script: map extents of AOIs. Extents are in ..data/extents.csv
##
## Author: Leonidas Liakos
##
## Date Created: 08/04/2020
##
##
## ---------------------------
##
## Notes: -
##
##
## ---------------------------
library(raster)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(ggsflabel)
library(ggsn)


####====== Read settings ==================== ####

cfg <- config::get(file = here::here("R", "config.yml"))
####========================================= ####

csv <-  read.csv(file=here::here("data",cfg$EXTENTS), header=TRUE, sep=";") %>%
    mutate(ID = row_number())

extents <-  csv %>%  dplyr::select(xmin, xmax, ymin, ymax)

f<- function(x) {
    m <- matrix(x, nrow = 2, ncol = 2,byrow=TRUE)
    colnames(m) <- c("xmin","ymin")
    rownames(m) <- c("x","y")
    ext <- raster::extent(m)
    e <- as(ext, "SpatialPolygons")
    return(e)
}

poly_list <- apply(extents, 1, f)

polys <-do.call(bind, poly_list)

proj4string(polys) <- sp::CRS("+init=epsg:2100")

p <- SpatialPolygonsDataFrame(polys, extents, match.ID = TRUE)
p_wgs84 <- p %>% st_as_sf() %>% sf::st_transform( 4326)


p1 <- spTransform(p, crs(sp::CRS("+init=epsg:4326")))
ext <- extent(bbox(p1))
p <- p %>% st_as_sf()


p_centroid<- p %>%
    st_transform(4326) %>%
    sf:: st_centroid() %>%
    mutate(ID = row_number()) %>%
    dplyr::inner_join(csv, by=c("ID"))



world <- ne_countries(scale = "large", returnclass = "sf")
label <-data.frame(x=22,y=39,NAME="Greece")
ggplot() +
    theme_set(theme_bw())+
    #theme_bw(base_size = 7)+
    #ggtitle("Regions of interest") +
    geom_sf(data = world, fill = "grey91", col = 'grey32', size=0.1) +
    theme(
        panel.grid.major = element_line(
            color = gray(.5),
            linetype = "dotted",
            size = 0.15
        ),
        panel.background = element_rect(fill = "white")
    ) +

    geom_sf(data = p,
            size=0.2,
            colour = "black",
            fill = NA) +
    coord_sf(
        xlim = c(ext@xmin - 1.5, ext@xmax + 1.5),
        ylim = c(ext@ymin - 0.7, ext@ymax + 0.7),
        expand = FALSE

    )+

        # Add one annotation
  geom_sf_label(data = world,
                aes(label = name),
                size = 3.5,
                #nudge_x = 1,
                nudge_y=0.5)+
                

    geom_sf_text_repel(data = p_centroid,
                 aes( label = NAME),
                 colour = "black" ,
                 #fontface = "bold",
                 point.padding=0.05,
                 box.padding=0.6,
                 #nudge_x=1,
                 #nudge_y=0.5,
                 # label size
                 size = 3.3,# 3.5 
                 # color of the line segments
                 segment.colour = "black",
                 # Draw an arrow from the label to the data point.
                 arrow = arrow(length = unit(0.01, 'npc')),
                 # line segment transparency
                 segment.alpha = 1,
                 # line segment thickness
                 segment.size = .3,
                 seed = 5,
                 # the force of repulsion between overlapping text labels
                 force = 50,
                 # maximum number of iterations to attempt to resolve overlaps
                 max.iter = 10e3,
                 # turn labels off in the legend
                 show.legend = FALSE)+
     theme(axis.title = element_blank())+
    theme(axis.text=element_text(size=8.0))+
    theme(axis.ticks.length=unit(.04, "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.margin=margin(t = 0, r = 0.5, b = 0, l = 0.2, "cm"),
          panel.margin=margin(t = 0, r = 0, b = 0, l = 0, "cm"))
  

# save as png
ggsave(
    'fig.2_map_extents.tif',
    plot = last_plot(),
    device = "tiff",
    path = here::here("output"),
    scale = 1,
    width = 9,
    height =8,
    units = c("cm"),
    dpi = 600,
    limitsize = TRUE
)
  
