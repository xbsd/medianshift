library(cluster)
library(dplyr)
library(rasterVis)
library(purrr)
library(ggplot2)
library(lunar)
library(raster)
library(lubridate)
library(sf)
library(parallel)
library(furrr)
library(suncalc)



source(here::here("R", 'myfunctions.R'))


####====== Read settings ==================== ####
name <- "Mykonos"
Sys.setenv(R_CONFIG_ACTIVE = name)
cfg <- config::get(file = here::here("R", "config.yml"))
####========================================= ####

# set some settings
options(scipen=999)


# generate extent object for csv data for AOI -----------------------------
# 
extents <-  read.csv(file=here::here("data",cfg$EXTENTS), header=TRUE, sep=";") %>%
    dplyr::filter(NAME==name) %>%
    dplyr::select(xmin, xmax, ymin, ymax) %>% unlist()

# name not in CSV file then stop
try(if(length(extents)==0) stop("AOI name is not available in csv file"))

m <- matrix(extents, nrow = 2, ncol = 2,byrow=TRUE)
colnames(m) <- c("xmin","ymin")
rownames(m) <- c("x","y")
ext <- raster::extent(m)





# Read/crop and save raster stacks for AOI --------------------------------

# Read
DNB_na.original  <- raster::stack(here::here("data","grd", cfg$DNB_geotiffs_DIR, cfg$original_DNB_grd))
DNB_na <- raster::stack(here::here("data","grd", cfg$DNB_geotiffs_DIR, cfg$DNB_grd))

# Crop
DNB_na.original  <- crop(DNB_na.original, ext) #original data
DNB_na  <- crop(DNB_na, ext) # lunar corrected (roman method)
mydates <- as.Date(names( raster::stack(DNB_na.original) ), "X%Y.%m.%d")


# φίλτρο για ένα έτος
filter_dates <-
  mydates >= sprintf("%s-01-01", cfg$YEAR) &
  mydates <= sprintf("%s-12-31", cfg$YEAR)

DNB_na.original   <- subset(DNB_na.original  , which(filter_dates)) #keep only 2018
DNB_na  <- subset(DNB_na  , which(filter_dates)) #keep only 2018
mydates <- mydates[filter_dates]


#  κράτα μόνο τα layers που έχουν πλήθος > από non NA ---------------------
# επειδη το cellStats σκάει με μεγάλα stacks φτίαχνω δική μου function με purrr
countPixelsNonNa <- purrr::map(as.list(DNB_na.original),
                               ~ {
                                 raster::cellStats(.x, function(x, ...)
                                   sum(!is.na(x))) %>% setNames(names(.x))
                               }) %>%
  unlist()


my_cellStats <- function(s,f,na.rm){
     purrr::map(as.list(s),
         ~ {raster::cellStats(.x, function(x, ...) do.call(f, list(x, na.rm=na.rm)))}) %>%
    unlist() %>%
    setNames(names(s))

}


nonNA_theshold_count  <- 30 # όριο πάνω από το οπόιο κρατάμε layers
DNB_na.original       <-  DNB_na.original[[which(countPixelsNonNa >= nonNA_theshold_count)]]
DNB_na                <-  DNB_na[[which(countPixelsNonNa >= nonNA_theshold_count)]]


### extract dates from layer names of raster stack
mydates <- as.Date(names(raster::stack(DNB_na.original)), "X%Y.%m.%d")
# .............................................................................




# Lunar Data
moon_df <- data.frame(mydates = mydates, month = lubridate::month(mydates), year=lubridate::year(mydates), phase=lunar.illumination(mydates))

newmoons <- moon_df %>%
    dplyr::group_by(year, month) %>%
    dplyr::mutate(phase = min(phase)) %>%
    dplyr::left_join( moon_df, by=(c("month"="month", "year"="year", "phase"="phase"))) %>%
    dplyr::distinct(mydates.y) %>%
    dplyr::pull(mydates.y)

fullmoons <- moon_df %>%
    dplyr::group_by(year, month) %>%
    dplyr::mutate(phase = max(phase)) %>%
    dplyr::left_join( moon_df, by=(c("month"="month", "year"="year", "phase"="phase"))) %>%
    dplyr::distinct(mydates.y) %>%
    dplyr::pull(mydates.y)


# διορθωμένα raster stacks
DNB_na.original_corrected <- correct_DNB_30Days(DNB_na.original, newmoons)

# υπολογισμός mean διορθωμένων raster raster stacks
DNB_na.original_stat <- cellStats(DNB_na.original, mean, na.rm = T)
DNB_na.original_corrected_stat <- cellStats(DNB_na.original_corrected, mean, na.rm = T)
DNB_na_stat <- cellStats(DNB_na, mean, na.rm = T)


# just a plot with all time series
plot(DNB_na.original_stat, type="l")
lines(DNB_na_stat, col="green" )
lines(DNB_na.original_corrected_stat, col="blue" )







#  Metric 1. mean difference of the average value of a scene in on --------


mean_diff_original  <- mean_diff(DNB_na.original_stat)
mean_diff_corrected <-mean_diff(DNB_na.original_corrected_stat)
mean_diff_corrected_brdf <-mean_diff(DNB_na_stat)

(v_mean_diff <-
  c(mean_diff_original,
    mean_diff_corrected_brdf,
    mean_diff_corrected
    ) %>% setNames(c(
      "Mean diff Original",
      "Mean diff BRDF Corrected",
      "Mean diff Median Shift Corrected"

    )))





# Metric 2. The deviation to the second order regression line, fit --------

# Δοκιμάζουμε fit στα αδιόρθωτα και διορθωμένα δεδομένα, second order regression.
# Στόχος είναι στα διορθωμένα δεδομένα να έχουμε καλυτερο R2 και μικρότερα κατάλοιπα(?)
# Ελέγχουμε για R2 (Multiple R-squared) και για Residuals.
############# second order regression #############
reg_df <- data.frame(dnb=DNB_na.original_stat, dnb_corrected=DNB_na.original_corrected_stat,dnb_corrected_roman=DNB_na_stat, time=c(1:length(mydates)))
fit_original <- lm(formula = dnb ~ poly(time, 2), reg_df)
fit_corrected <- lm(formula = dnb_corrected ~ poly(time, 2), reg_df)
fit_corrected_roman <- lm(formula = dnb_corrected_roman ~ poly(time, 2), reg_df)

par(mfrow=c(1,2))

plot(x = reg_df$time, y=reg_df$dnb, xlab = "Days", ylab = "Original DNB")
predictedcounts <- predict(fit_original,list(time=reg_df$time), time2=reg_df$time^2)
lines(reg_df$time, predictedcounts,  lwd = 3, xlab = "Time (s)", ylab = "Counts",col = "blue")


plot(x = reg_df$time, y=reg_df$dnb_corrected, xlab = "Days", ylab = "Corrected DNB")
predictedcounts <- predict(fit_corrected,list(time=reg_df$time), time2=reg_df$time^2)
lines(reg_df$time, predictedcounts,  lwd = 3, xlab = "Time (s)", ylab = "Counts",col = "blue")

plot(x = reg_df$time, y=reg_df$dnb_corrected_roman, xlab = "Days", ylab = "BRDF Corrected DNB")
predictedcounts <- predict(fit_corrected_roman,list(time=reg_df$time), time2=reg_df$time^2)
lines(reg_df$time, predictedcounts,  lwd = 3, xlab = "Time (s)", ylab = "Counts",col = "blue")

summary(fit_original)
summary(fit_corrected_roman)
summary(fit_corrected)



ylim=c(min(reg_df$dnb), max(reg_df$dnb))
g_secOrdReg_1 <- ggplot_g_second_order_regression(reg_df,"dnb", "Original DNB" )
g_secOrdReg_2 <- ggplot_g_second_order_regression(reg_df, "dnb_corrected", "Corrected DNB" )



#### Batch save ggplots
plots <- list(g_secOrdReg_1, g_secOrdReg_2)
invisible(
  lapply(
    seq_along(plots), 
    function(x) {
    ggsave(
      sprintf('fig.5_%s_ggplot_second_order_regression_%s_%s.tif',x, name, 2018),
      plot = plots[[x]],
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 9,
      height = 5,
      units = c("cm"),
      dpi = 600,
      limitsize = TRUE
)}))




###### End of second order regression #######################





# generate dates sequence for one year
days365<-seq(min(mydates), max(mydates), "days")

# combine mean (or SoL), full dates year
my_dataframe <-data.frame(mydates=days365) %>%
    dplyr::left_join(data.frame(mydates=mydates,
                                original=DNB_na.original_stat,
                                median_shift_corrected=DNB_na.original_corrected_stat
                                #na.intepolated_corrected=DNB_na.intepolated_corrected_stat
                                #na.first_previous_corrected=DNB_na.first_previous_corrected_stat
                                ))


# convert to long format for ggplot
data_long <- tidyr::gather(my_dataframe, dataset, value, original:median_shift_corrected) %>%
    dplyr::filter(!is.na(value))




# ggplot ------------------------------------------------------------------

days365 <- seq(min(mydates),max(mydates),by="day")
mi <- getMoonIllumination(date = days365, 
                                 keep = c("fraction", "phase", "angle")) %>% 
  mutate(state = ifelse(phase >=0.4 & phase<=0.6, "Full", NA))

runs<-rle(as.vector(mi$state))
myruns = which(runs$values == "Full")
runs.lengths.cumsum = cumsum(runs$lengths)
ends = runs.lengths.cumsum[myruns]
newindex = ifelse(myruns>1, myruns-1, 0)
starts = runs.lengths.cumsum[newindex] + 1
if (0 %in% newindex) starts = c(1,starts)
starts_moon_dates <- days365[starts]
starts_end_dates  <- days365[ends]

dateRanges_moons <- data.frame(from = starts_moon_dates,
                         to = starts_end_dates) %>% na.omit()


Sys.setlocale("LC_ALL", 'en_US.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")


# full year with original and corrected data and FULL MOONS on vetrical bars
(
  g_year_fullmoons  <-
    ggplot_Original_Corrected_Data(
      data_long,
      dateRanges_moons,
      "0.4<=moon phase<= 0.6,\n(0.5=Full Moon)"
    )
) 

  # a PERIOD of year with original and corrected data and  FULL MOONS on vetrical bars
(g_month_fullmoons <-
    g_year_fullmoons + coord_cartesian(xlim = (c(
      as.Date("2018-06-25"), as.Date("2018-07-31")
    )))) #zoom to date range




### save all plots as tiff

plots <- list(g_year_fullmoons)
invisible(
  lapply(
    seq_along(plots), 
    function(x) {
    ggsave(
      sprintf('fig.6_ggplot_original_corrected_%s_%s_%s.tif',name,x,2018),
      plot = plots[[x]],
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 25,
      height = 10,
      units = c("cm"),
      dpi = 1200,
      limitsize = TRUE
) }))

ggsave(
      sprintf('fig.4_ggplot_original_corrected_%s_%s_%s.tif',name,"month",2018),
      plot = g_month_fullmoons,
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 20,
      height = 10,
      units = c("cm"),
      dpi = 600,
      limitsize = TRUE
)

# ..............................................................................





# histograms New/Full moons -----------------------------------------------



newmoon_date <-newmoons[1]
g1 <- ggplot_hist(DNB_na.original, newmoon_date,sprintf("New Moon (original), %s",newmoon_date))
g2 <- ggplot_hist(DNB_na.original_corrected, newmoon_date,sprintf("New Moon (median shift corrected), %s",newmoon_date))


full_date <-fullmoons[1] #fullmoons[5]
g3 <- ggplot_hist(DNB_na.original, full_date,sprintf("Full Moon (original), %s",full_date))
g4 <- ggplot_hist(DNB_na.original_corrected, full_date,sprintf("Full Moon (median shift corrected), %s",full_date))



### Save all plots as tiff
plots <- list(g1, g2, g3, g4)
invisible(
  lapply(
    seq_along(plots), 
    function(x) {
   ggsave(
    sprintf('hist_g%s_%s.tif',x, 2018),
    plot = plots[[x]],
    device = "tiff",
    path = here::here("output", name),
    scale = 1,
    width = 9,
    height = 5,
    units = c("cm"),
    dpi = 600,
    limitsize = TRUE
)
       }))

