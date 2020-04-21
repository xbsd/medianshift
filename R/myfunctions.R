

generate_daily_dnb <- function(dnb_files){

    message("Generate Stack")
    #original_dnb_stack <- raster::stack(original_dnb_files)
    dnb_stack <- raster::stack(dnb_files)


    dates <- Reduce(c,lapply(X=dnb_files, FUN=function(x) {
      date = strsplit(strsplit(x,"#")[[1]][2],"-")[[1]]
      year = date[1]
      julianday = date[2]
      as.Date(as.integer(julianday)-1, origin=as.Date(sprintf("%s-01-01",year)))
    }))



    ## Find Duplicate dates
    duplicated_dates <- dates[duplicated(dates)]
    if (length(duplicated_dates)>0){
        duplicated_dates_indexes <- which(dates %in% duplicated_dates) #find indexes of duplicated dates
        duplicated_dnb_stack <- subset(dnb_stack, duplicated_dates_indexes) %>% setZ(dates[dates %in% duplicated_dates])

    #### mean composite for duplicate dates
    tryCatch({
      system.time({
        no_cores <- parallel::detectCores() - 1
        raster::beginCluster(no_cores)
        myFun <- function(x, ...) {
          mean(x, na.rm=T)
        }
        daily_DNB_apply <- raster::clusterR(duplicated_dnb_stack, raster::zApply, args = list(by=as.Date,fun = myFun))
        raster::endCluster()})
    }, error = function(e) {
      raster::endCluster()
      return(e)
    }, finally = {
      try(raster::endCluster())
    })
    # myFun <- function(x, ...) {
    #   mean(x, na.rm=T)
    # }
   #daily_DNB_apply <- raster::zApply(duplicated_dnb_stack, by=as.Date,fun = myFun)

    #merge unique dates with mean composite stack of duplicate dates
    dnb_stack_merged <- dnb_stack %>%
      dropLayer(duplicated_dates_indexes) %>% #drop duplicate tiffs
      raster::stack(daily_DNB_apply) %>% #merge with   mean composite stack
      setZ(as.Date(
        c(setdiff(dates,dates[duplicated_dates_indexes]), # unique dates
          duplicated_dates), #duplicate dates
        origin ="1970-01-01 UTC" ))


    #reorder layers and reassign dates to names and z values
    z_values <- getZ(dnb_stack_merged)
    daily_DNB <- dnb_stack_merged %>%
      subset(order(z_values)) %>%
      setZ(z_values[order(z_values)]) %>%
      setNames(z_values[order(z_values)])
    } else
        {
        daily_DNB <-dnb_stack %>%
            setZ(dates)%>%
            setNames(dates)
    }



    return(daily_DNB)
}





correct_DNB_30Days <- function(dnb_stack, newmoons_dates) {
    #browser()
    mydates <-
        as.Date(names(raster::stack(dnb_stack)), "X%Y.%m.%d")
    #phase <- lunar.phase(mydates)



    dnb_stack_df <-
        data.frame(mydates = mydates, month = lubridate::month(mydates)) %>% # month=lubridate::month(mydates)
        dplyr::mutate(meandnb = my_cellStats(dnb_stack, median, na.rm = T))  %>% as_tibble() %>% #%>% #todo: check median??
        mutate(isnewmoon = case_when(mydates %in% newmoons_dates ~ T))

    dnb_stack_df <-
        dnb_stack_df %>%  dplyr::filter(isnewmoon == T) %>% left_join(dnb_stack_df, by =
                                                                          ("month" = "month")) %>%
        dplyr::select(mydates.x, meandnb.x, meandnb.y) %>%
        dplyr::rename(
            mydates = mydates.x,
            newmoon_mediandnb = meandnb.x ,
            mediandnb = meandnb.y
        ) %>%
        dplyr::mutate(diff = mediandnb - newmoon_mediandnb) %>%
        dplyr::mutate(diff = replace(diff, diff < 0, 0)) %>%  #todo
        tibble::as_tibble()

    diff_filled_dnb <-dnb_stack_df$diff

    if (length(diff_filled_dnb)  != nlayers(dnb_stack)){
        stop("error at nlayers VS lenght")
    }
    dnb_stack_corrected <- dnb_stack - diff_filled_dnb #todo, is this correct?
    #dnb_stack_corrected[dnb_stack_corrected < 0] <- 0
    dnb_stack_corrected <- calc(dnb_stack_corrected,function(x) {x[x <0 ] <-0; return(x)}) %>% setNames(names(dnb_stack))


    return(dnb_stack_corrected)
}






mean_diff <- function(x){
    mean(abs(diff(x)))
}





ggplot_g_second_order_regression <- function(df, dnb_col, label){
  dnb_col <- sym(dnb_col)
  ggplot(df, aes(time, !!dnb_col)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(
      method = "lm",
      formula = y ~ poly(x, 2),
      colour = "black",
      size = 0.5
    ) +
    ylim(ylim) +
    ylab(label) +
    xlab("Days") +
    ggtitle(sprintf("%s, second order regression fit", label)) +
    theme_gray(base_size = 7)
}




ggplot_Original_Corrected_Data <- function(data_long,dateRanges, bar_legend ){
  ggplot2::ggplot(data = data_long) +
    scale_fill_manual('',
                      values = 'grey',
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    geom_rect(
        data = dateRanges,
        aes(
            xmin = from,
            xmax = to,
            ymin = -Inf,
            ymax = Inf#,
            #fill = bar_legend #uncomment to show in legend
        ),
        alpha = 0.5
    ) +
    #ggtitle(sprintf("%s, 2018, Mean VIIRS DNB Radiance",name)) +
    geom_line( data = data_long,aes(x = as.Date(mydates , "X%Y.%m.%d"), y =value,linetype = dataset), size=0.3) +
    scale_linetype_manual(name = NULL,#"DNB Datasets:",
                        values = c("median_shift_corrected" = "solid", "original" = "dashed"),
                        labels = c("Median shift", "Original")) +
    xlab("Date") +
    ylab("DNB") +
    scale_y_continuous(limit = c(0, max(data_long$value))) +theme_gray(base_size = 10)+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position=c(0.2,0.8))
}



ggplot_hist <- function(dnb, newmoon_date, label) {
    #browser()
    newmoonvals <- dnb[[match(newmoon_date, mydates)]][]
    #newmoonvals[is.na(newmoonvals)]<-0 #κάνω το NA, μηδέν, το απορρίψαμε με τον Σταθάκη.
    #newmoonvals <-
    #newmoonvals[newmoonvals <= 15] #todo κρατάω τιμές μέχρι 15 nano
    dat <- data.frame(dnb = newmoonvals)
    my_mean <- round(mean(dat$dnb , na.rm = T), 2)
    my_median <- round(median(dat$dnb, na.rm = T), 2)

    g <- ggplot(dat ,   aes (x = dnb))   +
        ggtitle(sprintf('%s, \nMean:%s, Median:%s', label, my_mean, my_median)) +
        geom_histogram (binwidth = 1 ,
                        colour = "black" ,
                        fill = "white", size = 0.3)   +
        # geom_vline (aes (xintercept = my_mean, group = "mean"),
        #             linetype = "dashed" ,
        #             size = 0.4) +
        geom_vline (
            aes (xintercept = my_median, group = "median"),
            linetype = "dotted" ,
            size = 0.5
        ) +
        labs(group = '') +
        scale_x_continuous(breaks = seq(0,40,5), limits = c(0,40), expand = c(0, 0))+
        scale_y_continuous(expand = c(0,0),limits = c(0,300), breaks = seq(0,300,50))+
        #coord_cartesian(ylim = c(0, 100), expand = c(0, 0)) + # αν θέλω να περιορίσω το ευρος που φτάνει το ylim (ή xlim)

        #xlim(0, 40)+
        #ylim(0, 100)+
        #
        #scale_x_continuous(breaks = round(seq(min(dat$dnb, na.rm = T), max(dat$dnb, na.rm = T), by = 5), 1)) +

        # scale_y_continuous(breaks = round(seq(0,255, by = 50), 1)) +
        # scale_x_continuous(breaks = round(seq(0,100, by = 5), 1)) +
        theme_gray(base_size = 7)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text.x =
                  element_text(
                      size  = 5,
                      angle = 0,
                      hjust = 1,
                      vjust = 1
                  ))

  #           scale_x_continuous(breaks = round(seq(min(dat$dnb, na.rm = T), max(dat$dnb, na.rm = T), by = 5), 1), labels = function(labels) {
  #   fixedLabels <- c()
  #   for (l in 1:length(labels)) {
  #     fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
  #   }
  #   return(fixedLabels)
  # })
    return(g)
}

moon_phases <- function(mydates){
  ##################################  *lunar phase, for visualization ##########################################################################

  days365<-seq(min(mydates), max(mydates), "days")
  mymoon<-lunar.phase(days365, shift = 0, name=T)


  #extract full moon days
  #http://masterr.org/r/how-to-find-consecutive-repeats-in-r/
  runs<-rle(as.vector(mymoon))
  myruns = which(runs$values == "Full")
  runs.lengths.cumsum = cumsum(runs$lengths)
  ends = runs.lengths.cumsum[myruns]
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)

  #extract NEW moon days
  #http://masterr.org/r/how-to-find-consecutive-repeats-in-r/
  runs2<-rle(as.vector(mymoon))
  myruns2 = which(runs2$values == "New")
  runs.lengths.cumsum2 = cumsum(runs2$lengths)
  ends2 = runs.lengths.cumsum2[myruns2]
  newindex2 = ifelse(myruns2>1, myruns2-1, 0)
  starts2 = runs.lengths.cumsum2[newindex2] + 1
  if (0 %in% newindex2) starts2 = c(1,starts2)
  my_list <- list(fullmoon_starts = starts,
                  fullmoon_ends = ends,
                  newmoon_starts = starts2,
                  newmoon_ends = ends2,
                  days365 =days365)
  return(my_list)

}



