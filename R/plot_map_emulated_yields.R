#' A function to save 4 core Impact response surfaces for paper,
#' given a specific functional form for a specific ensemble
#' and the corresponding coefficients for mu and sigma.
#'
#' @param crop The crop to calculate changes for, including looking up coefficients
#' @param coefpath The location of the coefficients, defaults to "data/emulatorcoefficients"
#' @param ctwpath The location of the driving CTW data, defaults to "data/climatedata/hadgem2es"
#' @param mapsavepath The location to save the output maps, defaults to "figures/maps"
#' @param plotyear The year to plot, defaults to 2085
#' @param plotirr The type of crop to plot, RFD or IRR. Defaults to RFD
#' @param func The name of the functional form for mu, for looking up coefficients
#' @param sdfunc The name of the functional form for sigma, for looking up coefficients
#' @param high If you want to return the high response case instead of the mean, defaults to FALSE
#' @param low If you want to return the low response case instead of the mean, defaults to FALSE
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr gather spread
#' @export

plot_map_emulated_yields <- function(crop,
                                     coefpath = "data/emulatorcoefficients",
                                     ctwpath = "data/climatedata/hadgem2es",
                                     mapsavepath = "figures/maps",
                                     plotyear = 2085, plotirr = "RFD",
                                     func, sdfunc,
                                     high = FALSE, low = FALSE){


  standardize.CTW.vals <- read.csv("data/c3mpdata/ctw_standardizingvalues.csv")

  ##########################################################################################################################################
  #### Helpful Plotting pieces
  # get small area mask for plotting
  read.csv("data/seasondata/smallAreaMask_allcropsirr.csv") %>%
    tibble::as_tibble() %>%
    mutate_if(is.factor, as.character) %>%
    filter(crop == crop) ->
    smallareamask

  # full set of latlons
  latlons <- tibble::as_tibble(read.csv("data/seasondata/latlon_fullgrid.csv"))


  # countries for mapping
  countries <- map_data("world")

  maxgroup <- max(countries$group)

  countries %>%
    mutate(group = if_else(long > 180, maxgroup + 1, group),
           long = if_else(long > 180, long - 360, long)) ->
    countries1


  #############################################################################################################
  # legend info for cohesive scale plots

  # color bar for % change
  col.vals <- read.csv("data/idmaps/brblHex.csv")
  cols <- paste0(col.vals$x[c(1,3,7,9,13,15,19,21,25,27,29,
                              32,
                              36,38,40,44,46,50,52,56,58,62,64, 64)])
  br <-  c(seq(-105,-5,length = 11),-1,1, seq(5,105, length=11), 600)
  maxInd <- length(br)
  vals <- plotrix::rescale(br, c(0,1))


  # colorbar - need to be able to subset the shared colors and breaks
  colorbar <- function(data, colorbreaks = br){

    which(br < min(data$pctYieldChange)) -> lowInd
    if (length(lowInd)>0){
      lowInd <- max(lowInd)
    }else {
      lowInd <- 1
    }

    which(br > max(data$pctYieldChange)) -> highInd
    if (length(highInd)>0){
      highInd <- min(highInd)
    }else {
      highInd <- maxInd
    }

    ### need to fix this here, below colorbar functions, and in plot irs
    br2 <- br[lowInd:highInd]
    cols2 <- cols[lowInd:highInd]
    vals2 <- vals[lowInd:highInd]

    # actually cut the data according to breaks
    brks1 <-  cut(data$pctYieldChange, breaks = colorbreaks)

    # the categegories we actually plot
    data$brks1 <- (brks1)

    # list return
    returndata <- list(data, cols2)

    return(returndata)
  }


  # colorbar - need to be able to subset the shared colors and breaks
  br_precip <- (1/100)* br + 1
  br_precip <- (br_precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd
  colorbar_precip <- function(data, colorbreaks = br_precip){

    # data %>%
    #   mutate(Precip1 = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd) ->
    #   data

    which(colorbreaks < min(data$Precip)) -> lowInd
    if (length(lowInd)>0){
      lowInd <- max(lowInd)
    }else {
      lowInd <- 1
    }

    which(colorbreaks > max(data$Precip)) -> highInd
    if (length(highInd)>0){
      highInd <- min(highInd)
    }else {
      highInd <- maxInd
    }

    br2 <- colorbreaks[lowInd:highInd]
    cols2 <- cols[lowInd:highInd]
    vals2 <- vals[lowInd:highInd]

    # actually cut the data according to breaks
    brksp <-  cut(data$Precip, breaks = colorbreaks)

    # the categegories we actually plot
    data$brksp <- (brksp)

    # list return
    returndata <- list(data, cols2)

    return(returndata)
  }




  # color bar for Temperature Change
  cols_temp <- rev(c("#B90000", "#BD0F0F", "#C11F1F", "#C62F2F", "#CA3F3F", "#CE4F4F", "#D35F5F",
                 "#D76F6F", "#DC7F7F", "#E08F8F", "#E49F9F", "#E9AFAF", "#EDBFBF", "#F1CFCF",
                 "#F6DFDF", "#FAEFEF", "#FFFFFF", "#DCDCDC", "#B9B9B9"))
  br_temp <- seq(-1, 8, by = 0.5)
  br_temp <- (br_temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd
  maxInd <- length(br_temp)

  colorbar_temp <- function(data, colorbreaks = br_temp){

    # data %>%
    #   mutate(Temp1 = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd) ->
    #   data

    which(colorbreaks < min(data$Temp)) -> lowInd
    if (length(lowInd)>0){
      lowInd <- max(lowInd)
    }else {
      lowInd <- 1
    }

    which(colorbreaks > max(data$Temp)) -> highInd
    if (length(highInd)>0){
      highInd <- min(highInd)
    }else {
      highInd <- maxInd
    }

    ### need to fix this here, below colorbar functions, and in plot irs
    br2 <- colorbreaks[lowInd:highInd]
    cols2 <- cols_temp[lowInd:highInd]
    vals2 <- vals[lowInd:highInd]

    # actually cut the data according to breaks
    brkst <-  cut(data$Temp, breaks = colorbreaks)

    # the categegories we actually plot
    data$brkst <- (brkst)

    # list return
    returndata <- list(data, cols2)

    return(returndata)
  }




  ##########################################################################################################################################
  #### CTW Data

  # read CO2 data
  read.csv("data/climatedata/magicc_rcp8p5_co2.csv") %>% # hits CO2 360 around 1998-1999 with Tgav 0.872029 - 0.896249
    tibble::as_tibble() %>%
    # mutate(Ca = Ca + 42.5) %>% # adjust so baseline CO2 approx 360 when global temp change about 0
    # filter(year %in% c(maxhistyear, future_years)) %>%
    select(year, value) %>%
    rename(CO2 = value) %>%
    filter(year >= 2070, year <= 2100) %>%
    summarize(CO2 = mean(CO2)) %>%
    mutate(year = plotyear)  %>%
    mutate(CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) %>%
    mutate_if(is.factor, as.character()) ->
    CO2

  baseCO2 <- (360 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd
  # CO2 %>% filter(year == maxhistyear) %>% dplyr::select(co2) %>% as.double -> baseCO2

  # # Plot and save CO2 time series
  # pCO2 <- ggplot(CO2, aes(x = year, y = CO2)) + geom_line() + ggtitle("Global CO2 concentration")
  # ggsave(paste0(mapsavepath, "/GlobalCO2.png"), width = 12, height = 9, units = "in")

  # create CTW for each crop, with some extra identifying info
  read.csv(paste0(ctwpath, "/TempPrecip_smoothed_future_years_", crop, "_2070_2100.csv")) %>%
    tibble::as_tibble() %>%
    mutate_if(is.factor, as.character()) %>%
    filter(irr == plotirr) %>%
    left_join(CO2, by = "year") %>%
    left_join(smallareamask, by = c("latgrid", "longrid", "crop", "irr")) %>%
    mutate(deltaT = deltaT * areamask,
           deltaP = if_else(areamask == 0, 1, deltaP),
           CO2 = if_else(areamask == 0, baseCO2, CO2)) %>%
    filter(areamask == 1,
           !((latgrid == 213 & longrid == 719) | (latgrid == 213 & longrid == 720))) %>%
    rename(Temp = deltaT,
           Precip = deltaP) %>%
    filter(year == plotyear) ->
    temp_precip1



  ##########################################################################################################################################
  #### Coeff Data
  coefflist  <- list.files(path = coefpath, full.names=TRUE, recursive=FALSE)
  indfiles1 <- grep(paste0(func, "_mu_", sdfunc, "_sigma_", crop, "_", plotirr), coefflist)
  coefflist  <- coefflist[indfiles1]

  indfiles2 <- grep("mu_coeff", coefflist)
  coefflist.mu  <- coefflist[indfiles2]

  indfiles3 <- grep("sigma_coeff", coefflist)
  coefflist.sigma  <- coefflist[indfiles3]


  coef <- tibble::tibble()
  invisible(foreach(f = coefflist.mu) %do% {
    # get crop irrigation, latitude info from the filename
    f %>%
      tibble::as_tibble() %>%
      tidyr::extract(value, into = c("sname"), ".+(\\/.+)$") %>%
      mutate(sname = substr(sname, 2, nchar(sname)-4)) %>%
      separate(sname, c("func","t1", "sdfunc", "t2", "crop", "irr", "latband", "t3"), sep = "_") %>%
      dplyr::select(func, sdfunc, crop, irr, latband) ->
      f_name_id_info


    read.csv(f) %>%
      tibble::as_tibble() %>%
      dplyr::select(-StdDev) %>%
      # spread(mu.param, Mean) %>%
      repeat_add_columns(f_name_id_info) %>%
      bind_rows(coef, .) ->
      coef
  })


  sdcoef <- tibble::tibble()

  invisible(foreach(f = coefflist.sigma) %do% {
    # get crop irrigation, latitude info from the filename
    f %>%
      tibble::as_tibble() %>%
      tidyr::extract(value, into = c("sname"), ".+(\\/.+)$") %>%
      mutate(sname = substr(sname, 2, nchar(sname)-4)) %>%
      separate(sname, c("func","t1", "sdfunc", "t2", "crop", "irr", "latband", "t3"), sep = "_") %>%
      dplyr::select(func, sdfunc, crop, irr, latband) ->
      f_name_id_info


    read.csv(f) %>%
      tibble::as_tibble() %>%
      dplyr::select(-StdDev) %>%
      #spread(sigma.param, Mean) %>%
      repeat_add_columns(f_name_id_info) %>%
      bind_rows(sdcoef, .) ->
      sdcoef
  })



  ##########################################################################################################################################
  ### Add latband information to the temp_precip data for joining to appropriate coefficients

  # overall lat bands info
  temp_precip1 %>%
    mutate(latcol = lati) %>%
    add_lat_bands(tropic_cut = 30, mid_cut = 70) %>%
    select(-latcol) %>%
    rename(latband = latBand) ->
    temp_precip


  if(unique(temp_precip$crop) == "Corn"){
    temp_precip %>%
      mutate(crop = "Maize") ->
      temp_precip
  }

  # Bound to paramater boundaries so our emulators work, then standardize
  temp_precip %>%
    mutate(Precip = if_else(Precip > 1.5, 1.5, if_else(Precip < 0.5, 0.5, Precip)),
           Temp = if_else(Temp > 8, 8, Temp)) %>%
    mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
           Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd) ->
    temp_precip


  ##########################################################################################################################################
  ## emulate


  # # emulate for each gridcell-latband-crop-irr

    func_form_evaluate(coeff = filter(coef, latband == "mid"), SDcoeff = filter(sdcoef, latband == "mid"),
                       ctwvalues = filter(temp_precip, latband == "mid"),
                       functionalform = func, sdfunctionalform = sdfunc,
                       ensembleName = name,  stddev = TRUE, formaps = TRUE) ->
      # mutate(Temp = Temp * standardize.CTW.vals$T.sd + standardize.CTW.vals$T.mean,
      #        Precip = Precip * standardize.CTW.vals$W.sd + standardize.CTW.vals$W.mean,
      #        CO2 = CO2 * standardize.CTW.vals$C.sd + standardize.CTW.vals$C.mean) ->
    yields.mid


    func_form_evaluate(coeff = filter(coef, latband == "tropic"), SDcoeff = filter(sdcoef, latband == "tropic"),
                       ctwvalues = filter(temp_precip, latband == "tropic"),
                       functionalform = func, sdfunctionalform = sdfunc,
                       ensembleName = name,  stddev = TRUE, formaps = TRUE) %>%
      bind_rows(yields.mid) ->
      yields1


  if(high){
    yields1 %>%
      mutate(pctYieldChange = pctYieldChange + sd) ->
      yields1
  }
  if(low){
    yields1 %>%
      mutate(pctYieldChange = pctYieldChange - sd) ->
      yields1
  }



  # colorbars
  yields1 %>% colorbar() -> A

  A[[1]] -> yields1
  A[[2]] %>% na.omit -> cols2

  temp_precip %>%
    dplyr::select(long, lati, latgrid, longrid, crop, irr, year, Temp, Precip)%>%
    colorbar_precip() ->
    A
  A[[1]] -> temp_precip
  A[[2]] -> colsp

  temp_precip %>%
    colorbar_temp() ->
    A
  A[[1]] -> temp_precip
  A[[2]] -> colst




  temp_precip ->
    CTW1


  # plot
  yields1 %>%
    # filter(year == plotyear, irr == plotirr) %>%
    mutate(pctYieldChange = if_else(pctYieldChange < -100, -100, pctYieldChange)) ->
    plotyields

  pyield <- ggplot(data = plotyields, aes(x = long, y = lati)) +
    geom_tile(aes(fill = brks1)) +
    scale_fill_manual(values = cols2) +
    geom_polygon(data = countries1, aes(x = long, y = lat, group = group), col = "black", lwd = 0.5, fill = NA) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position = "none", panel.background = element_rect(fill = "grey50", colour = "grey50"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    ylim(-65, 85)

    if(high){
      ggsave(paste0(mapsavepath, "/highYieldChange_", crop, "_", plotirr, "_", plotyear, ".png"), pyield, width = 17, height = 11, units = "in")
    }else if(low){
      ggsave(paste0(mapsavepath, "/lowYieldChange_", crop, "_", plotirr, "_", plotyear, ".png"), pyield, width = 17, height = 11, units = "in")
      } else{
        ggsave(paste0(mapsavepath, "/meanYieldChange_", crop, "_", plotirr, "_", plotyear, ".png"), pyield, width = 17, height = 11, units = "in")
      }






    pT <- ggplot(data = CTW1, aes(x = long, y = lati)) +
      geom_tile(aes(fill = brkst)) +
      scale_fill_manual(values = colst) +
      geom_polygon(data = countries1, aes(x = long, y = lat, group = group), col = "black", lwd = 0.5, fill = NA) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0))  +
      theme(legend.position = "none", panel.background = element_rect(fill = "grey50", colour = "grey50"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ylim(-65, 85)

    ggsave(paste0(mapsavepath, "/tempChange_", crop, "_", plotirr, "_", plotyear, ".png"), pT, width = 17, height = 11, units = "in")


    pP <- ggplot(data = CTW1, aes(x = long, y = lati)) +
      geom_tile(aes(fill = brksp)) +
      scale_fill_manual(values = colsp) +
      geom_polygon(data = countries1, aes(x = long, y = lat, group = group), col = "black", lwd = 0.5, fill = NA) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0))  +
      theme(legend.position = "none", panel.background = element_rect(fill = "grey50", colour = "grey50"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ylim(-65, 85)


    ggsave(paste0(mapsavepath, "/precipChange_", crop, "_", plotirr, "_", plotyear, ".png"), pP, width = 17, height = 11, units = "in")




    ### save temperature colorbar
    legendsave <- FALSE

    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}


    if(legendsave){
      # create a dummy plot with the full legend
      temp <- tibble::tibble(Temp = seq(-1,8, by = 0.5))


      # COLOR BAR
      temp %>% colorbar_temp() -> A
      A[[1]] -> temp
      A[[2]] -> colst



      p_fake <- ggplot(data=temp, aes(x = Temp, y = 1)) +
        geom_tile(aes(fill = brkst)) +
        scale_fill_manual(values = colst)+
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        theme(legend.position = "right") +
        guides(fill = guide_legend(ncol = 1, reverse = TRUE, title=""))
      p_fake



      # extract legend from fake figure
      # could do the legend from single fake figure; this allows more flexibility in position though
      plot_legend <- g_legend(p_fake)

      ggsave(paste0(mapsavepath, "/","legend_temp.png"), plot_legend, width=2, height=8, units="in")
    }


} # end function
