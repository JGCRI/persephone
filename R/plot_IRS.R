#' A function to save 4 core Impact response surfaces for paper,
#' given a specific functional form for a specific ensemble
#' and the corresponding coefficients for mu and sigma.
#'
#' @param coefficients The point values of coefficients for mu
#' @param SDcoefficients The point values of coefficients for sigma
#' @param functionalform The name functional form for the mean of the likelihood
#' @param sdfunctionalform The name of the functional form for the standard deviation of the likelihood
#' @param baseCO2 The baseline value of CO2 concentration
#' @param baseTemp The baseline value of temperature changes
#' @param basePrecip The baseline value of precipitation changes
#' @param name The name of the ensemble under consideration
#' @param mesh Define how many steps in each direction to make the CTW grid for IRS plotting, defaults to 101
#' @param figpath The location to write the IRSplots, defaults to "figures/irs"
#' @param legendsave Whether to save the master legend of IRS changes, defaults to FALSE
#' @param high If you want to return the high response case instead of the mean, defaults to FALSE
#' @param low If you want to return the low response case instead of the mean, defaults to FALSE
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr gather spread
#' @export

plot_IRS <- function(coefficients, SDcoefficients, functionalform, sdfunctionalform,
         baseCO2, baseTemp, basePrecip, name, mesh = 101,
         figpath = "figures/irs", legendsave = FALSE, high = FALSE, low = FALSE){


  standardize.CTW.vals <- read.csv("data/c3mpdata/ctw_standardizingvalues.csv")




  #############################################################################################################
  # legend info for cohesive scale plots

  # color bar for yields
  col.vals <- read.csv("data/idmaps/brblHex.csv")
  cols <- paste0(col.vals$x[c(1,3,7,9,13,15,19,21,25,27,29,
                              32,
                              36,38,40,44,46,50,52,56,58,62,64)])
  br <-  c(seq(-105,-5,length = 11),-1.5,1.5, seq(5,105, length=11))
  maxInd <- length(br)
  vals <- plotrix::rescale(br, c(0,1))


  # black and white bar for sd
  # colsSD <- c("#000000","#191919","#333333","#4c4c4c","#666666","#7f7f7f",
  #              "#999999", "#b2b2b2","#cccccc","#e5e5e5", "#FFFFFF", "#FFFFFF")
  # brSD <- seq(0,55, length = 12)
  colsSD <- c("#000000","#333333","#666666",
              "#999999", "#b2b2b2",  "#FFFFFF", "#FFFFFF")
  brSD <- seq(0,60, length = 7)
  maxInd2 <- length(brSD)


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


  # grayscale  need to be able to subset teh shared colors and breaks for the SD grayscale
  grayscale <- function(data, grayscalebreaks = brSD) {

    which(brSD < min(data$sd)) -> lowInd
    if (length(lowInd)>0){
      lowInd <- max(lowInd)
    }else {
      lowInd <- 1
    }

    which(brSD > max(data$sd)) -> highInd
    if (length(highInd)>0){
      highInd <- min(highInd)
    }else {
      highInd <- maxInd2
    }

    cols2sd <- colsSD[lowInd:highInd]

    # actually cut the sd data according to breaks
    brks2 <-  cut(data$sd, breaks = grayscalebreaks)

    # the categegories we actually plot
    data$brks2 <- (brks2)

    # list return
    returndata <- list(data, cols2sd)

    return(returndata)

  }


  # precip labels
  precip_labels <- c("-50%", "-25%", "0%", "25%", "50%")
  precip_cuts <- c((0.5 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                   (0.75 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                   (1 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                   (1.25 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                   (1.5 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd)
  temp_labels <- c("-1", "0", "2", "4", "6", "8")
  temp_cuts <- c((-1 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 (0 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 (2 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 (4 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 (6 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 (8 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd)
  co2_labels <- c("330", "360", "450", "675", "900")
  co2_cuts <- c((330 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd,
                (360 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd,
                (450 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd,
                (675 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd,
                (900 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd)


  #############################################################################################################
  #############################################################################################################
  # Fixed CO2
  temp <- seq(-1,8, length = mesh)
  precip <- seq(0.5, 1.5, length = mesh)
  co2 <- baseCO2
  CTWgrid <- expand.grid(temp, precip, co2)
  colnames(CTWgrid) <- c("Temp", "Precip", "CO2")

  # evaluate the emulator on the CTW grid
  coefficients %>%
    func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
                       baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
                       functionalform = functionalform, sdfunctionalform = sdfunctionalform,
                       ensembleName = name,  stddev = TRUE) ->
    CTWgridE

  if(high){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange + sd) ->
      CTWgridE
  }
  if(low){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange - sd) ->
      CTWgridE
  }

  # COLOR BAR
  CTWgridE %>% colorbar() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2


  # GRAY SCALE
  CTWgridE %>% grayscale() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2sd


  # plot
  p2 <- ggplot(data = CTWgridE, aes(Temp, Precip, z = pctYieldChange)) +
    geom_tile(aes(fill = brks1)) +
    scale_fill_manual(values = cols2) +
    scale_x_continuous(expand=c(0,0), breaks = temp_cuts, labels = temp_labels) +
    scale_y_continuous(expand=c(0,0), breaks = precip_cuts, labels = precip_labels) +
    annotate("point",size = 4, shape = 22, color = "#ff0000", fill = "#ff0000", x = ((0 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd), y = ((1 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd))

  if(!high & !low){
    p2 <- p2 + geom_contour(data = CTWgridE, aes(Temp, Precip, z = sd, colour = brks2), bins = 3) +
      #scale_colour_gradient( low = "black", high = "white") +
      scale_colour_manual(values = cols2sd)
  }

  p2 <- p2 + theme(legend.position="none")+
    ggtitle(paste0(name, "; Fixed CO2 = baseline"))

  ggsave(paste0(figpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_", name, "_fixC.png"), p2, width=8, height=8, units="in")


  # if(individsites){
  #   # gridded site IRS
  #   sites <- unique(coefficients$site)
  #   figs <- list()
  #
  #   for(i in 1:length(sites)){
  #
  #     coefficients %>%
  #       filter(site == sites[i]) %>%
  #       mutate(site = "ensemble") %>%
  #       func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
  #                          baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
  #                          functionalform = unique(coef$func), stddev = FALSE) ->
  #       CTWgridE
  #
  #     # COLOR BAR
  #     CTWgridE %>% colorbar() -> A
  #     A[[1]] -> CTWgridE
  #     A[[2]] -> cols2
  #
  #
  #       ggplot(data = CTWgridE, aes(Temp, Precip, z = pctYieldChange)) +
  #            geom_tile(aes(fill = brks1)) +
  #            scale_fill_manual(values = cols2) +
  #            scale_x_continuous(expand=c(0,0), breaks = temp_cuts, labels = temp_labels) +
  #            scale_y_continuous(expand=c(0,0), breaks = precip_cuts, labels = precip_labels) +
  #            annotate("point",size = 4, shape = 22, color = "#ff0000", fill = "#ff0000", x = ((0 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd), y = ((1 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd))+
  #            ggtitle(paste0(sites[i], name, "; Fixed CO2 = baseline")) + theme(legend.position="none") ->
  #       p1
  #
  #
  #       ggsave(paste0(figpath, "/", unique(coef$func), "_", name, "_site", sites[i], "_fixC.png"), p1, width=8, height=8, units="in")
  #
  #   }
  # }


  #############################################################################################################
  # Fixed CO2 = 812
  temp <- seq(-1,8, length = mesh)
  precip <- seq(0.5, 1.5, length = mesh)
  co2 <- 812
  CTWgrid <- expand.grid(temp, precip, co2)
  colnames(CTWgrid) <- c("Temp", "Precip", "CO2")

  # evaluate the emulator on the CTW grid
  coefficients %>%
    func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
                       baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
                       functionalform = functionalform, sdfunctionalform = sdfunctionalform,
                       ensembleName = name,  stddev = TRUE) ->
    CTWgridE


  if(high){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange + sd) ->
      CTWgridE
  }
  if(low){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange - sd) ->
      CTWgridE
  }


  # COLOR BAR
  CTWgridE %>% colorbar() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2


  # GRAY SCALE
  CTWgridE %>% grayscale() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2sd


  # plot
  p2 <- ggplot(data = CTWgridE, aes(Temp, Precip, z = pctYieldChange)) +
    geom_tile(aes(fill = brks1)) +
    scale_fill_manual(values = cols2) +
    scale_x_continuous(expand=c(0,0), breaks = temp_cuts, labels = temp_labels) +
    scale_y_continuous(expand=c(0,0), breaks = precip_cuts, labels = precip_labels)

  if(!high & !low){
    p2 <- p2 + geom_contour(data = CTWgridE, aes(Temp, Precip, z = sd, colour = brks2), bins = 3) +
      #scale_colour_gradient( low = "black", high = "white") +
      scale_colour_manual(values = cols2sd)
    }

  p2 <- p2 + theme(legend.position="none")+
    ggtitle(paste0(name, "; Fixed CO2 = 812 ppm"))

  ggsave(paste0(figpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_",  name, "_fixC812.png"), p2, width=8, height=8, units="in")



  # if(individsites){
  #   # gridded site IRS
  #   sites <- unique(coefficients$site)
  #   figs <- list()
  #
  #   for(i in 1:length(sites)){
  #
  #     coefficients %>%
  #       filter(site == sites[i]) %>%
  #       mutate(site = "ensemble") %>%
  #       func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
  #                          baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
  #                          functionalform = unique(coef$func), stddev = FALSE) ->
  #       CTWgridE
  #
  #     # COLOR BAR
  #     CTWgridE %>% colorbar() -> A
  #     A[[1]] -> CTWgridE
  #     A[[2]] -> cols2
  #
  #
  #     ggplot(data = CTWgridE, aes(Temp, Precip, z = pctYieldChange)) +
  #       geom_tile(aes(fill = brks1)) +
  #       scale_fill_manual(values = cols2) +
  #       scale_x_continuous(expand=c(0,0), breaks = temp_cuts, labels = temp_labels) +
  #       scale_y_continuous(expand=c(0,0), breaks = precip_cuts, labels = precip_labels) +
  #       ggtitle(paste0(sites[i], name, "; Fixed CO2 = 812ppm")) + theme(legend.position="none") ->
  #       p1
  #
  #
  #     ggsave(paste0(figpath, "/", unique(coef$func), "_", name, "_site", sites[i], "_fixC812.png"), p1, width=8, height=8, units="in")
  #
  #   }
  # }


  ############################################################################################################
  # Fixed Temp
  temp <- baseTemp
  precip <- seq(0.5, 1.5, length = mesh)
  co2 <- seq(330,900, length = mesh)
  CTWgrid <- expand.grid(temp, precip, co2)
  colnames(CTWgrid) <- c("Temp", "Precip", "CO2")


  # evaluate the emulator on the CTW grid
  coefficients %>%
    func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
                       baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
                       functionalform = functionalform, sdfunctionalform = sdfunctionalform,
                       ensembleName = name,  stddev = TRUE) ->
    CTWgridE


  if(high){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange + sd) ->
      CTWgridE
  }
  if(low){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange - sd) ->
      CTWgridE
  }


  # COLOR BAR
  CTWgridE %>% colorbar() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2


  # GRAY SCALE
  CTWgridE %>% grayscale() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2sd


  # plot
  p1 <-ggplot(data = CTWgridE, aes(Precip, CO2, z = pctYieldChange)) +
    geom_tile(aes(fill = brks1)) +
    scale_fill_manual(values = cols2)+
    scale_x_continuous(expand=c(0,0), breaks = precip_cuts, labels = precip_labels) +
    scale_y_continuous(expand=c(0,0), breaks = co2_cuts, labels = co2_labels) +
    annotate("point",size = 4, shape = 22, color = "#ff0000", fill = "#ff0000", x = ((1 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd), y = ((360 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd))

  if(!high & !low){
    p1 <- p1 + geom_contour(data = CTWgridE, aes(Precip, CO2, z = sd, colour = brks2), bins = 3) +
      scale_colour_manual(values = cols2sd)
  }

  p1 <- p1 + theme(legend.position="none")  +
    ggtitle(paste0(name, "; Fixed Temp = baseline"))

  ggsave(paste0(figpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_",  name, "_fixT.png"), p1, width=8, height=8, units="in")


  # if(individsites){
  #   # gridded site IRS
  #   sites <- unique(coefficients$site)
  #   figs <- list()
  #
  #   for(i in 1:length(sites)){
  #
  #     coefficients %>%
  #       filter(site == sites[i]) %>%
  #       mutate(site = "ensemble") %>%
  #       func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
  #                          baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
  #                          functionalform = unique(coef$func), stddev = F) ->
  #       CTWgridE
  #
  #     # COLOR BAR
  #     CTWgridE %>% colorbar() -> A
  #     A[[1]] -> CTWgridE
  #     A[[2]] -> cols2
  #
  #
  #     ggplot(data = CTWgridE, aes(Precip, CO2, z = pctYieldChange)) +
  #       geom_tile(aes(fill = brks1)) +
  #       scale_fill_manual(values = cols2)+
  #       scale_x_continuous(expand=c(0,0), breaks = precip_cuts, labels = precip_labels) +
  #       scale_y_continuous(expand=c(0,0), breaks = co2_cuts, labels = co2_labels) +
  #       annotate("point",size = 4, shape = 22, color = "#ff0000", fill = "#ff0000", x = ((1 - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd), y = ((360 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd))+
  #       ggtitle(paste0(sites[i], name, "; Fixed CO2 = baseline")) + theme(legend.position="none") ->
  #       p1
  #
  #
  #     ggsave(paste0(figpath, "/", unique(coef$func), "_", name, "_site", sites[i], "_fixT.png"), p1, width=8, height=8, units="in")
  #
  #
  #   }
  # }

  ############################################################################################################
  # Fixed Precip
  temp <- seq(-1,8, length = mesh)
  precip <- basePrecip
  co2 <- seq(330,900,length=mesh)
  CTWgrid <- expand.grid(temp, precip, co2)
  colnames(CTWgrid) <- c("Temp", "Precip", "CO2")

  # evaluate the emulator on the CTW grid
  coefficients %>%
    func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
                       baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
                       functionalform = functionalform, sdfunctionalform = sdfunctionalform,
                       ensembleName = name,  stddev = TRUE) ->
    CTWgridE



  if(high){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange + sd) ->
      CTWgridE
  }
  if(low){
    CTWgridE %>%
      mutate(pctYieldChange = pctYieldChange - sd) ->
      CTWgridE
  }


  # COLOR BAR
  CTWgridE %>% colorbar() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2


  # GRAY SCALE
  CTWgridE %>% grayscale() -> A
  A[[1]] -> CTWgridE
  A[[2]] -> cols2sd


  # plot
  p3 <- ggplot(data = CTWgridE, aes(Temp, CO2, z = pctYieldChange)) +
    geom_tile(aes(fill = brks1)) +
    scale_fill_manual(values = cols2)+
    scale_x_continuous(expand=c(0,0), breaks = temp_cuts, labels = temp_labels) +
    scale_y_continuous(expand=c(0,0), breaks = co2_cuts, labels = co2_labels) +
    # add baseline point:
    annotate("point",size = 2, shape = 22, color = "#ff0000", fill = "#ff0000", x = ((0 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd), y = ((360 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd))

  if(!high & !low){
    p3 <- p3 + geom_contour(data = CTWgridE, aes(Temp, CO2, z = sd, colour = brks2), bins = 3) +
      scale_colour_manual(values = cols2sd)
    }

  p3 <- p3 + theme(legend.position="none")+
    ggtitle(paste0(name, "; Fixed Precip = baseline"))

  ggsave(paste0(figpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_",  name, "_fixW.png"), p3, width=8, height=8, units="in")


  # if(individsites){
  #   # gridded site IRS
  #   sites <- unique(coefficients$site)
  #   figs <- list()
  #
  #   for(i in 1:length(sites)){
  #
  #     coefficients %>%
  #       filter(site == sites[i]) %>%
  #       mutate(site = "ensemble") %>%
  #       func_form_evaluate(coeff = ., SDcoeff = SDcoefficients, ctwvalues = CTWgrid,
  #                          baseC = baseCO2, baseT = baseTemp, baseP = basePrecip,
  #                          functionalform = unique(coef$func), stddev = FALSE) ->
  #       CTWgridE
  #
  #     # COLOR BAR
  #     CTWgridE %>% colorbar() -> A
  #     A[[1]] -> CTWgridE
  #     A[[2]] -> cols2
  #
  #
  #     ggplot(data = CTWgridE, aes(Temp, CO2, z = pctYieldChange)) +
  #       geom_tile(aes(fill = brks1)) +
  #       scale_fill_manual(values = cols2)+
  #       scale_x_continuous(expand=c(0,0), breaks = temp_cuts, labels = temp_labels) +
  #       scale_y_continuous(expand=c(0,0), breaks = co2_cuts, labels = co2_labels) +
  #       # add baseline point:
  #       annotate("point",size = 2, shape = 22, color = "#ff0000", fill = "#ff0000", x = ((0 - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd), y = ((360 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd)) +
  #     ggtitle(paste0(sites[i], name, "; Fixed CO2 = baseline")) + theme(legend.position="none") ->
  #       p1
  #
  #
  #     ggsave(paste0(figpath, "/", unique(coef$func), "_", name, "_site", sites[i], "_fixW.png"), p1, width=8, height=8, units="in")
  #
  #
  #   }
  # }

  ############################################################################################################
  ###### Making shared full legend

  # function from stack exchange to extract a legend from a plot
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}


  if(legendsave){
    # create a dummy plot with the full legend
    temp <- seq(-1,8, length = mesh)
    precip <- seq(0.5, 1.5, length = mesh)
    dum <- expand.grid(temp, precip)
    dum %>%
      tibble::as_tibble() %>%
      mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
             Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd) ->
      dum

    pctYieldChange <- as.data.frame(seq(-104,104, length = nrow(dum)))
    sd <- as.data.frame(seq(1 , max(brSD)-1, length = nrow(dum)))
    dum <- bind_cols(dum, pctYieldChange, sd)
    colnames(dum) <- c("Temp", "Precip", "pctYieldChange", "sd")

    # COLOR BAR
    dum %>% colorbar() -> A
    A[[1]] -> dum
    A[[2]] -> cols2

    # GRAYSCALE
    dum %>% grayscale() -> A
    A[[1]] -> dum
    A[[2]] -> cols2sd


    p_fake <- ggplot(data=dum, aes(Temp, Precip, z = pctYieldChange)) +
      geom_tile(aes(fill = brks1)) +
      scale_fill_manual(values = cols2)+
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme(legend.position = "right") +
      guides(fill = guide_legend(ncol = 1, reverse = TRUE, title="Ensemble % change in Y from baseline CTW "))
    p_fake

    p_fake2 <- ggplot() + geom_contour(data = dum, aes(Temp, Precip, z = sd, colour = brks2)) +
      scale_colour_manual(values = cols2sd)+
      guides(colour = guide_legend(reverse = TRUE, override.aes=list(size=10), keyheight=0, title="sd across sites for fixed CTW"))+
      theme(legend.key=element_blank())
    p_fake2

    # extract legend from fake figure
    # could do the legend from single fake figure; this allows more flexibility in position though
    plot_legend <- g_legend(p_fake)
    plot_legend2 <- g_legend(p_fake2)

    # arrange
    p_legend <- gridExtra::grid.arrange(plot_legend, plot_legend2, ncol=2)

    ggsave(paste0(figpath, "/","legend.png"), p_legend, width=8, height=8, units="in")
  }

} # end function
