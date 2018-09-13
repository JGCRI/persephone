#' A function to convert input changes in T, changes in W, and values of C to be standardized with C3MP standaradization
#' @param ctwdataframe A dataframe of C, T, W values to be converted
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @importFrom tidyr gather spread
#' @importFrom rethinking map precis
#' @export

convert_ctw <- function(ctwdataframe){

  # Some inputs to create the converting function
  # CTW sensitivity tests
  tests <- read.csv(file="data/c3mpdata/tests.csv",head = TRUE, sep = ",")
  tests$Test <- paste0("test",seq(from=1, to=99, length=99))



  # baseline values
  CO2b <- 360
  tempb <- 0
  precipb <- 1


  tests %>%
    summarize(T.mean = mean(Temp),
              T.sd = sd(Temp),
              W.mean = mean(Precip),
              W.sd = sd(Precip),
              C.mean = mean(CO2),
              C.sd = sd(CO2)) %>%
    mutate(baseT.s = (tempb - T.mean) / T.sd,
           baseW.s = (precipb - W.mean) / W.sd,
           baseC.s = (CO2b - C.mean) / C.sd) ->
    standardize.CTW.vals



  # converted base values
  Tb <- standardize.CTW.vals$baseT.s
  Pb <- standardize.CTW.vals$baseW.s
  Cb <- standardize.CTW.vals$baseC.s

  ctwdataframe %>%
    #standardize
    mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
           Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
           CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) %>%
    mutate(T.1 = Temp - Tb,
           T.2 = (Temp - Tb)^2,
           W.1 = Precip - Pb,
           W.2 = (Precip - Pb)^2,
           C.1 = CO2 - Cb,
           C.2 = (CO2 - Cb)^2,
           T.W = (Temp - Tb) * (Precip - Pb),
           T.C = (Temp - Tb) * (CO2 - Cb),
           W.C = (Precip - Pb) * (CO2 - Cb),
           T.W.C = (Temp - Tb) * (Precip - Pb) * (CO2 - Cb),
           T.W.2 = (Temp - Tb) * (Precip - Pb)^2,
           T.C.2 = (Temp - Tb) * (CO2 - Cb)^2,
           T.2.W = (Temp - Tb)^2 * (Precip - Pb),
           T.2.C = (Temp - Tb)^2* (CO2 - Cb),
           W.2.C = (Precip - Pb)^2 * (CO2 - Cb),
           W.C.2 = (Precip - Pb) * (CO2 - Cb)^2,
           T.3 = (Temp - Tb)^3,
           W.3 = (Precip - Pb)^3,
           C.3 = (CO2 - Cb)^3)

}
