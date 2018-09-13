#' A function to evaluate yield changes in response to input CTW values,
#' given a specific functional form for a specific ensemble
#' and the corresponding coefficients for mu and sigma.
#'
#' @param coeff The point values of coefficients for mu
#' @param SDcoeff The point values of coefficients for sigma
#' @param ctwvalues The CTW values at which to evaluate
#' @param baseC The baseline value of CO2 concentration, defaults to C3MP baseline of 360 ppm
#' @param baseT The baseline value of temperature changes, defaults to C3MP baseline of 0C change
#' @param base P The baseline value of precipitation changes, defaults to C3MP baseline of 0% change
#' @param functionalform The name functional form for the mean of the likelihood
#' @param sdfunctionalform The name of the functional form for the standard deviation of the likelihood
#' @param ensembleName The name of the ensemble under consideration
#' @param stddev Whether to return the standard deviation of the yield estimate (ie whether to also evaluate sigma), defaults to TRUE
#' @param formaps Whether the outputs will be used for plotting maps, defaults to TRUE. Just changes some formatting.
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @importFrom tidyr gather spread
#' @export


func_form_evaluate <- function(coeff, SDcoeff, ctwvalues,  baseC = 360, baseT = 0, baseP = 1,
         functionalform, sdfunctionalform, ensembleName,  stddev = TRUE, formaps = FALSE){

  standardize.CTW.vals <- read.csv("data/c3mpdata/ctw_standardizingvalues.csv")


  if(!formaps) {
    ctwvalues %>% convert_ctw() -> ctwvalues
  }
  if(formaps) {
    Tb <- standardize.CTW.vals$baseT.s
    Pb <- standardize.CTW.vals$baseW.s
    Cb <- standardize.CTW.vals$baseC.s

    ctwvalues %>%
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
             C.3 = (CO2 - Cb)^3) ->
      ctwvalues
  }

  coeff %>%
   # filter(site == "ensemble") %>%
    select(mu.param, Mean) %>%
    spread(mu.param, Mean) %>%
    repeat_add_columns(ctwvalues) %>%
    mutate_if(is.factor, as.character) ->
    coefficients



  if(functionalform == "c3mp"){

    coefficients %>%
      mutate(pctYieldChange =  t * T.1
             + t2 * T.2
             + w * W.1
             + w2 * W.2
             + c * C.1
             + c2 * C.2
             + tw * T.W
             + tc * T.C
             + wc * W.C
             + twc * T.W.C) %>%
      # select(Temp, Precip, CO2,  pctYieldChange) %>%
      #distinct %>%
      mutate(func = paste0(functionalform)) ->
      outdata





    if(stddev & sdfunctionalform == "quadratic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C )) %>%
        select(Temp, Precip, CO2, sd) %>%
        #distinct %>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


    if(stddev & sdfunctionalform == "cubic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C
                        + twc * T.W.C
                        + tw2 * T.W.2
                        + tc2 * T.C.2
                        + t2w * T.2.W
                        + t2c * T.2.C
                        + w2c * W.2.C
                        + wc2 * W.C.2
                        + t3 * T.3
                        + w3 * W.3
                        + c3 * C.3)) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


  }# end If c3mp


  if(functionalform == "quadratic"){

    coefficients %>%
      mutate(pctYieldChange = t * T.1
             + t2 * T.2
             + w * W.1
             + w2 * W.2
             + c * C.1
             + c2 * C.2
             + tw * T.W
             + tc * T.C
             + wc * W.C) %>%
      # select(Temp, Precip, CO2,  pctYieldChange, tes) %>%
      # distinct%>%
      mutate(func = paste0(functionalform)) ->
      outdata




    if(stddev & sdfunctionalform == "quadratic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C )) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


    if(stddev & sdfunctionalform == "cubic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C
                        + twc * T.W.C
                        + tw2 * T.W.2
                        + tc2 * T.C.2
                        + t2w * T.2.W
                        + t2c * T.2.C
                        + w2c * W.2.C
                        + wc2 * W.C.2
                        + t3 * T.3
                        + w3 * W.3
                        + c3 * C.3)) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


  }# end If quadratic


  if(functionalform == "cubic"){
    coefficients %>%
      mutate(pctYieldChange = t * T.1
             + t2 * T.2
             + w * W.1
             + w2 * W.2
             + c * C.1
             + c2 * C.2
             + tw * T.W
             + tc * T.C
             + wc * W.C
             + twc * T.W.C
             + tw2 * T.W.2
             + tc2 * T.C.2
             + t2w * T.2.W
             + t2c * T.2.C
             + w2c * W.2.C
             + wc2 * W.C.2
             + t3 * T.3
             + w3 * W.3
             + c3 * C.3) %>%
      # select(Temp, Precip, CO2, pctYieldChange) %>%
      # distinct%>%
      mutate(func = paste0(functionalform)) ->
      outdata




    if(stddev & sdfunctionalform == "quadratic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C )) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


    if(stddev & sdfunctionalform == "cubic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C
                        + twc * T.W.C
                        + tw2 * T.W.2
                        + tc2 * T.C.2
                        + t2w * T.2.W
                        + t2c * T.2.C
                        + w2c * W.2.C
                        + wc2 * W.C.2
                        + t3 * T.3
                        + w3 * W.3
                        + c3 * C.3)) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


  }# end If cubic


  if(functionalform == "cross"){
    coefficients %>%
      mutate(pctYieldChange = t * T.1
             + t2 * T.2
             + w * W.1
             + w2 * W.2
             + c * C.1
             + c2 * C.2
             + tw * T.W
             + tc * T.C
             + wc * W.C
             + twc * T.W.C
             + tw2 * T.W.2
             + tc2 * T.C.2
             + t2w * T.2.W
             + t2c * T.2.C
             + w2c * W.2.C
             + wc2 * W.C.2) %>%
      # select(Temp, Precip, CO2, pctYieldChange) %>%
      # distinct%>%
      mutate(func = paste0(functionalform)) ->
      outdata




    if(stddev & sdfunctionalform == "quadratic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C )) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


    if(stddev & sdfunctionalform == "cubic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C
                        + twc * T.W.C
                        + tw2 * T.W.2
                        + tc2 * T.C.2
                        + t2w * T.2.W
                        + t2c * T.2.C
                        + w2c * W.2.C
                        + wc2 * W.C.2
                        + t3 * T.3
                        + w3 * W.3
                        + c3 * C.3)) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


  }# end If cross


  if(functionalform == "pure"){
    coefficients %>%
      mutate(pctYieldChange = t * T.1
             + t2 * T.2
             + w * W.1
             + w2 * W.2
             + c * C.1
             + c2 * C.2
             + tw * T.W
             + tc * T.C
             + wc * W.C
             + t3 * T.3
             + w3 * W.3
             + c3 * C.3) %>%
      # select(Temp, Precip, CO2, pctYieldChange) %>%
      # distinct%>%
      mutate(func = paste0(functionalform)) ->
      outdata





    if(stddev & sdfunctionalform == "quadratic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C )) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }


    if(stddev & sdfunctionalform == "cubic"){


      SDcoeff %>%
        mutate_if(is.factor, as.character) %>%
        select(sigma.param, Mean) %>%
        spread(sigma.param, Mean) %>%
        repeat_add_columns(ctwvalues) %>%
        mutate(sd = abs(const + t * T.1
                        + t2 * T.2
                        + w * W.1
                        + w2 * W.2
                        + c * C.1
                        + c2 * C.2
                        + tw * T.W
                        + tc * T.C
                        + wc * W.C
                        + twc * T.W.C
                        + tw2 * T.W.2
                        + tc2 * T.C.2
                        + t2w * T.2.W
                        + t2c * T.2.C
                        + w2c * W.2.C
                        + wc2 * W.C.2
                        + t3 * T.3
                        + w3 * W.3
                        + c3 * C.3)) %>%
        select(Temp, Precip, CO2, sd) %>%
        # distinct%>%
        mutate(func = paste0(functionalform)) ->
        SDtmp


      outdata %>%
        left_join(SDtmp, by = c("func", "Temp", "Precip", "CO2")) ->
        outdata

    }



  }# end If pure


  outdata %>%
    select(Temp, Precip, CO2, pctYieldChange, sd) %>%
    distinct ->
    outdata1

  if(formaps){
    outdata %>%
      select(long, lati, latgrid, longrid, crop, irr, year, latband, Temp, Precip, CO2, pctYieldChange, sd) %>%
      distinct ->
      outdata1
  }

 return(tibble::as_tibble(outdata1))

} # end function
