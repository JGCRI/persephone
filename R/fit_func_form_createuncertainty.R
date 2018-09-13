# 11-14-2017
# Evaluates the baseline ensemble functional form via func_form_evalue_plusstddev_4.R, for several CTW grids (like IRS plotting)
# and gets std dev across sites for each.
# Then form "harsh" and "mellow" ctw-yield sets, feed into fit_func_form_ensemble3.R and get resulting coefficients that
# define our harsh and mellow response surfaces



# Can probably be smarter in terms of what coefficients are included. -> instead of repeating for each functional form,
# have conditions like If coefficient$twc exists, multiply by twc, else fill in 0.



function(coeff, ctwmesh = 101,
         baseCO2 = 360, co2lo = 330, co2hi = 900,
         baseTemp = 0, templo = -1, temphi = 1,
         basePrecip = 1,preciplo = 0.5, preciphi = 1.5,
         functionalform, stddevMult = 0.5,
         nme){

  # build out large variety of ctw and corresponding yields and standard deviations, to use to create optimistic and pessimisstic
  # ctwY tables to fit coefficients to.
        ctwY <- tibble::tibble()


        #############################################################################################################
        # Fixed CO2 = baseline
        temp <- seq(-1,8, length = ctwmesh)
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- baseCO2
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,  baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        #############################################################################################################
        # Fixed CO2 = co2lo
        temp <- seq(-1,8, length = ctwmesh)
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- co2lo
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,  baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        #############################################################################################################
        # Fixed CO2 = co2hi
        temp <- seq(-1,8, length = ctwmesh)
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- co2hi
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,  baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        #############################################################################################################
        # Fixed CO2 = baseCO2 + (co2hi - baseCO2) / 3
        temp <- seq(-1,8, length = ctwmesh)
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- baseCO2 + (co2hi - baseCO2) / 3
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,  baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        #############################################################################################################
        # Fixed CO2 = baseCO2 + 2 * (co2hi - baseCO2) / 3
        temp <- seq(-1,8, length = ctwmesh)
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- baseCO2 + 2 * (co2hi - baseCO2) / 3
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Temp = baseTemp
        temp <- baseTemp
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- seq(330,900, length = ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY



        ############################################################################################################
        # Fixed Temp = templo
        temp <- templo
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- seq(330,900, length = ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Temp = temphi
        temp <- temphi
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- seq(330,900, length = ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY



        ############################################################################################################
        # Fixed Temp = baseTemp + (temphi - baseTemp) /3
        temp <- baseTemp + (temphi - baseTemp) /3
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- seq(330,900, length = ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Temp = baseTemp + 2 * (temphi - baseTemp) /3
        temp <- baseTemp + 2 * (temphi - baseTemp) /3
        precip <- seq(0.5, 1.5, length = ctwmesh)
        co2 <- seq(330,900, length = ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Precip = baseline
        temp <- seq(-1,8, length = ctwmesh)
        precip <- basePrecip
        co2 <- seq(330,900,length= ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Precip = preciphi
        temp <- seq(-1,8, length = ctwmesh)
        precip <- preciphi
        co2 <- seq(330,900,length= ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Precip = preciplo
        temp <- seq(-1,8, length = ctwmesh)
        precip <- preciplo
        co2 <- seq(330,900,length= ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Precip = basePrecip + (preciplo - basePrecip) / 2
        temp <- seq(-1,8, length = ctwmesh)
        precip <- basePrecip + (preciplo - basePrecip) / 2
        co2 <- seq(330,900,length= ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,   baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


        ############################################################################################################
        # Fixed Precip = basePrecip + (preciphi - basePrecip) / 2
        temp <- seq(-1,8, length = ctwmesh)
        precip <- basePrecip + (preciphi - basePrecip) / 2
        co2 <- seq(330,900,length= ctwmesh)
        CTWgrid <- expand.grid(temp, precip, co2)
        colnames(CTWgrid) <- c("Temp", "Precip", "CO2")
        CTWgrid %>%
          tibble::as_tibble() %>%
          mutate(Temp = (Temp - standardize.CTW.vals$T.mean)/standardize.CTW.vals$T.sd,
                 Precip = (Precip - standardize.CTW.vals$W.mean) / standardize.CTW.vals$W.sd,
                 CO2 = (CO2 - standardize.CTW.vals$C.mean) / standardize.CTW.vals$C.sd) ->
          CTWgrid

        # evaluate the emulator on the CTW grid
        coeff %>%
          func_form_evaluate(coeff = ., ctwvalues = CTWgrid,  baseC = baseCO2, baseT = baseTemp,
                             baseP = basePrecip, functionalform = functionalform, stddev = TRUE) %>%
          bind_rows(ctwY) ->
          ctwY


  ############################################################################################################
  ############################################################################################################
  # Form the Harsh/pessimistic CTW Y tablea and fit
        ctwY %>%
          # mutate(yield = if_else(pctYieldChange < 0, pctYieldChange - stddevMult * sd, pctYieldChange + stddevMult * sd)) %>%
          mutate(yield = pctYieldChange + sign(pctYieldChange) * stddevMult * sd) %>%
          select(func, site, Temp, Precip, CO2, yield) %>%
          fit_func_form_general() %>%
          mutate(site = "enhanced") ->
          ctwY_pess




  # Form the mellow/optimisitc CTW Y table
        ctwY %>%
          # mutate(yield = if_else(pctYieldChange < 0, pctYieldChange + stddevMult * sd, pctYieldChange - stddevMult * sd)) %>%
          mutate(yield = pctYieldChange - sign(pctYieldChange) * stddevMult * sd) %>%
          select(func, site, Temp, Precip, CO2, yield) %>%
          fit_func_form_general() %>%
          mutate(site = "moderated") ->
          ctwY_opt

        return(list(ctwY_pess, ctwY_opt))


} # end function
