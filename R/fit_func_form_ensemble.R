#' A function to fit the Bayesian model and write the resulting rethinking model,
#' point estimates and the standard deviation of the coefficients for mu,
#' point estimates and the standard deviation of the coefficients for sigma,
#' and a plot of the above data.
#'
#' @param ensSiteData A list of sites for ensemble under consideration
#' @param testnmes The names of the CTW tests
#' @param ctwtests The actual ctw tests in full
#' @param individuSiteDataPath The location of the raw C3MP yield data, defaults to "data/c3mpdata/csvs/"
#' @param coefpath The location to write the outputs, defaults to "data/emulatorcoefficients"
#' @param functionalform The name functional form for the mean of the likelihood
#' @param sdfunctionalform The name of the functional form for the standard deviation of the likelihood
#' @param figurepath The location to write the plots of the coefficients + uncertainty, defaults to "figures/bayeseval"
#' @param ensemblename The name of the ensemble under consideration
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @import rstan
#' @import rethinking
#' @importFrom tidyr gather spread
#' @export



fit_func_form_ensemble <- function(ensSiteData, testnmes, ctwtests,
                                   individSiteDataPath = "data/c3mpdata/csvs/",
                                   coefpath = "data/emulatorcoefficients",
                                   functionalform, sdfunctionalform,
                                   figurepath = "figures/bayeseval",
                                   ensemblename){

### form the vector of ensemble member sites
simSites <- ensSiteData$site2

### read in and bind data for simulation sites to form the ensemble.
data <- tibble::tibble()
for (i in 1:length(simSites)){
    read.csv(file=paste0(individSiteDataPath, "/", simSites[i], ".csv"),head = FALSE, sep = ",")[1:30,1:99] %>%
    mutate(site = paste0(simSites[i])) %>%
    bind_rows(data) ->
    data
}


colnames(data) <- c(testnmes, "site")

data %>%
  gather(Test, Yield, -site) ->
  data2

# add the ctw tests and standardize the variables
ctwtests %>%
  left_join(data2, by = "Test") %>%
  left_join(select(siteData, site2, baseY), by = c("site" = "site2")) %>%
  group_by(Test, Temp, Precip, CO2, site) %>%
  summarize(Yield = mean(Yield),
            baseY = unique(baseY)) %>%
  ungroup %>%
  mutate(Yield = 100 * (Yield - baseY)/baseY ) %>%
  as.data.frame() ->
  x

x %>%
  convert_ctw(.) %>%
  as.data.frame() ->
  data3


if(functionalform == "c3mp" & sdfunctionalform == "quadratic"){

  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      twc = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05, # *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
      + t2 * T.2
      + w * W.1
      + w2 * W.2
      + c * C.1
      + c2 * C.2
      + tw * T.W
      + tc * T.C
      + wc * W.C
      + twc * T.W.C,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )
  }# end If c3mp


if(functionalform == "c3mp" & sdfunctionalform == "cubic"){

  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      twc = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05, # *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0,
                      sigma.twc = 0,
                      sigma.tw2 = 0,
                      sigma.tc2 = 0,
                      sigma.t2w = 0,
                      sigma.t2c = 0,
                      sigma.w2c = 0,
                      sigma.wc2 = 0,
                      sigma.t3 = 0,
                      sigma.w3 = 0,
                      sigma.c3 = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
      + t2 * T.2
      + w * W.1
      + w2 * W.2
      + c * C.1
      + c2 * C.2
      + tw * T.W
      + tc * T.C
      + wc * W.C
      + twc * T.W.C,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C
      + sigma.twc * T.W.C
      + sigma.tw2 * T.W.2
      + sigma.tc2 * T.C.2
      + sigma.t2w * T.2.W
      + sigma.t2c * T.2.C
      + sigma.w2c * W.2.C
      + sigma.wc2 * W.C.2
      + sigma.t3 * T.3
      + sigma.w3 * W.3
      + sigma.c3 * C.3,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )
}# end If c3mp





if(functionalform == "quadratic" & sdfunctionalform == "quadratic"){

  start_list <- list( #t = mean(data3$Yield),
                      t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05,#  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
      + t2 * T.2
      + w * W.1
      + w2 * W.2
      + c * C.1
      + c2 * C.2
      + tw * T.W
      + tc * T.C
      + wc * W.C,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 50000)
  )


}# end If quadratic


if(functionalform == "quadratic" & sdfunctionalform == "cubic"){

  start_list <- list( #t = mean(data3$Yield),
    t = 0,
    t2 = 0,
    w = 0,
    w2 = 0,
    c = 0,
    c2 = 0,
    tw = 0,
    tc = 0,
    wc = 0,
    # again, for solvers to work and avoid dividing by 0, we have to include a constant term
    # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
    # with high confidence at baseline.
    sigma.const  = 0.05,#  sd(data3$Yield),
    sigma.t = 0, #sd(data3$Yield) ,
    sigma.t2 = 0,
    sigma.w = 0,
    sigma.w2 = 0,
    sigma.c = 0,
    sigma.c2 = 0,
    sigma.tw = 0,
    sigma.tc = 0,
    sigma.wc = 0,
    sigma.twc = 0,
    sigma.tw2 = 0,
    sigma.tc2 = 0,
    sigma.t2w = 0,
    sigma.t2c = 0,
    sigma.w2c = 0,
    sigma.wc2 = 0,
    sigma.t3 = 0,
    sigma.w3 = 0,
    sigma.c3 = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
      + t2 * T.2
      + w * W.1
      + w2 * W.2
      + c * C.1
      + c2 * C.2
      + tw * T.W
      + tc * T.C
      + wc * W.C,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C
      + sigma.twc * T.W.C
      + sigma.tw2 * T.W.2
      + sigma.tc2 * T.C.2
      + sigma.t2w * T.2.W
      + sigma.t2c * T.2.C
      + sigma.w2c * W.2.C
      + sigma.wc2 * W.C.2
      + sigma.t3 * T.3
      + sigma.w3 * W.3
      + sigma.c3 * C.3,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 50000)
  )


}# end If quadratic





if(functionalform == "cubic" & sdfunctionalform == "quadratic"){
  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      twc = 0,
                      tw2 = 0,
                      tc2 = 0,
                      t2w = 0,
                      t2c = 0,
                      w2c = 0,
                      wc2 = 0,
                      t3 = 0,
                      w3 = 0,
                      c3 = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05,# *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0)

  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
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
      + c3 * C.3,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )


}# end If cubic


if(functionalform == "cubic" & sdfunctionalform == "cubic"){
  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      twc = 0,
                      tw2 = 0,
                      tc2 = 0,
                      t2w = 0,
                      t2c = 0,
                      w2c = 0,
                      wc2 = 0,
                      t3 = 0,
                      w3 = 0,
                      c3 = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05,# *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0,
                      sigma.twc = 0,
                      sigma.tw2 = 0,
                      sigma.tc2 = 0,
                      sigma.t2w = 0,
                      sigma.t2c = 0,
                      sigma.w2c = 0,
                      sigma.wc2 = 0,
                      sigma.t3 = 0,
                      sigma.w3 = 0,
                      sigma.c3 = 0
                      )

  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
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
      + c3 * C.3,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C
      + sigma.twc * T.W.C
      + sigma.tw2 * T.W.2
      + sigma.tc2 * T.C.2
      + sigma.t2w * T.2.W
      + sigma.t2c * T.2.C
      + sigma.w2c * W.2.C
      + sigma.wc2 * W.C.2
      + sigma.t3 * T.3
      + sigma.w3 * W.3
      + sigma.c3 * C.3,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )


}# end If cubic





if(functionalform == "cross" & sdfunctionalform == "quadratic"){
  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      twc = 0,
                      tw2 = 0,
                      tc2 = 0,
                      t2w = 0,
                      t2c = 0,
                      w2c = 0,
                      wc2 = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05,# *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
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
      + wc2 * W.C.2,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )


}# end If cross


if(functionalform == "cross" & sdfunctionalform == "cubic"){
  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      twc = 0,
                      tw2 = 0,
                      tc2 = 0,
                      t2w = 0,
                      t2c = 0,
                      w2c = 0,
                      wc2 = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05,# *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0,
                      sigma.twc = 0,
                      sigma.tw2 = 0,
                      sigma.tc2 = 0,
                      sigma.t2w = 0,
                      sigma.t2c = 0,
                      sigma.w2c = 0,
                      sigma.wc2 = 0,
                      sigma.t3 = 0,
                      sigma.w3 = 0,
                      sigma.c3 = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
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
      + wc2 * W.C.2,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C
      + sigma.twc * T.W.C
      + sigma.tw2 * T.W.2
      + sigma.tc2 * T.C.2
      + sigma.t2w * T.2.W
      + sigma.t2c * T.2.C
      + sigma.w2c * W.2.C
      + sigma.wc2 * W.C.2
      + sigma.t3 * T.3
      + sigma.w3 * W.3
      + sigma.c3 * C.3,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )


}# end If cross





if(functionalform == "pure" & sdfunctionalform == "quadratic"){
  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      t3 = 0,
                      w3 = 0,
                      c3 = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05, # *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
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
      + c3 * C.3,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )

}# end If pure


if(functionalform == "pure" & sdfunctionalform == "cubic"){
  start_list <- list( t = mean(data3$Yield),
                      #t = 0,
                      t2 = 0,
                      w = 0,
                      w2 = 0,
                      c = 0,
                      c2 = 0,
                      tw = 0,
                      tc = 0,
                      wc = 0,
                      t3 = 0,
                      w3 = 0,
                      c3 = 0,
                      # again, for solvers to work and avoid dividing by 0, we have to include a constant term
                      # for sigma. We use our initial guess and a prior to constrain it to be very close to 0
                      # with high confidence at baseline.
                      sigma.const  = 0.05, # *  sd(data3$Yield),
                      sigma.t = 0, #sd(data3$Yield) ,
                      sigma.t2 = 0,
                      sigma.w = 0,
                      sigma.w2 = 0,
                      sigma.c = 0,
                      sigma.c2 = 0,
                      sigma.tw = 0,
                      sigma.tc = 0,
                      sigma.wc = 0,
                      sigma.twc = 0,
                      sigma.tw2 = 0,
                      sigma.tc2 = 0,
                      sigma.t2w = 0,
                      sigma.t2c = 0,
                      sigma.w2c = 0,
                      sigma.wc2 = 0,
                      sigma.t3 = 0,
                      sigma.w3 = 0,
                      sigma.c3 = 0)


  m1 <-  map(
    alist(
      Yield ~ dnorm(mu, abs(sig)),
      mu <-   t * T.1
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
      + c3 * C.3,
      sig <-  sigma.const + sigma.t * T.1
      + sigma.t2 * T.2
      + sigma.w * W.1
      + sigma.w2 * W.2
      + sigma.c * C.1
      + sigma.c2 * C.2
      + sigma.tw * T.W
      + sigma.tc * T.C
      + sigma.wc * W.C
      + sigma.twc * T.W.C
      + sigma.tw2 * T.W.2
      + sigma.tc2 * T.C.2
      + sigma.t2w * T.2.W
      + sigma.t2c * T.2.C
      + sigma.w2c * W.2.C
      + sigma.wc2 * W.C.2
      + sigma.t3 * T.3
      + sigma.w3 * W.3
      + sigma.c3 * C.3,
      sigma.const ~ dnorm(0, 0.01)
    ),
    data = data3, start = start_list, control = list(maxit = 10000)
  )




}# end If pure





####################
# output


precis(m1, digits = 12)@output %>%
  mutate(param = row.names(.)) %>%
  select(param, Mean, StdDev) %>%
  separate(param, c("mu.param", "sigma.param")) ->
  coeff

coeff %>%
  filter(mu.param != "sigma") %>%
  select(-sigma.param) ->
  mu.coeff

coeff %>%
  filter(mu.param == "sigma") %>%
  select(-mu.param) ->
  sigma.coeff


write.csv(mu.coeff, paste0(coefpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_", ensemblename,"_mu_coeff.csv"), row.names = FALSE)

write.csv(sigma.coeff, paste0(coefpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_", ensemblename,"_sigma_coeff.csv"), row.names = FALSE)

save(m1, file = paste0(coefpath, "/", functionalform, "_mu_", sdfunctionalform, "_sigma_", ensemblename,"_model.rda"))



} # end function
