#' A function for leave one out cross validation on a predefined partition of data
#' @param ensSiteData A list of sites for ensemble under consideration
#' @param testnmes The names of the CTW tests
#' @param ctwtests The actual ctw tests in full
#' @param functionalform The name functional form for the mean of the likelihood
#' @param sdfunctionalform The name of the functional form for the standard deviation of the likelihood
#' @param ensemblename The name of the ensemble under consideration
#' @param testpartition The names of the CTW tests that will be excluded from training to be used for testing.
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @importFrom tidyr gather spread
#' @importFrom rethinking map precis
#' @export



cross_validate_LOO <- function(ensSiteData, testnmes, ctwtests, individSiteDataPath, functionalform, sdfunctionalform, ensemblename, testpartition){

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
  fulldata


fulldata %>%
  filter(!(Test %in% testpartition)) ->
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

fulldata %>%
  filter(Test %in% testpartition) ->
  fulltestdata


fulltestdata %>%
  select(Test, Yield) %>%
  group_by(Test) %>%
  summarise(meanY = mean(Yield),
            highY = mean(Yield) + sd(Yield),
            lowY  = mean(Yield) - sd(Yield)) %>%
  ungroup ->
  testdata

ctwtests %>%
  filter(Test %in% testpartition) ->
  testCTW


precis(m1, digits = 12)@output %>%
  mutate(param = row.names(.)) %>%
  select(param, Mean) %>%
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


func_form_evaluate(coeff = mu.coeff, SDcoeff = sigma.coeff, ctwvalues = testCTW,  baseC = 360, baseT = 0, baseP = 1,
                   functionalform = functionalform, sdfunctionalform = sdfunctionalform,
                   ensembleName = ensemblename,stddev = TRUE) %>%
  #select(-func) %>%
  summarise(emuMeanY = pctYieldChange,
            emuHighY = pctYieldChange + sd,
            emuLowY  = pctYieldChange - sd) %>%
  bind_cols(tibble::tibble(Test = testCTW$Test)) %>%
  left_join(testdata, by = "Test") %>%
  summarise(mean_diff = emuMeanY - meanY,
            high_diff = emuHighY - highY,
            low_diff = emuLowY - lowY) %>%
  ungroup ->
  compareY


tibble::tibble(Ensemble = ensemblename,
               form = paste0(func, "_", sdfunc)) %>%
  bind_cols(compareY) ->
  output
return(output)



} # end function
