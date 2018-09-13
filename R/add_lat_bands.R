#' A function to categorize latitudes
#'
#' @param data An input data frame with a column, latcol, of latitude data to be categorized into tropic, mid, and high
#' @param tropic_cut The boundary between tropic and mid latitudes, defaults to 30
#' @param mid_cut The boundary between mid and high latitudes, defaults to 70
#'
#' @import dplyr
#' @export

add_lat_bands <- function(data, tropic_cut = 30, mid_cut = 70){
  data %>%
    mutate(latBand = if_else(-tropic_cut < latcol & latcol < tropic_cut, "tropic" ,
                             if_else(tropic_cut <= latcol & latcol <= mid_cut, "mid",
                                     if_else(-mid_cut <= latcol & latcol <= -tropic_cut, "mid",
                                             if_else(latcol > mid_cut, "high",
                                                     if_else(latcol < -mid_cut, "high", "err"))))))
}
