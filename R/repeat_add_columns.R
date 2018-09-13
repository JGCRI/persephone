#' A function to repeat for added columns
#' Cribbed from gcamdata
#'
#' @param x dataframe of columns to be repeated for each added column of y
#' @param y dataframe of columns to add to x.
#'
#' @importFrom tibble tibble
#' @import dplyr
#' @importFrom tidyr gather spread
#' @export


repeat_add_columns <- function(x, y) {
  UNIQUE_JOIN_FIELD <- NULL           # silence package checks.

  x %>%
    mutate(UNIQUE_JOIN_FIELD = 1) %>%
    full_join(mutate(y, UNIQUE_JOIN_FIELD = 1), by = "UNIQUE_JOIN_FIELD") %>%
    select(-UNIQUE_JOIN_FIELD)
}
