#' @title Calculate HC05
#'
#' @description The HC05 is calculated on each specified column of a specified data frame.
#'
#' @details The HC05 is the 5th percentile hazardous concentration of a set of
#' related 95th percentile extirpation concentration values (XC95). The HC05 is
#' the estimated value of the stressor at which 5% of species will be extirpated.
#'
#' @param df_data Data frame of taxa data with XC95 values in one or more columns.
#' @param col_XC95  Column names for XC95 values to calculate the HC05.
#'
#' @return A vector of values with names corresponding to the provided column names.
#'
#' @examples
#' #' \dontrun{
#' # Calculate XC95 values using fish.wt.cdf example
#' # data
#' data(dta.do)
#' data(ss.sites)
#' # run function (~20 seconds)
#' df_tv_do <- fish.wt.cdf(datafile = dta.do, ss = ss.sites, plot = T, dogam = T
#'                        , SampleID = "Station_Date", tag = "wt", sortvect = NULL
#'                        , np = 61, nt = 25, addtrend = T
#'                        , wd = getwd(), groups = c("BigHUC","ECOREGL3","WS_AREA")
#'                        , xvar = "cond")
#' # Calculate HC05 values
#' hc05_tv_do <- hc05(df_tv_do, c("XC95.cdf", "XC95.gam"))
#' print(hc05_tv_do)
#' }
#'
#' @export
hc05 <- function(df_data, col_XC95){##FUNCTION.START
  #
  hc05_calc <- apply(df_data[, c(col_XC95)], 2
                    , quantile, probs=0.05, na.rm=TRUE, type=6)
  #
  return(hc05_calc)
  #
}##FUNCTION.END
