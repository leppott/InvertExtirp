#' A function to determine the curve shape of a binomial probability
#'
#' Adapted from Lester Yuan's code.  Used in taxon.response.sort().
#'
#' @param mnr mean predicted probability
#' @param ubnd upper boundary of mean predicted probability
#' @param lbnd lower boundary of mean predicted probability
#' @return returns one of "Unimodal", "Increasing", "Decreasing", "Concave up", or NA.
#' @keywords logistic regression, quantiles, xc95, hc05, cdf, gam, taxon response
#' @examples
#' mean.resp <- 1:1000
#' up.bound <- mean.resp+1
#' low.bound <- mean.resp-1
#' curve.shape(mean.resp, up.bound,low.bound)
#' @export
curve.shape <- function(mnr, ubnd, lbnd) {##FUNCTION.curve.shape.START
  ######
  # mnr, ubnd, and lbnd are mean, upper, and lower prediction confidence interval
  #####
    # Find the maximum and minimum predicted mean probabilities
    lmax <- max(mnr)
    lmin <- min(mnr)
    # Find index locations for these probabilities
    imax <- match(lmax, mnr)
    imin <- match(lmin, mnr)
    x.out <- F
    y.out <- F

    # Compare mean predicted probability to the left of maximum point
    # with upper confidence bound. Store a T in x.out if
    # any point in the mean response deviates from the
    # upper confidence limit
    if (imax > 1) {
      x.out <- sum(lmax == pmax(lmax, ubnd[1:(imax-1)])) > 0
    }

    # Store a T in y.out if any point in the mean probability
    # to the right of the maximum point deviates from the upper
    # confidence limit
    if (imax < length(ubnd)) {
      y.out <- sum(lmax == pmax(lmax,
      ubnd[(imax+1):length(ubnd)])) > 0
    }

    # Perform same set of tests for lower confidence limit
      a.out <- F
      b.out <- F
      if (imin > 1) {
        a.out <- sum(lmin == pmin(lmin, lbnd[1:(imin-1)])) > 0
      }
      if (imin < length(lbnd)) {
        b.out <- sum(lmin == pmin(lmin,
      lbnd[(imin+1):length(lbnd)])) > 0
      }
    # The information on where the mean curve deviates from the
    # confidence limits tells us its curve shape...


      if (x.out & y.out) {
        return("Unimodal")
      } else if (a.out & b.out) {
        return("Concave up")
      } else if (x.out | b.out) {
        return("Increasing")
      } else if (y.out | a.out) {
        return("Decreasing")
      } else {
        return(NA)
      }
      # else # if (! (x.out | y.out | a.out | b.out)) {
      #   return(NA)
      # #}
}##FUNCTION.curve.shape.END
