#' A function to compute area under a curve.
#'
#' Used in other functions in this library (taxon.response(), taxon.response.sort(), and tolerance()).
#' Different versions for response and tolerance.
#'
#' @param xrange y
#' @param mod x
#' @param dense.N the number of areas to cut into in the calculation of area under the curve; default = 201
#' @param taus determine the output the percentile of env variable; default = c(0,95,100) for response.
#' @param auc.version calculation version, "response" or "tolerance"; default = "response".
#' @return returns area under the specified curve.
#' @keywords area under curve,
#' @examples
#' xrange <- 1
#' lrm1 <- 2
#'
#' # response version
#'
#' dense.N <- 201
#' taus=c(0,95,100)
#' #' lrm1.95f <- auc.response(xrange, lrm1, dense.N, taus, "response")
#' # tolerance version
#'
#' taus = c(50,95)
#' lrm195f <- auc(xrange, lrm1, dense.N, taus, "tolerance")
#' @export
auc <- function(xrange, mod, dense.N, taus, auc.version = "response") {##FUNCTION.auc.START
  #
  x <- seq(min(xrange), max(xrange), length = dense.N - 1)
  s.area <-rep(NA, dense.N - 1)
  y <- predict(mod, newdata = data.frame(dose = x), type = "response")
  for (index in 1:dense.N-1) {
    s.area[index] <- (y[index] + y[index + 1])/2*(x[index+1]-x[index])
  }
  tsum <- sum(s.area, na.rm=T)
  jj=1
  csum = sum(s.area[1:jj])
  #
  auc.version <- tolower(auc.version)
  if (auc.version=="response") {##IF.version.START
    #
    while(csum < taus[2]/100*tsum) {
      jj = jj + 1
      csum <- sum(s.area[1:jj], na.rm=T)
    }
    #
    xc95 <- (x[jj]+ x[jj-1])/2
    yc95 <- (y[jj]+ y[jj-1])/2
    #
    return(c(xc95,yc95))
    #
  } else if (auc.version=="tolerance") {
    #
    xc95 <- yc95 <- rep(NA, length(taus))
    for(ii in 1:length(taus)) {
      while(csum < taus[ii]/100*tsum) {
        jj = jj + 1
        csum <- sum(s.area[1:jj], na.rm=T)
      }
      if(jj == 1) {
        xc95[ii] <- x[jj]
        yc95[ii] <- y[jj]
      } else {
        xc95[ii] <- (x[jj]+ x[jj-1])/2
        yc95[ii] <- (y[jj]+ y[jj-1])/2
      }
    }
    return(c(xc95))
    #
  } else {
    msg <- "Incorrect input; 'auc.version' has to be 'response' or 'tolerance'."
    stop(msg)
  }
  #
}##FUNCTION.auc.END

