### Dec. 8, 2010 simple weighted cdf
################################### changed weight of stations
weightcdf <- function(df1 = my.samp, ss = ss, SampID = "Sample.ID", xvar = "cond", nt = 25, 
        log = TRUE, np = 61) {

  df1 <- merge(df1, ss, by = SampID)
  tnames <- names(ss)[-1]

  df1 <- df1[, c(SampID, tnames, xvar)]
  #df1 <- na.omit(df1)

  numocc <- apply(df1[, tnames], 2, function(x) sum((x>0), na.rm =T))
  tnames.sav <- names(numocc)[numocc >= nt]

  cutp <- seq(from = min(df1[,xvar]), to = max(df1[,xvar]), length = np)

  cutm <- 0.5*(cutp[-1] + cutp[-np])
  df1$cutf <- cut(df1[,xvar], cutp, include.lowest = T, labels = 1:(np-1))

  wt <- 1/table(df1$cutf)
  wt[is.infinite(wt)] <- NA
  wt <- wt/sum(wt, na.rm = T)

  dftemp1 <- data.frame(cutf = names(wt), wt = as.vector(wt))

  df1 <- merge(df1, dftemp1, by = "cutf")

  tolval.e <- rep(NA, times = length(tnames.sav))

  cnew <- seq(from = min(df1[,xvar]), to = max(df1[,xvar]), length = 100)
  new.data <- data.frame(cond = cnew)

  cond.u <- sort(unique(df1[,xvar]))

  for (i in 1:length(tnames.sav)) {
      resp <- df1[, tnames.sav[i]] > 0
      df2 <- df1[resp,]
      
      tolval.e[i] <- wtd.quantile(df2[,xvar], df2$wt, normwt = TRUE, prob = 0.95)   ### R buildin ecdf function
  }
  print(tolval.e)
  if(log) XC95.all = round(10^tolval.e)
  else   XC95.all = round(tolval.e, 1)
  dftv <- data.frame(taxaname = tnames.sav, XC95.all = XC95.all,
                     N = as.vector(numocc[numocc >= nt]), stringsAsFactors = F)
  dftv <- dftv[order(dftv$XC95.all),]
  return(dftv)
 }
