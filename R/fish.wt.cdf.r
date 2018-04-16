  ##~~~~~~~~~~~~~~  4)  fish condutivity weighted method for appalanchian region
  # Chem.eco3 ; fish
##
##  weighted.cdf function to calculate XC95 values
##  datafile  environmental data
##  ss, species crosstabed data
##  plot, a boolean to choose if plot cdf and gam plots
##  dogam, a booleen to choose if a gam fit is calculated
##  sortvect: to provide a vect of species list so plots will be sorted according to the list
##  nt: minimum number of occurence
##  addtrend: if a trend should be added ( = ">" etc) in the output
##  np: number of bins
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#' A function to calculate a weighted CDF for XC95 values.
#'
#' Weighted cdf function to calculate XC95 values
#'#'
#' @param datafile environmental data, default = "datafile" from global environment.
#' @param ss Species crosstabed data; default = "ss" from global environment.
#' @param plot A boolean to choose if plot cdf and gam plots; default = T.
#' @param dogam A booleen to choose if a gam fit is calculated; default = T.
#' @param SampleID Site/sample id column; default = "Station_Date"
#' @param tag Default = "".
#' @param sortvect to provide a vector of species list so plots will be sorted according to the list; default = NULL.
#' @param np Number of bins; default = 61.
#' @param nt Minimum number of occurence; default = 25.
#' @param addtrend A booleen if a trend should be added ( = ">" etc) in the output (T or F); default = T.
#' @param wd Working directory for saving files.
#' @param groups column names in datafile used for grouping the data; HUC (BigHUC), Ecoregion (ECOREGL3), and Watershed Area (WS_Area)
#' @param xvar variable on which to base calculations; default  = "cond"
#' @return A dataframe of XC95 values and TIFF files of gam and/or cdf plots in the subdirectory "Results" of the directory specified by "wd".
#'
#' @examples
#' # data
#' data(dta.do)
#' data(ss.sites)
#' # run function (~20 seconds)
#' dftv.do <- fish.wt.cdf(datafile = dta.do, ss = ss.sites, plot = T, dogam = T
#'                       , SampleID = "Station_Date", tag = "wt", sortvect = NULL
#'                       , np = 61, nt = 25, addtrend = T
#'                       , wd = getwd(), groups = c("BigHUC","ECOREGL3","WS_AREA")
#'                       , xvar = "cond")
#' View(dftv.do)
#~~~~~~~~~~~~~~~~~~~
# QC
# data(dta.do)
# data(ss.sites)
# datafile = dta.do
# ss = ss.sites
# plot = T
# dogam = T
# SampleID = "Station_Date"
# tag = "wt"
# sortvect = NULL
# np = 61
# nt = 25
# addtrend = T
# wd = getwd()
# groups = c("BigHUC","ECOREGL3","WS_AREA")
# xvar = "cond"
# index <- 1
#~~~~~~~~~~~~~~~~~~~~~
#' @export
fish.wt.cdf <- function(datafile = datafile, ss = ss, plot = T, dogam = T,
              SampleID = "Station_Date", tag = "", sortvect = NULL, np = 61, nt = 25, addtrend = T,
              wd = getwd(), groups = c("BigHUC","ECOREGL3","WS_AREA"), xvar = "cond") {##FUNCTION.fish.wt.cdf.START
  #
  # 20170424, change df1$cond to df1[,xvar]
  #
  # 20170424, add directory check
  region <- "Results"
  dir.check.add(wd,region)
  dir.check.add(file.path(wd,region),"cdf")
  dir.check.add(file.path(wd,region),"gam")
  #wd <- file.path(wd,region)
  #
  #my.ss <- merge(datafile[c(SampleID, "HUC", "cond", "BigHUC","ECOREGL3","WS_AREA")], ss); dim(my.ss)   # merge env data with species data
  # 20170418
  my.ss <- merge(datafile,ss); dim(my.ss)

  # 20170418, size of datafile
  tcol.start <- ncol(datafile)+1

  tnames.count <- apply(my.ss[tcol.start:ncol(my.ss)] > 0, 2, sum)            # preliminary count of number of occurences
  tnames.sav <- tnames.count[tnames.count >= nt ]; length(tnames.sav)
  if( is.null(sortvect)) {##IF.sortvect.START
    mod.tnames <- names(tnames.sav)
  } else mod.tnames <- sortvect
  ##IF.sortvect.END

  tolval.cdf <- rep(NA, times = length(mod.tnames))  # xc95 from cdf
  tolval.gam<- rep(NA, times = length(mod.tnames))   # xc95 from gam
  total.n <- rep(NA, times = length(mod.tnames))    # Xc95 for number of occurence sites
  total.N <- rep(NA, times = length(mod.tnames))    # xc95 for total number of sites
  trend <- rep(NA, times = length(mod.tnames))      ##  > or ~ or <
  wtshd.df <- rep(NA, times = length(mod.tnames))
  eco3.df <- rep(NA,times = length(mod.tnames))
  huc.df <- rep(NA,times = length(mod.tnames))

#  row.names(wtshd.df) <- mod.tnames
#  if(plotcdf) {        # device setup



#    pdf(file = paste(wd, "/Results/",tag,"cdf.fish.wtshed.pdf", sep = ""), width = 9, height = 6, pointsize = 11)
#    par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#  } else if (plotgam) {
#    pdf(file = paste(wd, "/Results/",tag,"gam.fish.wtshed.pdf", sep = ""), width = 9, height = 6, pointsize = 11)
#    par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#  }

  # parameter-ize (EWL, 20170424)
  groups1 <- groups[1] # "BigHUC"
  groups2 <- groups[2] # "ECOREGL3"
  groups3 <- groups[3] # "WS_AREA"


  for (index in 1:length(mod.tnames)) {##IF.index.START
    ## for mod.taxa list
    imatch <- match(mod.tnames[index], names(my.ss)); imatch      # find colum index
    #HUCs <- unique(my.ss[my.ss[,imatch] > 0 ,"BigHUC"])   # find all Site HUCs
    HUCs <- unique(my.ss[my.ss[,imatch] > 0 ,groups1])   # find all Site HUCs
    huc.list <- paste(sort(HUCs), collapse ="_")
    #eco3 <- paste(unique(my.ss[my.ss[,imatch]>0,"ECOREGL3"]), collapse="_")
    eco3 <- paste(unique(my.ss[my.ss[,imatch]>0,groups2]), collapse="_")
    #wtshed <- paste(range(my.ss[my.ss[,imatch]>0,"WS_AREA"],na.rm =T), collapse="_")
    wtshed <- paste(range(my.ss[my.ss[,imatch]>0,groups3],na.rm =T), collapse="_")
    df1 <- subset(my.ss, BigHUC %in% HUCs); dim(df1)
    #print(paste(nrow(my.ss), nrow(df1)))

    # cnew <- seq(from = min(df1$cond), to = max(df1$cond), length = 100)      # for plot
    # new.data <- data.frame(cond = cnew)
    # cond.u <- sort(unique(df1$cond))
    cnew <- seq(from = min(df1[,xvar]), to = max(df1[,xvar]), length = 100)      # for plot
    new.data <- data.frame(cond = cnew)
    cond.u <- sort(unique(df1[,xvar]))

    #cutp <- seq(from = min(df1$cond), to = max(df1$cond), length = np)  # cut entire gradient into np bins
    cutp <- seq(from = min(df1[,xvar]), to = max(df1[,xvar]), length = np)  # cut entire gradient into np bins
    cutm <- 0.5*(cutp[-1] + cutp[-np])                   # find middle point
    #df1$cutf <- cut(df1$cond, cutp, include.lowest = T)   # define bins
    df1$cutf <- cut(df1[,xvar], cutp, include.lowest = T)   # define bins


    wt <- 1/table(df1$cutf)                 # weight for each bin
    wt[is.infinite(wt)] <- NA
    wt <- wt/sum(wt, na.rm = T)

    dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
    df1 <- merge(df1, dftemp, by = "cutf"); dim(df1); names(df1)  # final data

    #  ECDF of observed conductivity
    resp <- df1[, mod.tnames[index]] > 0              # response variable
    df2 <- df1[resp,]

    #tolval.cdf[index] <- Hmisc::wtd.quantile(df2$cond, df2$wt, normwt = TRUE, prob = 0.95) # xc95 calculation
    tolval.cdf[index] <- Hmisc::wtd.quantile(df2[,xvar], df2$wt, normwt = TRUE, prob = 0.95) # xc95 calculation
    total.n[index] <- nrow(df2)
    total.N[index] <- nrow(df1)
    eco3.df[index] <- eco3
    wtshd.df[index] <- wtshed
    huc.df[index] <- huc.list
    #wtshd.df[index,] <- c(mywtshed[1], mywtshed[2], huc.list)
    if (dogam) {##IF.dogam.START
      ### only if dogam is selected, calculate gam based xc95 and plot
      mod <- mgcv::gam(resp ~ s(cond, k = 3), data = df1, family = "binomial")
      predresp <- predict(mod, new.data, type = "link", se.fit = T)
      # Compute upper and lower 90% confidence limits
      up.bound.link <- predresp$fit + qnorm(0.95) * predresp$se.fit
      low.bound.link <- predresp$fit - qnorm(0.95) * predresp$se.fit

      # Convert from logit transformed values to probability.
      up.bound <- exp(up.bound.link)/(1 + exp(up.bound.link))
      low.bound <- exp(low.bound.link)/(1 + exp(low.bound.link))

      predresp.vline <- predict(mod, newdata = data.frame(cond = max(df1$cond, na.rm=T)), type = "link", se.fit = T)
      low500.lk <- predresp.vline$fit - qnorm(0.95) * predresp.vline$se.fit
      low500.resp <- exp(low500.lk)/(1 + exp(low500.lk))
      mean500.lk <- predresp.vline$fit
      mean500.resp <- exp(mean500.lk)/(1 + exp(mean500.lk))
      mn <- plogis(predresp$fit)
      crit <- 0.05 * max(mn)           # arbitrarily define a criterion to assign trend
      if(mean500.resp <= crit) {##IF.mean500.resp.START
        trend[index] <- "="
      }  else  if(low500.resp <= crit & mean500.resp> crit ) {
        trend[index] <- "~"
      }  else {
        trend[index] <- ">"
      }##IF.mean500.resp.END


      bararea <- (cnew[2] - cnew[1])* mn         ### determine gam xc95 values
      tot <- sum(bararea)
      ii <- 1
      while(sum(bararea[1:ii]) < 0.95 * tot) ii <- ii + 1
      tolval.gam[index] <- cnew[ii]
    }##IF.dogam.END

    if(plot) {##IF.plot.START
       if(index%%6==1 ) {##IF.index.START
       pgs <- (index-1)%/%6 + 1
       tiff(file = paste(wd,"/Results/cdf/", pgs, ".taxon.cdf.tiff",sep=""),
           width = 650, height = 450, pointsize = 13)
       par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))

       tiff(file = paste(wd,"/Results/gam/", pgs, ".taxon.gam.tiff",sep=""),
           width = 650, height = 450, pointsize = 13)
        par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
        }##IF.index.END
     dev.set(2)
     Hmisc::Ecdf(df2[,"cond"], weights = df2$wt, ylim = c(0,1), col = "blue", pch = 1, axes = F,
           main = bquote(italic(.(mod.tnames[index]))), ylab = "Proportion of occurrence ",
           xlab = expression(paste("Conductivity ( ", mu, "S/cm)")))
      abline(h = 0.95, col = "red", lty = 2)
      abline(v = tolval.cdf[index], col = "red", lty = 2)
      box(bty = "l")
      min.pow <- floor(min(df1$cond)); max.pow <- ceiling(max(df1$cond))
      at0 <- min.pow:max.pow              # add ticks at log10 = even
      lab0 <- 10^(min.pow:max.pow)
      axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
      axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
           labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
      xtick <- at0 <= max(df1[,"cond"]) & at0 >= min(df1[,"cond"])          # major labels
      ### if only two or fewer major lables, then add tick labels at 2 and 5 log
      if(sum(xtick)<=2) {##IF.xtick.START
        axis(1, at = log10(c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 1,
             labels = (c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))))
      }##IF.xtick.END
      axis(2)

      dev.set(3)
      val <- tapply(resp, df1$cutf, mean)
      plot(cutm, val, type = "n", ylim = range(c(val, plogis(up.bound.link), plogis(low.bound.link)),na.rm = T), axes= F,
           xlab = "" , ylab = "")
      points(cutm, val)
      lines(new.data$cond, plogis(up.bound.link), lty = 2)
      lines(new.data$cond, plogis(low.bound.link), lty = 2)
      lines(new.data$cond, mn)
      min.pow <- floor(min(df1$cond)); max.pow <- ceiling(max(df1$cond))
      at0 <- min.pow:max.pow              # add ticks at log10 = even
      lab0 <- 10^(min.pow:max.pow)
      axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
      axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
           labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
      xtick <- at0 <= max(df1[,"cond"]) & at0 >= min(df1[,"cond"])          # major labels
      ### if only two or fewer major lables, then add tick labels at 2 and 5 log
      if(sum(xtick)<=2) {
        axis(1, at = log10(c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 1,
             labels = (c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))))
      }
      axis(2)
      mtext(expression(paste("Specific conductivity ( ", mu, "S/cm)")), side = 1, line = 2.3, cex = 0.8)
      mtext("Probability of observing", side = 2, line = 2.3, cex =0.8)
      mtext(bquote(italic(.(mod.tnames[index]))), side = 3.5)
      box()
      abline(v = tolval.cdf[index], col = "red", lty = 2)

      if (index%%6==0 |index== length(mod.tnames)) graphics.off()
    }##IF.plot.END
  }##IF.index.START

  if(addtrend) {##IF.addtrend.START
    dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf), trend = trend,
                       XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
                       HUC = huc.df, eco3 = eco3.df, wtshed = wtshd.df, stringsAsFactors = F)
  } else {
    dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf),
                       XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
                       HUC = huc.df, eco3 = eco3.df, wtshed = wtshd.df, stringsAsFactors = F)
  }##IF.addtrend.END
  dftv <- dftv[order(dftv$XC95.cdf),]
#  write.table(wtshd.df,paste(wd, "/Results/wtshed.size.csv", sep= ""),sep=",", row.names = T)
  return(dftv)                     # return a data.frame with xc95s
}##FUNCTION.fish.wt.cdf.END
