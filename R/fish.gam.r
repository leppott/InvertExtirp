#' GAM
#'
#'
#'
#' @param datafile1 Environmental data file 1.
#' @param datafile2 Environmental data file 2.
#' @param ss1 Species data file 1.
#' @param ss2 Species data file 2.
#' @param taxon Species data frame taxon name.
#' @param SampleID Site/sample id column; default = "FishInvertFID"
#' @param wtHUC Should data be weighted by HUC; default = TRUE
#' @param SampleID2 Sample id field; default = "Station_Date"
#' @param fit = TRUE
#' @param np Number of bins; default = 61.
#' @param nt Minimum number of occurence; default = 25.
#' @return **Need something here**
#' @keywords logistic regression, quantiles, xc95, hc05, cdf, gam, taxon response
#' @examples
#' #
#' @export
fish.gam <- function(datafile1, datafile2, ss1, ss2, taxon, SampleID = "FishInvertFID", wtHUC = TRUE,
    SampleID2 = "Station_Date", fit = TRUE, np = 61, nt = 25) {##FUNCTION.fish.gam.START

  my.ss <- merge(datafile1[c(SampleID, "HUC", "cond", "HUC04")], ss1)#; dim(my.ss)   # merge env data with species data
  my.ss2 <- merge(datafile2[c(SampleID2, "HUC", "cond", "BigHUC","WS_AREA")], ss2)#; dim(my.ss2)   # merge env data with species data

    if(wtHUC) {##IF.wtHUC.START
    HUCs <- unique(my.ss[my.ss[, taxon] > 0 ,"HUC04"])   # find all Site HUCs
    df1 <- subset(my.ss, HUC04 %in% HUCs)#; dim(df1)
    HUCs2 <- unique(my.ss2[my.ss2[,taxon] > 0 ,"BigHUC"])   # find all Site HUCs
    df2 <- subset(my.ss2, BigHUC %in% HUCs2)#; dim(df2)
    } else {
      df1 <- my.ss
      df2 <- my.ss2
    }##IF.wtHUC.END

  ###### data set 1, MN database
    cnew <- seq(from = min(df1$cond), to = max(df1$cond), length = 100)      # for plot
    new.data <- data.frame(cond = cnew)
    cond.u <- sort(unique(df1$cond))

    cutp <- seq(from = min(df1$cond), to = max(df1$cond), length = np)  # cut entire gradient into np bins
    cutm <- 0.5*(cutp[-1] + cutp[-np])                   # find middle point
    df1$cutf <- cut(df1$cond, cutp, include.lowest = T)   # define bins

    wt <- 1/table(df1$cutf)                 # weight for each bin
    wt[is.infinite(wt)] <- NA
    wt <- wt/sum(wt, na.rm = T)

    dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
    df1 <- merge(df1, dftemp, by = "cutf")#; dim(df1); names(df1)  # final data

    #  ECDF of observed conductivity
    resp <- df1[, taxon] > 0              # response variable
    df12 <- df1[resp,]

    tolval.cdf <- wtd.quantile(df12$cond, df12$wt, normwt = TRUE, prob = 0.95) # xc95 calculation

    mod <- gam(resp ~ s(cond, k = 3), data = df1, family = "binomial")
    predresp <- predict(mod, new.data, type = "link", se.fit = T)

  # Compute upper and lower 90% confidence limits
    up.bound.link <- predresp$fit + qnorm(0.95) * predresp$se.fit
    low.bound.link <- predresp$fit - qnorm(0.95) * predresp$se.fit

      # Convert from logit transformed values to probability.
      up.bound <- exp(up.bound.link)/(1 + exp(up.bound.link))
      low.bound <- exp(low.bound.link)/(1 + exp(low.bound.link))
      mn <- plogis(predresp$fit)
      bararea <- (cnew[2] - cnew[1])* mn         ### determine gam xc95 values
      tot <- sum(bararea)
      ii <- 1
      while(sum(bararea[1:ii]) < 0.95 * tot) ii <- ii + 1
      tolval.gam <- cnew[ii]

      val <- tapply(resp, df1$cutf, mean)

################# dataset 2

    cnew2 <- seq(from = min(df2$cond), to = max(df2$cond), length = 100)      # for plot
    new.data2 <- data.frame(cond = cnew2)
    cond.u <- sort(unique(df2$cond))

    cutp <- seq(from = min(df2$cond), to = max(df2$cond), length = np)  # cut entire gradient into np bins
    cutm2 <- 0.5*(cutp[-1] + cutp[-np])                   # find middle point
    df2$cutf <- cut(df2$cond, cutp, include.lowest = T)   # define bins

    wt <- 1/table(df2$cutf)                 # weight for each bin
    wt[is.infinite(wt)] <- NA
    wt <- wt/sum(wt, na.rm = T)

    dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
    df2 <- merge(df2, dftemp, by = "cutf")#; dim(df2); names(df2)  # final data

    #  ECDF of observed conductivity
    resp2 <- df2[, taxon] > 0              # response variable
    df22 <- df2[resp2,]

    tolval.cdf2 <- wtd.quantile(df22$cond, df22$wt, normwt = TRUE, prob = 0.95) # xc95 calculation
    mod2 <- gam(resp2 ~ s(cond, k = 3), data = df2, family = "binomial")
    predresp2 <- predict(mod2, newdata = new.data2, type = "link", se.fit = T)
      # Compute upper and lower 90% confidence limits
      up.bound.link2 <- predresp2$fit + qnorm(0.95) * predresp2$se.fit
      low.bound.link2 <- predresp2$fit - qnorm(0.95) * predresp2$se.fit

      # Convert from logit transformed values to probability.
      up.bound2 <- exp(up.bound.link2)/(1 + exp(up.bound.link2))
      low.bound2 <- exp(low.bound.link2)/(1 + exp(low.bound.link2))
      mn2 <- plogis(predresp2$fit)
      bararea <- (cnew2[2] - cnew2[1])* mn2         ### determine gam xc95 values
      tot <- sum(bararea)
      ii <- 1
      while(sum(bararea[1:ii]) < 0.95 * tot) ii <- ii + 1
      tolval.gam2 <- cnew[ii]
      val2 <- tapply(resp2, df2$cutf, mean)

  ############### plot MN dataset #########
   plot(cutm, val, type = "n", ylim = range(c(val,val2, up.bound, low.bound, up.bound2, low.bound2), na.rm = T),
     xlim = range(df1$cond, df2$cond, na.rm =T),
     axes= F, xlab = "" , ylab = "")
      min.pow <- floor(min(df1$cond, df2$cond)); max.pow <- ceiling(max(df1$cond, df2$cond))
      at0 <- min.pow:max.pow              # add ticks at log10 = even
      lab0 <- 10^(min.pow:max.pow)
      axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
      axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
           labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
      xtick <- at0 <= max(df1$cond, df2$cond) & at0 >= min(df1$cond, df2$cond)          # major labels
      ### if only two or fewer major lables, then add tick labels at 2 and 5 log
      if(sum(xtick)<=2) {##IF.xtick.START
        axis(1, at = log10(c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 1,
             labels = (c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))))
      }##IF.xtick.END
      axis(2)
      mtext(expression(paste("Specific conductivity ( ", mu, "S/cm)")), side = 1, line = 2.3, cex = 0.8)
      mtext("Probability of observing", side = 2, line = 2.3, cex =0.8)
      mtext(bquote(italic(.(taxon))), side = 3.5)
      box()
      abline(v = tolval.cdf, col = 1, lty = 2)
      abline(v = tolval.cdf2, col = "skyblue", lty = 2)
  ### number Applianchia dataset
  if(fit) {##IF.fit.START
     points(cutm, val)
      lines(new.data$cond, plogis(up.bound.link), lty = 2)
      lines(new.data$cond, plogis(low.bound.link), lty = 2)
      lines(new.data$cond, mn)

      points(cutm2, val2, pch = 2, col ="skyblue")
      lines(new.data2$cond, up.bound2, lty = 2, col = "skyblue")
      lines(new.data2$cond, low.bound2, lty = 2, col = "skyblue")
      lines(new.data2$cond, lty =1, mn2, col = "skyblue")

      } else {
      polygon(c(new.data$cond,rev(new.data$cond)), c(mn,rep(0, length(mn))), col = "black", border = "black", density = 20, angle = -45)
      polygon(c(new.data2$cond,rev(new.data2$cond)), c(mn2,rep(0, length(mn2))), col = "skyblue", border = "skyblue", density = 20, angle = 45)
      box(bty="l")
  }##IF.fit.END
}##FUNCTION.fish.gam.END
