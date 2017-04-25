#' Plot Taxa Response Curve and Calculate Species Extirpation Concentrations
#'
#' The output of this function to return 1.Weighted Average, 2. cdf_Abundance based, 3. cdf_ presence/absence based;
#' 4. ecdf weighted, 5. cdf weight new; 6. Linear logistic regression, 7. quadratic logistic 8. GAM  5~7 using full data range;
#' 9~11. repeat 6~8 but uses observed range for each single taxon; 12 Count. 13. Raw quantiles.  Requires Hmisc for wtd.quantile() and Ecdf() and mgcv to gam().
#'
#' @param df1 data frame
#' @param sort.vect when plot, sort the taxa list according to a vector called file called sort.vec
#' @param xvar. xvariable, could be column index or name
#' @param cutoff a required minimum sample size for calculation
#' @param region a subregion code to name the final output files
#' @param mtype could be 1 to 3, indicating which regression model to use; default = 3.
#' @param dense.N is the number of areas to cut into in the calculation of area under the curve
#' @param plot.pdf to decide if we want species vs. env plots options "none", "pdf", "tiff"
#' @param add.map to decide if a map should be added before plots.
#' @param maintext title of the multiplots area
#' @param nbin number of bins for logits
#' @param log.x if xvar should be logtransformated
#' @param rounder xvar rounder, default = 0
#' @param taus determine the output the percentile of env variable
#' @param wd Working directory for saving files.
#' @return Output to the screen for each taxon as it is completed.  CDF and GAM plots are saved to the specified directory in subfolders ("cdf" and "gam").
#' 1.Weighted Average, 2. cdf_Abundance based, 3. cdf_ presence/absence based;
#' 4. ecdf weighted, 5. cdf weight new; 6. Linear logistic regression, 7. quadratic logistic 8. GAM  5~7 using full data range;
#' 9~11. repeat 6~8 but uses observed range for each single taxon; 12 Count. 13. Raw quantiles
#' @keywords logistic regression, quantiles, xc95, hc05, cdf, gam, taxon response
#' @examples
#' switch0 <- 1
#' ecolab <- ifelse (switch0 ==1, "eco69", "eco70")
#' unitlab <- expression(paste("SO"[4]^{2-phantom()}," + HCO"[3]^{-phantom()}," (mg/L)"))
#' full.results <- taxon.response.sort(df1 = df1, xvar = "lgSO4HCO3", cutoff = 25, region = ecolab
#' , mtype = 3, dense.N = 201, plot.pdf = T, xlabs = unitlab, add.map = F, , maintext = ""
#' , GIS.cord = c("Long_DD", "Lat_DD"), log.x = TRUE, rounder = 0, taus = c(0,95,100), nbin = 61, sort.vect = taxalist
#' , wd=getwd())
#' # view results
#' View(full.results)
#' @export
taxon.response.sort <- function(df1 = df1, xvar="Conductivity", cutoff = 25, region = "all",
        mtype = 3, dense.N = 201, plot.pdf=F, xlabs="Specific conductivity (uS/cm)", add.map = FALSE,
        GIS.cord = c("LONG_DD", "LAT_DD"), extirpation = NULL,
        maintext = "Macroinvertebrates response to specific conductivity",
        log.x=TRUE, rounder=0, taus=c(0,95,100), nbin = 61, sort.vect = sort.vect, wd=getwd()) {##FUNCTION.taxon.response.sort.START

    tnames <- sort.vect
    ntaxa <- length(tnames)
    xlims <- range(df1[,xvar], na.rm = TRUE)
    df1[,xvar] <- (df1[,xvar] - xlims[1])/diff(xlims)
    if(is.null(extirpation)) extirpation <- xlims[2]
    extirpation <- (extirpation - xlims[1])/diff(xlims)

    df1 <- df1[order(df1[,xvar]),]
    stress <- df1[,xvar]
    stress.u <- sort(unique(stress))

    xrange <- range(stress)
    xnew <- seq(from = xrange[1], to = xrange[2], length = 100)
    ###################### A weighted function to calculate cumsum
    ecdf.w <- function (x,w) {##FUNCTION.ecdf.w.START
      w1 <- round(100000*w)
      x0 <- numeric(0)
      for (i in 1:length(x)) {
        x0 <- c(x0, rep(x[i], times = w1[i]))
      }
      x0 <- sort(x0)
      n <- length(x0)
      if (n < 1)
        stop("'x' must have 1 or more non-missing values")
      vals <- unique(x0)
      rval <- approxfun(vals, cumsum(tabulate(match(x0, vals)))/n,
                      method = "constant", yleft = 0, yright = 1,
                      f = 0, ties = "ordered")
      class(rval) <- c("ecdf", "stepfun", class(rval))
      attr(rval, "call") <- sys.call()
      rval
    }##FUNCTION.ecdf.w.END

    ###### calculate weight for cdf

    cutp <- seq(from = min(stress), to = max(stress), length = nbin)
    cutm <- 0.5*(cutp[-1] + cutp[-nbin])
    df1$cutf <- cut(df1[,xvar], cutp, include.lowest = T)

    wt <- 1/table(df1$cutf)
    wt[is.infinite(wt)] <- NA
    wt <- wt/sum(wt, na.rm = T)
    dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
    df1 <- merge(df1, dftemp, by = "cutf")#[-1]
    df1 <- df1[order(df1[,xvar]),]
#    print(names(df1))
   #### end
   # ##### function to compute area under the curve
   #  auc <- function(xrange, mod, dense.N) {##FUNCTION.auc.START
   #      x <- seq(min(xrange), max(xrange), length = dense.N - 1)
   #      s.area <-rep(NA, dense.N - 1)
   #      y <- predict(mod, newdata = data.frame(dose = x), type = "response")
   #      for (index in 1:dense.N-1) {
   #        s.area[index] <- (y[index] + y[index + 1])/2*(x[index+1]-x[index])
   #        }
   #      tsum <- sum(s.area, na.rm=T)
   #      jj=1
   #      csum = sum(s.area[1:jj])
   #      while(csum < taus[2]/100*tsum) {
   #         jj = jj + 1
   #         csum <- sum(s.area[1:jj], na.rm=T)
   #      }
   #
   #      xc95 <- (x[jj]+ x[jj-1])/2
   #      yc95 <- (y[jj]+ y[jj-1])/2
   #      return(c(xc95,yc95))
   #  }##FUNCTION.auc.END

    # plot pdf option
    simpleCap <- function(x) {##FUNCTION.simpleCap.START
               s <- strsplit(x, " ")[[1]]
               paste(toupper(substring(s, 1,1)), substring(s, 2),
               sep="", collapse=" ")
    }##FUNCTION.simpleCap.END

     totpage <- ceiling((length(tnames)+1)/6)
     optsave <- matrix(NA, ncol = length(taus) + 13, nrow  = ntaxa)
     tolcl <- rep(NA, times = ntaxa)
     roc <- rep(NA, times = ntaxa)
     trend <- rep(NA, length(ntaxa))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     for (i in 1:ntaxa) {##FOR.i.START
      isel <- match(tnames[i], names(df1))       # selected taxa i
        print(i)

    ## (1)  Weighted averaging
        WA <- sum(df1[,isel] * df1[,xvar])/sum(df1[,isel])
    ### (2)  cdf or cumsum no weighting  both abundance based and p/a based
        csum1 <- cumsum(df1[,isel])/sum(df1[,isel])          # cumlative vectors devided total abundanc
        ic1 <- 1
        while(csum1[ic1] < taus[2]/100) ic1<- ic1 + 1
    #### (3) cdf p/a
        csum2 <-cumsum(df1[,isel])/sum(df1[,isel])          # cumlative vectors devided total abundance
        ic2 <- 1
        while(csum2[ic2] < taus[2]/100) ic2<- ic2 + 1

    ##### (4) cdf weighted
        resp <- df1[, isel] > 0
        sel <- df1[resp,]
        eout <- ecdf.w(sel[,xvar], sel$wt)
        eout2 <- Hmisc::wtd.quantile(sel[,xvar], sel$wt, normwt = TRUE, prob = taus[2]/100)
        ic3 <- 1; ic4 <- 1
        while(eout(stress.u[ic3]) < taus[2]/100) ic3 <- ic3 + 1
        while(eout(stress.u[ic4]) < 0.5) ic4 <- ic4 + 1

    ######### (5)(6)(7) logistic regression model

    #      resp <- df1[, isel] >0
        dose <- df1[, xvar]

        lrm1 <- glm(resp ~ dose, family="binomial")
        lrm2 <- glm(resp ~ dose + I(dose^2), family="binomial")
        lrm3 <- mgcv::gam(resp ~ s(dose,k=3), family= "binomial")
        if(mtype == 1) model <- lrm1
        if(mtype == 2) model <- lrm2
        if(mtype == 3) model <- lrm3
        predresp<- predict(model, newdata = data.frame(dose=xnew), type="link", se.fit=T)
        # Compute upper and lower 90% confidence limits
        up.bound.link <- predresp$fit + qnorm(0.95) * predresp$se.fit
        low.bound.link <- predresp$fit - qnorm(0.95) * predresp$se.fit
        mean.resp.link <- predresp$fit

         # Convert from logit transformed values to probability.
        up.bound <- exp(up.bound.link)/(1 + exp(up.bound.link))
        low.bound <- exp(low.bound.link)/(1 + exp(low.bound.link))
        mean.resp <- exp(mean.resp.link)/(1 + exp(mean.resp.link))

        ##### predict
        predresp.vline <- predict(model, newdata = data.frame(dose = extirpation), type = "link", se.fit = T)
        low.lk <- predresp.vline$fit - qnorm(0.95) * predresp.vline$se.fit
        low.extirp <- exp(low.lk)/(1 + exp(low.lk))
        mean.lk <- predresp.vline$fit
        mean.extirp <- exp(mean.lk)/(1 + exp(mean.lk))

        crit <- 0.01 * max(mean.resp)            # probably used 0.05 for 500 and 300
        if(mean.extirp <= crit) {##IF.mean.extirp.START
        trend[i] <- "="
        }  else  if(low.extirp <= crit & mean.extirp > crit ) {
        trend[i] <- "~"
        }  else {
        trend[i] <- ">"
        }##IF.mean.extirp.END

    # ########  assess model performance
        predout <- predict(model, type = "response")

        xx <- predout[resp]
        yy <- predout[! resp]
        rocmat <- matrix(NA, nrow = length(xx), ncol = length(yy))
        for (j in 1:length(xx)) {
         rocmat[j,] <- as.numeric(xx[j] > yy)
        }
    # Summarize all comparisons to compute area under ROC
        roc[i] <- sum(rocmat)/(length(xx)*length(yy))
    ### curve shape
        tolcl[i] <- curve.shape(mean.resp, up.bound,low.bound)
    #### 5,6,7 full range
        auc.version <- "response"
        ######## find modeled xc95 full data range
        lrm1.95f <- auc(xrange = xrange, mod = lrm1, dense.N = dense.N, taus, auc.version)
        lrm2.95f <- auc(xrange = xrange, mod = lrm2, dense.N = dense.N, taus, auc.version)
        lrm3.95f <- auc(xrange = xrange, mod = lrm3, dense.N = dense.N, taus, auc.version)
    #### 8,9,10 observed range
        ######## find modeled xc95 observed data range
        shortrange <- range(df1[df1[,isel]>0,xvar])
        lrm1.95p <- auc(xrange = shortrange, mod = lrm1, dense.N = dense.N, taus, auc.version)
        lrm2.95p <- auc(xrange = shortrange, mod = lrm2, dense.N = dense.N, taus, auc.version)
        lrm3.95p <- auc(xrange = shortrange, mod = lrm3, dense.N = dense.N, taus, auc.version)
    ##### 11  n values  12 quantiles
        samplen <- sum(df1[,isel]>0)
        limits<- quantile(df1[df1[,isel]>0, xvar], probs=taus/100)   # quantile calculation

    ###### Save results

       if(log.x) {##IF.log.x.START
        optsave[i,1:(length(taus)+ 12)] <- round(10^(c(WA,df1[ic1,xvar],df1[ic2,xvar],stress.u[ic3],eout2,
           stress.u[ic4],lrm1.95f[1],lrm2.95f[1],lrm3.95f[1],lrm1.95p[1],lrm2.95p[1],lrm3.95p[1], limits)*
           diff(xlims) + xlims[1]), rounder)             # WA optimum
        } else {
        optsave[i,1:(length(taus)+12)] <- round(c(WA,df1[ic1,xvar],df1[ic2,xvar],stress.u[ic3],eout2,
           stress.u[ic4],lrm1.95f[1],lrm2.95f[1],lrm3.95f[1],lrm1.95p[1],lrm2.95p[1],lrm3.95p[1], limits)*
           diff(xlims) + xlims[1], rounder)
        }##IF.log.x.END
        optsave[i,(length(taus) + 13)] <- c(samplen)

        bvals <- tapply(df1[,isel] >0, df1$cutf, mean)
        bvals[is.na(bvals)] <- 0

    # optiion for plot single species env graph
        if(plot.pdf) {##IF.plot.pdf.START
          if(i%%6==1 ) {
            pgs <- (i-1)%/%6 + 1

           # 20170406, add directory check
            dir.check.add(wd,region)
            dir.check.add(file.path(wd,region),"cdf")
            dir.check.add(file.path(wd,region),"gam")
           #
           tiff(file = paste(wd,"/", region,"/cdf/", pgs, ".taxon.cdf.tiff",sep=""),
                 width = 650, height = 450, pointsize = 13)
            par(mfrow = c(2, 3), pty = "m", mar = c(4, 4, 3, 1))

            tiff(file = paste(wd,"/", region,"/gam/", pgs, ".taxon.gam.tiff",sep=""),
                 width = 650, height = 450, pointsize = 13)
            par(mfrow = c(2, 3), pty = "m", mar = c(4, 4, 3, 1))
          }
         name.lab <- simpleCap(as.vector(tnames[i]) )
         dev.set(3)
         plot(xnew, mean.resp, axes = FALSE, type = "l", xlab = xlabs,
              ylab = "Probability of observing", xlim=c(0,1),
                    ylim = range(low.bound, up.bound, bvals, na.rm =T), col=1)
         points(xnew, low.bound, type="l", col=1, lty="dashed")
         points(xnew, up.bound, type="l", col=1,lty="dashed")
         points(cutm, bvals, pch=21, col= 1, cex = 0.7, bg ="gray")                      # probability

         abline(v= eout2, lty= "dashed", col="red")

         max.pow <- ceiling(max(xlims)); min.pow <- floor(min(xlims))  ## add x range for log formation

         if(log.x) {##IF.log.x.START
              at0 <- (min.pow:max.pow - xlims[1])/diff(xlims)             # add ticks at log10 = even
              lab1 <- 10^(min.pow:max.pow)
              axis(1, at = at0, labels = lab1, lwd.ticks = 1)    # major ticks with labels
              axis(1, at = (log10(1:10 * rep(lab1[-1]/10, each = 10)) -xlims[1])/diff(xlims), tcl = -0.3,
                  labels = FALSE, lwd.ticks = 0.8)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
              xtick <- at0 <= max(df1[,xvar]) & at0 >= min(df1[,xvar])          # major labels
                         ### if only two or fewer major lables, then add tick labels at 2 and 5 log
              if(sum(xtick)<=2) {
                axis(1, at = (log10(c(2, 5) * rep(lab1[-1]/10, each = 2)) -xlims[1])/diff(xlims), tcl = -0.4,lwd.ticks= 1,
                   labels = (c(2, 5) * rep(lab1[-1]/10, each = 2)))
              }
          } else  {
              at0 <- (pretty(xlims) - xlims[1])/diff(xlims)             # add ticks
  #            lab1 <- round(at0*diff(xlims) + xlims[1], digits = rounder)    # scale back to orginial value
              axis(1, at = at0, labels = pretty(xlims), lwd.ticks = 1.2, las = 1)    # major ticks with labels
          }##IF.log.x.END
         axis(2)
         mtext(bquote(italic(.(tnames[i]))), side = 3, line = 0.5)
         box()

    #### plot cdf
         dev.set(2)
         Hmisc::Ecdf(sel[,xvar], weights = sel$wt, xlim = c(0,1), col = "black", pch = 1,
              axes = F, main = bquote(italic(.(name.lab))), xlab = xlabs, ylab ="Cumulative proportion")
        max.pow <- ceiling(max(xlims)); min.pow <- floor(min(xlims))  ## add x range for log formation

        if(log.x) {##IF.log.x.START
            at0 <- (min.pow:max.pow - xlims[1])/diff(xlims)             # add ticks at log10 = even
            lab1 <- 10^(min.pow:max.pow)
            axis(1, at = at0, labels = lab1, lwd.ticks = 1)    # major ticks with labels
            axis(1, at = (log10(1:10 * rep(lab1[-1]/10, each = 10)) -xlims[1])/diff(xlims), tcl = -0.3,
                    labels = FALSE, lwd.ticks = 0.8)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
            xtick <- at0 <= max(df1[,xvar]) & at0 >= min(df1[,xvar])          # major labels
        if(sum(xtick)<=2) {
            axis(1, at = (log10(c(2, 5) * rep(lab1[-1]/10, each = 2)) -xlims[1])/diff(xlims), tcl = -0.4,lwd.ticks= 1,
              labels = (c(2, 5) * rep(lab1[-1]/10, each = 2)))
           }
        } else  {
           at0 <- (pretty(xlims) - xlims[1])/diff(xlims)             # add ticks
           axis(1, at = at0, labels = pretty(xlims), lwd.ticks = 1.2, las = 1)    # major ticks with labels
        }##IF.log.x.END
         axis(2)
         abline(h = c(taus[2]/100), col = "red", lty = 2)
         abline(v = eout2, col = "red", lty = 2)
         box(bty = "l")
        if(i%%6==0 |i== ntaxa) graphics.off()
      }##IF.plot.pdf.END
    }##FOR.i.END
    optout <- data.frame(tnames, optsave, roc, tolcl, trend)
    names(optout) <- c("taxaname", "WAopt","CDF_ABD","CDF_PA", "CDF_WT", "XC95", "CDF_WT05",
       "LRM_Full", "QLRM_Full","GAM_Full", "LRM_Ob", "QLRM_Ob", "GAM_Ob",
        paste("Observ",taus, sep=""),"N", "ROC", "CurveShape", "Trend")
    return(optout)
}##FUNCTION.taxon.response.sort.END
