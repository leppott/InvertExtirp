#  tolerance function was modified from taxon.response function 20130801
# a function to calculate species response curve using various models
# Get optimum value from taxon-environment relationship
# 1.WA, 2. cdf_Abundance 3. cdf_p/a. 4. ecdf weighted
# 5. Linear logistic 6. quadratic logistic 7. GAM  5~7 using full data range;
# 11. Count. 12. Raw quantiles
## spdata has
## envdata
# mtype = 1 linear, 2 quadratic 3 gam model for ploting purpose
# dense.N is the number of areas to cut into in the calculation of area under the curve
# plot.pdf to decide if we want specise vs. env plots
# log.x to decide if we want to log transform and plot data
# rounder how many digits
# bad is a boolean variabl to determine if the variable is bad for senstivie taxa
# taus determine the output the percentile of env variable
# r is a proportional ratio for predicted model output, as 0.95 is the 95% of optimum probability
#  coord c("LONG_DD", "LAT_DD")
#  if lim == "GAM", add gam plot xc95 otherwise, add   "CDF"
#  log.x to use log1o transformed data, plus to use log10(+1) transformed data, sqrt
  #graphics.off()
  #library(reshape)  #not sure if got all references, EWL, 20170329
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#' @title Tolerance
#'
#' @description A function to calculate species response curves using various models.
#' Get optimum value from taxon-environment relationship
#' 1.WA, 2. cdf_Abundance 3. cdf_p/a. 4. ecdf weighted
#' 5. Linear logistic 6. quadratic logistic 7. GAM  5~7 using full data range;
#'  11. Count. 12. Raw quantiles
#'  Requires the "reshape", "mgcv", and "maps" libraries.
#'
#' @param spdata Species data.
#' @param envdata Environmental data.
#' @param sp.siteid Species data frame column name for site/sample id; default = "Sample.ID"
#' @param species Species data frame column name for taxa name; default = ""GENUS"
#' @param covar covariance; default = NULL
#' @param sp.abndid Species data frame column name for species abundance; default = "RA"
#' @param env.siteid Environmental data frame column name for site/sample id.  Default is "Sample.ID"
#' @param xvar Environmental data frame column name for analyte used in analysis; default = "cond"
#' @param cutoff cutoff 1; default = 20
#' @param cutoff2 cutoff 2; Min number of samples for taxa occurrence.  Only those above cut off are plotted.  default = 10
#' @param region region.  A folder will be created in the working directory (wd) for the region; default = "all"
#' @param lim "GAM" or "CDF" ; default = "GAM"  if lim == "GAM", add gam plot xc95 otherwise, add   "CDF"
#' @param coord Map coordinates of samples; c("LONG_DD", "LAT_DD").  Default is NULL.
#' @param mtype 1 linear, 2 quadratic 3 gam model for ploting purpose; default = 3.
#' @param dense.N the number of areas to cut into in the calculation of area under the curve; default = 201.
#' @param cast ; default = TRUE,
#' @param plot.pdf to decide if we want specise vs. env plots; default = F
#' @param add.map boolean for if a map should be plotted.  Requires statename coord; default = F
#' @param statename State name (e.g., Maryland) used for map; default = NULL
#' @param add.lab boolean for adding labels to map; default = F
#' @param add.abund boolean for adding abundance to map; default = T.
#' @param main = "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient",
#' @param mar Plot margins; default = c(5, 4, 3, 4)
#' @param xlabs X-axis labels for plots; default = expression(paste("Conductivity ( ", mu, "S/cm)"))
#' @param log.x Whether to log transform and plot data; default = TRUE.  log.x to use log10 transformed data, plus to use log10(+1) transformed data, sqrt
#' @param plus ; Use with log.x=TRUE if want log10(+1) transformation.  default = F
#' @param rounder How many digits to round results; default = 0.
#' @param taus determine the output the percentile of env variable; default = c(50,95)
#' @param nbin number of bins; default = 61
#' @param wd Working directory for saving files.
#
#' @return Returns a data frame (cond.opt) and generates a PDF (results.taxon.gam.pdf) of plots (within the working directory in a subfolder for the region).
#
#' @keywords logistic regression, quantiles, xc95, hc05, cdf, gam, taxon response
#
#' @details Species data file can be long or wide format.  
#' If in long format use cast=TRUE to transform.
#' 
#' Cutoff2 will limit the analysis based on taxon occurrence across all sites.  
#' Only those above this values will have plots generated.
#' 
#' Cutoff is used for co-occurrence of taxa and the environmental variable.  
#' Taxa below this threshold will not have curves added to the plots.
#' 
#' Plots will be produced in a PDF and saved to the specified folder (region) based on user defined variables.
#' 
#' The output data frame will have the following fields.  
#' The variable 'taus' determines the quantiles for some calculations.
#' With the default values of 50 and 95 the 50th and 95th fields are shown below.
#' 
#' * tnames = taxa names derived from input data frame spdata
#' 
#' * N = Number of samples where taxon and environmental variable co-occur. samplen
#' 
#' * min_ob = minimum observation of environmental variable where taxon occurs. limits.
#' 
#' * X50_th_ob = 50th quantile observation of environmental variable where taxon occurs. limits.
#' 
#' * X95_th_ob = 95th quantile observation of environmental variable where taxon occurs.  limits.
#' 
#' * max_ob	= minimum observation of environmental variable where taxon occurs. limits.
#' 
#' * Opt_WA = Results model 1, weighted averaging, WA
#' 
#' * Tol_WA =  Results model 1, weighted averaging, tol
#' 
#' * CDF_50_th_Abund and CDF_95_th_Abund= Results model 2, cdf_Abundance,  df1[ic1,xvar]
#' 
#' * CDF_50_th_PA and CDF_95_th_PA= Results model 3, cdf_p/a,  df1[ic2,xvar]
#' 
#' * CDF_wt_50_th and CDF_wt_95_th = Results model 4, ecdf weighted, eout
#' 
#' * LRM_50_th and LRM_95_th= Results model 5, linear logistic,  lrm1.95f
#' 
#' * QLRM_50_th and QLRM_95_th = Results model 6, quadratic logistic, lrm2.95f
#' 
#' * Opt_qlrm = Results model 6, quadratic logistic, lrm2.opt
#' 
#' * Tol_qlrm = Results model 6, quadratic logistic, lrm2.tol
#' 
#' * GAM_50_th and GAM_95_th = Results model 7, GAM, lrm3.95f
#' 
#' * ROC = Summarizes all comparisons to compute area under ROC, roc
#
#' @examples
#' region = "results"
#' env.cond <- tol.env.cond
#' ss <- tol.ss
#' condunit <- expression(paste("Conductivity ( ", mu, "S/cm)", sep = ""))
#' cond.opt <- tolerance(spdata = ss, envdata = env.cond,  sp.siteid = "NewSampID", species = "Genus", covar=NULL,
#'                       sp.abndid = "RA", env.siteid = "NewSampID", xvar = "cond", cutoff = 25, cutoff2 = 10,
#'                       region = "results", lim ="GAM", coord = NULL, mtype = 3, dense.N = 201, cast = FALSE,
#'                       plot.pdf = T, add.map = F, statename = NULL,  add.lab = F, add.abund=T,
#'                       main = "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient",
#'                       mar = c(5,4,3,4), xlabs = condunit, log.x = T,
#'                       plus = F, rounder = 3, taus = c(50,95), nbin = 61, wd=getwd())
#' # view returned data frame
#' View(cond.opt)
#' # open saved PDF
#' system(paste0('open "',file.path(getwd(),region,"results.taxon.gam.pdf"),'"'))
#
#' @export
tolerance <- function(spdata, envdata,  sp.siteid = "Sample.ID", species = "GENUS", covar = NULL,
        sp.abndid="RA", env.siteid = "Sample.ID", xvar = "cond", cutoff = 20, cutoff2 = 10,
        region = "all", lim ="GAM", coord = NULL, mtype = 3, dense.N = 201, cast = TRUE,
        plot.pdf = F, add.map = F,statename = NULL,  add.lab = F, add.abund = T,
        main = "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient",
        mar = c(5,4,3,4), xlabs = expression(paste("Conductivity (", mu, "S/cm)")),
        log.x = TRUE, plus = F, rounder = 0, taus = c(50,95), nbin = 61, wd = getwd()) {##FUNCTION.tolerance.START

  # 20170413, add directory check
  dir.check.add(wd,region)
  #dir.check.add(file.path(wd,region),"cdf")
  #dir.check.add(file.path(wd,region),"gam")
  wd <- file.path(wd,region)

  dfenv <- envdata[!is.na(envdata[,xvar]), c(env.siteid, xvar,covar, coord)]

  if(cast) {##IF.cast.START
    df.sp <- data.frame(SITEID = spdata[, sp.siteid], species = spdata[, species],
                          abund = spdata[, sp.abndid])
    df.sp <- subset(df.sp, !is.na(species) & !is.na(abund))
    df.tmp <- merge(df.sp, dfenv, by.x = "SITEID", by.y = env.siteid )
    ss <- reshape::cast(df.tmp, SITEID ~ species, sum, value = "abund")
    names(ss)[1] <- sp.siteid
  } else ss <- spdata
  ##IF.cast.END

  taxa.count <- apply(ss[-1] > 0, 2, sum)  # same as colSums(ss[,-1]>0)
  tnames <- names(taxa.count)[taxa.count >= cutoff2]
  ntaxa <- length(tnames)                           # length of taxa list for taxonomic level i

    df1 <- merge(ss, dfenv, by.x = sp.siteid, by.y = env.siteid)   # merged crosstab abundance file
    nc <- ncol(df1)            # 1+ ntaxa + env
    nv <- ncol(dfenv) - 1     # nunmber of env var included in the data
    if(!is.null(covar)) {##IF.covar.START
      df1$cutcov <- cut(df1[,covar], quantile(df1[,covar], prob = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
          include.lowest =T)
      levels(df1$cutcov) <- 1:5
    } else {
      #20170414, EWL, add condition for if covar is NULL
      df1$cutcov <- NA
    }##IF.covar.END

    df2 <- df1 <- df1[order(df1[,xvar]),]
    df2[2:(nc- nv)] <- apply(df2[2:(nc - nv )] > 0, 2, as.numeric)                 # p/a file

    stress <- df1[, xvar]
    stress.u <- sort(unique(stress))

    xrange <- range(stress)
    xnew <- seq(from = xrange[1], to = xrange[2], length = 100)

    ###### calculate weight for cdf

    cutp <- seq(from = min(stress), to = max(stress), length = nbin)
    cutm <- 0.5*(cutp[-1] + cutp[-nbin])
    df2$cutf <- cut(df2[,xvar], cutp, include.lowest = T)

    wt <- 1/table(df2$cutf)
    wt[is.infinite(wt)] <- NA
    wt <- wt/sum(wt, na.rm = T)
    dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
    df2 <- merge(df2, dftemp, by = "cutf")[-1]
    df2 <- df2[order(df2[,xvar]),]

   #### end


#    ##### function to compute area under the curve
#     auc <- function(xrange, mod, dense.N) {##FUNCTION.auc.START
#         x <- seq(min(xrange), max(xrange), length = dense.N - 1)
#         s.area <-rep(NA, dense.N - 1)
#         y <- predict(mod, newdata = data.frame(dose = x), type = "response")
#         for (index in 1:dense.N-1) {
#           s.area[index] <- (y[index] + y[index + 1])/2*(x[index+1]-x[index])
#           }
#         tsum <- sum(s.area, na.rm=T)         # all area
#         jj=1
#         csum = sum(s.area[1:jj])              # cum area
#         xc95 <- yc95 <- rep(NA, length(taus))
#         for(ii in 1:length(taus)) {
#         while(csum < taus[ii]/100*tsum) {
#            jj = jj + 1
#            csum <- sum(s.area[1:jj], na.rm=T)
#           }
#         if(jj == 1) {
#           xc95[ii] <- x[jj]
#           yc95[ii] <- y[jj]
#           } else {
#          xc95[ii] <- (x[jj]+ x[jj-1])/2
#          yc95[ii] <- (y[jj]+ y[jj-1])/2
#          }
#         }
#         return(c(xc95))
# #      print(xc95)
#     }##FUNCTION.auc.END

# plot pdf option
    if(plot.pdf) {##IF.plot.cdf.START
    pdf(file = paste(wd,"/", region, ".taxon.gam.pdf",sep=""),
           width = 9, height = 6.5, pointsize = 12)
    par(mfrow = c(2, 3), pty = "m", mar =mar, oma=c(1.5, 0.5, 2,0.5) )
    totpage <- ceiling((length(tnames)+1)/6)

    if(add.map) {##IF.add.map.START
       simpleCap <- function(x) {##FUNCTION.simpleCap.START
               s <- strsplit(x, " ")[[1]]
               paste(toupper(substring(s, 1,1)), substring(s, 2),
               sep="", collapse=" ")
       }##FUNCTION.simpleCap.END
      med = 1
      while(is.na(df1[med, coord[1]])) med = 1 + med

      if(is.null(statename)) {
      statename <<- maps::map.where("state", df1[med, coord[1]], df1[med, coord[2]])
      statename <- simpleCap(statename)
      }

       if(add.lab) {
       maps::map.text("state", region = statename, mar = mar)
       } else {
       maps::map("state", regions = statename, mar = mar)
       }

       text(df1[,coord[1]],df1[,coord[2]],".")
#       mtext(paste("Ecoregion",region), side =3,line = 0.5, cex=1.5)
     }##IF.add.map.END
    }##IF.plot.cdf.END

    varnames <- c("N", "min_ob", paste(taus, "th_ob",sep="_"), "max_ob", "Opt_WA",
      "Tol_WA", paste("CDF", taus,"th","Abund", sep ="_"),
      paste("CDF", taus,"th","PA", sep ="_"),paste("CDF_wt", taus,"th",sep ="_"),
      paste("LRM", taus,"th", sep ="_"),paste("QLRM", taus,"th", sep ="_"),
      "Opt_qlrm", "Tol_qlrm", paste("GAM",taus, "th", sep="_"), "ROC")

    optsave <- matrix(NA, ncol = length(taus)*7 + 8, nrow  = ntaxa)
    colnames(optsave) <- varnames
    rownames(optsave) <- tnames

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:ntaxa) {##IF.i.START
      isel <- match(tnames[i], names(df1))       # selected taxa i
    ### (0) observed range
      samplen <- length(df2[df2[,isel]>0, xvar])
      limits <- quantile(df2[df2[,isel]>0, xvar], prob = c(0, taus/100, 1))

    ## (1)  Weighted averaging
      WA <- sum(df1[,isel]*df1[,xvar])/sum(df1[,isel])
      tol <- sqrt(sum(df1[,isel] * (df1[,xvar]- WA)^2) /sum(df1[,isel]))
  #    WTA <- sum(df1[,isel]*df1[,xvar]*df2$wt)/sum(df1[,isel]*df2$wt)     # sample weighted weighted average

    ### (2)  cdf or cumsum no weighting  both abundance based and p/a based
      csum1 <-cumsum(df1[,isel])/sum(df1[,isel])          # cumlative vectors divided total abundance
      ic1 <- rep(1, length(taus))
      print(c(i, ic1) )
      for(ii in 1:length(taus)) {##FOR.ii.START
        if(!is.na(sum(csum1))==TRUE) {
          while(csum1[ic1[ii]] < taus[ii]/100) ic1[ii] <- ic1[ii] + 1
        }
      }##FOR.ii.END
      # 20170414, EWL, need condition for NA or NaN, else the while fails inside the for loop above
      # added if !is.na for sum of csum1

    #### (3) cdf p/a
      csum2 <-cumsum(df2[,isel])/sum(df2[,isel])    # cumlative vectors devided total abundance
      ic2 <- rep(1, length(taus))
      for(ii in 1:length(taus)) {
        if(!is.na(sum(csum2))==TRUE) {
          while(csum1[ic2[ii]] < taus[ii]/100) ic2[ii]<- ic2[ii] + 1
        }
      }
      # 20170414, EWL, repeat extra check for same error as above


    ##### (4) cdf weighted
      pres <- df2[, isel] > 0
      sel <- df2[pres,]
      eout2 <- rep(NA, length(taus))
      for(ii in 1:length(taus)) {
        eout2[ii] <- Hmisc::wtd.quantile(sel[,xvar], sel$wt, normwt = TRUE, prob = taus[ii]/100)
      }

    ######### (5)(6)(7) logistic regression model
      resp <- df2[, isel]
      dose <- df2[, xvar]
      covwt <- as.numeric(df2[, "cutcov"])
      lrm1 <- glm(resp ~ dose, family="binomial")
      lrm2 <- glm(resp ~ dose + I(dose^2), family="binomial")
      #### u = -b1/(2b2)  t = 1/sqrt(-2b2) for lrm2
      # Remove the 2 lines below, per Lei, 20170407
      # also removes 2 columns in output
      #lrm2.opt <- (-lrm2$coef[2]/(2*lrm2$coef[3])  )
      #lrm2.tol <- (1/sqrt(-2*lrm2$coef[3]))
      lrm2.opt <- NA
      lrm2.tol <- NA
      lrm3 <- mgcv::gam(resp ~ s(dose,k = 3), family= "binomial")
      if(mtype == 1) model <- lrm1
      if(mtype == 2) model <- lrm2
      if(mtype == 3) model <- lrm3

     # Compute mean predicted probability of occurrence
      predlk <- predict(model, type = "link")
      predout <- exp(predlk)/(1 + exp(predlk))
     # Generate logical vector corresponding to presence/absence
     # Divide predicted probabilities into sites where
     # species is present ( x ) and sites where the species is
     # absent ( y ).
         x <- predout[pres]
         y <- predout[!pres]
     # Now perform all pairwise comparisons of x vs. y
     # and store results in a matrix
         rocmat <- matrix(NA, nrow = length(x), ncol = length(y))
         for (j in 1:length(x)) {
           if(!is.na(sum(csum2))==TRUE) {
              rocmat[j,] <- as.numeric(x[j] > y)
           }
         }
         # 20170414, add error condition



      # Summarize all comparisons to compute area under ROC
         roc <- sum(rocmat)/(length(x)*length(y))
    #        wilcox <- wilcox.test(x, y, conf.int = FALSE)
    #        AUC <- wilcox$statistic / (length(x)*length(y))  # same as roc
    #        print(roc)
    #        print(AUC)
      predresp <- predict(model, newdata = data.frame(dose=xnew), type="link", se.fit=T)

      # Compute upper and lower 90% confidence limits
      up.bound.link <- predresp$fit + qnorm(0.95)*predresp$se.fit
      low.bound.link <- predresp$fit + qnorm(0.05)*predresp$se.fit
      mean.resp.link <- predresp$fit

      # Convert from logit transformed values to probability.
      up.bound <- exp(up.bound.link)/(1+exp(up.bound.link))
      low.bound <- exp(low.bound.link)/(1+exp(low.bound.link))
      mean.resp <- exp(mean.resp.link)/(1+exp(mean.resp.link))

    #### 5,6,7 full range
      auc.version <- "tolerance"
      ######## find modeled xc95 full data range
      lrm1.95f <- auc(xrange = xrange, mod = lrm1, dense.N = dense.N, taus, auc.version)
      lrm2.95f <- auc(xrange = xrange, mod = lrm2, dense.N = dense.N, taus, auc.version)
      lrm3.95f <- auc(xrange = xrange, mod = lrm3, dense.N = dense.N, taus, auc.version)

    ###### Save results

      optsave[i,] <- c(samplen, round(c(limits, WA, tol, df1[ic1,xvar],
         df1[ic2,xvar], eout2, lrm1.95f,lrm2.95f, lrm2.opt, lrm2.tol, lrm3.95f),rounder), roc)
         #df1[ic2,xvar], eout2, lrm1.95f,lrm2.95f, lrm3.95f),rounder), roc)
          #20170414, EWL, remove 2 variables that Lei had me comment out above

      #### The following code calcualte probabilities of occurrence
      binf <- cut(df2[,xvar], cutp, include.lowest = T)
      bvals <- tapply(df2[, isel] > 0, binf, mean, na.rm =T)

    # optiion for plot single species env graph
     if(plot.pdf) {##IF.plot.pdf.START

       typ <- ifelse(samplen >= cutoff, "l", "n")
       plot(xnew, mean.resp, axes = FALSE, type = typ, xlab = "", ylab = "", xlim= xrange,
                  ylim = range(low.bound, up.bound, bvals, na.rm =T), col=1)

       if(samplen >= cutoff) {##IF.samplen.START
         points(xnew, low.bound, type="l", col=1, lty="dashed")
         points(xnew, up.bound, type="l", col=1,lty="dashed")
         if(lim =="GAM") {##IF.lim.START
            abline(v= lrm3.95f, lty= "dashed", col="red")
            } else if(lim =="CDF"){
            abline(v= eout2, lty= "dashed", col="blue")
            }##IF.lim.END
       }##IF.samplen.END

       max.pow <- ceiling(max(xrange)); min.pow <- floor(min(xrange))  ## add x range for log formation

       if(log.x) {##IF.log.x.START
            at0 <- min.pow:max.pow              # add ticks at log10 = even
            lab1 <- 10^at0
            # if log (x + 1) transformed then need to minus 1 from the orignial location
            axis(1, at = log10(lab1 + plus), labels = lab1, lwd.ticks = 1)    # major ticks with labels
            axis(1, at = log10((1:10 * rep(lab1[-1]/10, each = 10)) + plus), tcl = -0.4,
                labels = FALSE, lwd.ticks = 0.8)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
            xtick <- at0 <= max(df1[,xvar]) & at0 >= min(df1[,xvar])          # major labels
                       ### if only two or fewer major lables, then add tick labels at 2 and 5 log
            if(sum(xtick)==2) {
              at1 <- log10(3 * rep(lab1[-1]/10, each = 1) + plus)
              axis(1, at = at1, tcl = -0.3, lwd.ticks= 1, labels = (3 * rep(lab1[-1]/10, each = 1)),
               tcl = - 0.4)
            }   else if (sum(xtick)<2) {
              at2 <- log10(c(2, 5) * rep(lab1[-1]/10, each = 2) + plus)
              axis(1, at = at2, tcl = -0.3, lwd.ticks= 1, labels = (c(2, 5) *
                 rep(lab1[-1]/10, each = 2)))
            }
          } else  {
            axis(1, at = pretty(xrange), labels = pretty(xrange), lwd.ticks = 1, las = 1)    # major ticks with labels
       }##IF.log.x.END

       axis(2)
       y1range <- range(low.bound, up.bound, bvals, na.rm =T)
       y2range <- range(df1[, isel])

       if(add.abund) {##IF.add.abund.START
          y2new <- (df1[,isel] - min(y2range))*diff(y1range)/diff(y2range)+ min(y1range)    # convert y2 to y1 scale
          y1new.at <- (axTicks(2) - min(y1range))*diff(y2range)/diff(y1range)+ min(df1[, isel])    # convert to y2 scale
          y2.lab <- round(y1new.at, 2)    # add labels at original y2 scale
          axis(4, at = axTicks(2), labels = y2.lab, tcl =-0.4, lwd.ticks = 1)    # major ticks with labels
          if(is.null(covar)) {##IF.covar.START
            points(dose, y2new, pch = 21, col = 1, cex = .4, bg = "lightgray")
          } else {
            points(dose[y2new>0], y2new[y2new>0], pch = 21, col = 1, cex = as.numeric(covwt)*0.2, bg = "lightgray")
          }##IF.covar.END
        mtext("Relative Abundance", side = 4, line = 2.5, col= 1 , cex = 0.8)
        } else {
          points(cutm, bvals, pch=21, col= 1, cex = 0.7, bg ="gray")                      # probability
        }##IF.add.abund.END

       mtext(xlabs, side = 1, line = 2.3, cex = 0.9)
       mtext("Capture Probability", side = 2, line = 2.3, cex = 0.9, col=1)
       mtext(bquote(italic(.(tnames[i]))), side = 3, line = 0.5)
       box()

       if((i+5)/6-round((i+5)/6)==0) {##IF.i+5.START
          title(main= main, outer=TRUE, cex.main= 1.8)
          mtext(paste("Page", (i+5)/6, "of", totpage), outer = TRUE, side =1)
        }##IF.i+5.END

      }##IF.plot.pdf.END
   }##IF.i.END
    if(plot.pdf) graphics.off()
    return(optout <- data.frame(tnames, optsave))
}##FUNCTION.tolerance.END







