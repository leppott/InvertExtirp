#' Plot Taxa Response Curve and Calculate Species Extirpation Concentrations
#'
#' The output of this function to return 1.Weighted Average, 2. cdf_Abundance based, 3. cdf_ presence/absence based;
#' 4. ecdf weighted, 5. cdf weight new; 6. Linear logistic regression, 7. quadratic logistic 8. GAM  5~7 using full data range;
#' 9~11. repeat 6~8 but uses observed range for each single taxon; 12 Count. 13. Raw quantiles
#'
#' @param spdata Species data.
#' @param envdata Environmental data.
#' @param sp.siteid Site/sample id column; default = "Sample.ID"
#' @param species ; default = "GENUS"
#' @param sp.abndid ; default = "RA"
#' @param env.siteid ; default = "Sample.ID"
#' @param xvar ; default = "COND"
#' @param cutoff ; default = 30
#' @param region ; default = "all"
#' @param lim ; default ="CDF"
#' @param coord ; default = NULL
#' @param mtype ; default = 3
#' @param dense.N ; default = 201
#' @param plot.pdf ; default = F
#' @param add.map ; default =F
#' @param statename ; default = NULL
#' @param add.lab ; default = F
#' @param main Plot title (main); default = "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient"
#' @param mar ; default = c(5,4,3, 2)
#' @param xlabs ; default = expression(paste("Conductivity ( ", mu, "S/cm)"))
#' @param log.x ; default = TRUE
#' @param rounder xvar rounder, default = 0.
#' @param taus ; default = c(0,95,100)
#' @param nbin Number of bins; default = 61.
#~~~~~~~~~~~~
# Lei had the ones below but they don't match the function
# @param df1, data frame
# @param sort.vect when plot, sort the taxa list according to a vector called file called sort.vec
# @param xvar. xvariable, could be column index or name
# @param cutoff a required minimum sample size for calculation
# @param region a subregion code to name the final output files
# @param mtype could be 1 to 3, indicating which regression model to use.
# @param dense.N is the number of areas to cut into in the calculation of area under the curve
# @param plot.pdf to decide if we want species vs. env plots options "none", "pdf", "tiff"
# @param add.map to decide if a map should be added before plots.
# @param maintext title of the multiplots area
# @param nbin number of bins for logits
# @param log.x if xvar should be logtransformated
# @param rounder xvar rounder, default is 0
# @param taus determine the output the percentile of env variable
# @param lim if lim == "GAM", add gam plot xc95 otherwise, add   "CDF"
#~~~~~~~~~~~~~~
#' @return **Need something here**
#' @keywords logistic regression, quantiles, xc95, hc05, cdf, gam, taxon response
#' @examples
#' #
#' @export
taxon.response <- function(spdata, envdata,  sp.siteid="Sample.ID", species="GENUS",
        sp.abndid="RA", env.siteid="Sample.ID", xvar="COND", cutoff=30, region = "all", lim ="CDF",
        coord = NULL, mtype = 3, dense.N = 201, plot.pdf = F, add.map=F,statename = NULL,  add.lab = F,
        main = "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient", mar = c(5,4,3, 2),
        xlabs=expression(paste("Conductivity ( ", mu, "S/cm)")), log.x=TRUE, rounder=0, taus=c(0,95,100), nbin = 61) {##FUNCTION.taxon.response.START

  graphics.off()
  dfenv <- envdata[, c(env.siteid, xvar, coord)]
  xlims <- range(dfenv[,xvar], na.rm = TRUE)
  dfenv[,xvar] <- (dfenv[,xvar] - xlims[1])/diff(xlims)
#  xrange <- dfenv[,xvar]
  bcnt0 <- merge(spdata, envdata, by.x = sp.siteid, by.y = env.siteid)       # big table merged togetther
 # print(dim(bcnt0))
  bcnt0 <- (bcnt0[bcnt0[, sp.abndid]>0,])
 # print(dim(bcnt0))
  # assign station names to sitenames
  sitenames.u <- sort(unique(unlist(bcnt0[,sp.siteid])))                        # unique site names

  nsamp <- length(sitenames.u)           # count total number of stations
  luniq <- function(x) length(unique(x))     # function to count unique id

  numocc <- tapply(bcnt0[, sp.siteid], bcnt0[,species], luniq)  # No. of occurance(unique stations) for a list of taxa
#  numocc2 <- nsamp - numocc                                       # total stations- vector of stations for all taxa
#  numocc3 <- pmin(numocc, numocc2, na.rm= TRUE)                   # find the min of present and absent's sites
  numocc <- numocc[!is.na(numocc)]
  incvec <- numocc >= cutoff                                    # presence/absence min >= cutoff value

  tnames <- sort(names(numocc)[incvec])                          # selected taxa list
  ntaxa <- length(tnames)                           # length of taxa list for taxonomic level i

# Build site species matrix for only taxa with enough observations

    abund.all <- bcnt0[, sp.abndid]                           # abundance only
    ss1 <- matrix(0, ncol = ntaxa, nrow = nsamp)              # build a matrix to store crosstab table
    for (j in 1:ntaxa) {##FOR.j.START                         # all taxa in j list
      incvec <- bcnt0[, species] == tnames[j]        # match a taxon name at the selected taxonomic level with tnames.loc[j]
      incvec[is.na(incvec)] <- FALSE
      sitenames.g <- bcnt0[incvec,sp.siteid]           # selected rows/station for taxon j
      abund.g <- bcnt0[incvec,sp.abndid]               # selected abundance data for taxon j
      abund.s <- tapply(abund.g, sitenames.g, sum)       # crosstab total abundance for taxon j
      abund.s[is.na(abund.s)] <- 0
      sitenames.s <- names(abund.s)
         for (k in 1:length(sitenames.s)) {
         ss1[match(sitenames.s[k], sitenames.u), j] <- abund.s[k]
         }
    }##FOR.j.END

    ss <- data.frame(sitenames.u, ss1)
    names(ss) <- c("SITEID", tnames)
    df1 <- merge(ss, dfenv, by.x = "SITEID", by.y = env.siteid)   # merged crosstab abundance file

    df1[is.na(df1[2:(ntaxa + 1)]),2:(ntaxa + 1)] <- 0
    nc <- ncol(df1)
    df1 <- df1[order(df1[,xvar]),]
    df2 <- df1
    df2[2:(ntaxa + 1)][df2[2:(ntaxa + 1)] > 0] <- 1                 # p/a file
    stress <- df1[,xvar][!is.na(df1[,xvar])]
    stress.u <- sort(unique(stress))

    xrange <- range(stress)
    xnew <- seq(from = xrange[1], to = xrange[2], length = 100)

    ###################### A weighted function to calculate cumsum
    ecdf.w <- function (x,w) {##FUNCTION.ecdf.w.START
      w1 <- round(1000000*w)
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
    df2$cutf <- cut(df2[,xvar], cutp, include.lowest = T)

    wt <- 1/table(df2$cutf)
    wt[is.infinite(wt)] <- NA
    wt <- wt/sum(wt, na.rm = T)
    dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
    df2 <- merge(df2, dftemp, by = "cutf")[-1]
    df2 <- df2[order(df2[,xvar]),]

   #### end
   ##### function to compute area under the curve
    auc <- function(xrange, mod, dense.N) {##FUNCTION.auc.START
        x <- seq(min(xrange), max(xrange), length = dense.N - 1)
        s.area <-rep(NA, dense.N - 1)
        y <- predict(mod, newdata = data.frame(dose = x), type = "response")
        for (index in 1:dense.N-1) {
          s.area[index] <- (y[index] + y[index + 1])/2*(x[index+1]-x[index])
          }
        tsum <- sum(s.area, na.rm=T)
        jj=1
        csum = sum(s.area[1:jj])
        while(csum < taus[2]/100*tsum) {
           jj = jj + 1
           csum <- sum(s.area[1:jj], na.rm=T)
        }

        xc95 <- (x[jj]+ x[jj-1])/2
        yc95 <- (y[jj]+ y[jj-1])/2
        return(c(xc95,yc95))
    }##FUNCTION.auc.END

  #    print(dim(df1))
  # plot pdf option
    if(plot.pdf) {##IF.plot.pdf.START
    pdf(file = paste(wd,"/", region, ".taxon.cdf.pdf",sep=""),
           width = 9, height = 6, pointsize = 12)
    par(mfrow = c(2, 3), pty = "m", mar = mar, oma=c(1.5, 0.5, 2,0.5) )
    pdf(file = paste(wd,"/", region, ".taxon.gam.pdf",sep=""),
           width = 9, height = 6.5, pointsize = 12)
    par(mfrow = c(2, 3), pty = "m", mar =mar, oma=c(1.5, 0.5, 2,0.5) )
    totpage <- ceiling((length(tnames)+1)/6)

    if(add.map) {
       simpleCap <- function(x) {
               s <- strsplit(x, " ")[[1]]
               paste(toupper(substring(s, 1,1)), substring(s, 2),
               sep="", collapse=" ")
              }
      med = 1
      while(is.na(df1[med, coord[1]])) med = 1 + med

      if(is.null(statename)) {
      statename <<- map.where("state", df1[med, coord[1]], df1[med, coord[2]])
      statename <- simpleCap(statename)
      }
       if(add.lab) {
       map.text("state", region = statename, mar = mar)
       } else {
       map("state", regions = statename, mar = mar)
       }

       text(df1[,coord[1]],df1[,coord[2]],".")
#       mtext(paste("Ecoregion",region), side =3,line = 0.5, cex=1.5)
     }
    }##IF.plot.pdf.END

    optsave <- matrix(NA, ncol = length(taus) + 12, nrow  = ntaxa)

    ##~~~~~~~~~~~~~~~
    for (i in 1:ntaxa) {##FOR.i.START
      isel <- match(tnames[i], names(df1))       # selected taxa i

  ## (1)  Weighted averaging
      WA <- sum(df1[,isel]*df1[,xvar])/sum(df1[,isel])
      tol <- sqrt(sum(df1[,isel] * (df1[,xvar]- WA)^2) /sum(df1[,isel]))
  #    WTA <- sum(df1[,isel]*df1[,xvar]*df2$wt)/sum(df1[,isel]*df2$wt)     # sample weighted weighted average
  ### (2)  cdf or cumsum no weighting  both abundance based and p/a based
      csum1 <-cumsum(df1[,isel])/sum(df1[,isel])          # cumlative vectors devided total abundance
      ic1 <- 1
      while(csum1[ic1] < taus[2]/100) ic1<- ic1 + 1

  #### (3) cdf p/a
      csum2 <-cumsum(df2[,isel])/sum(df2[,isel])          # cumlative vectors devided total abundance
      ic2 <- 1
      while(csum2[ic2] < taus[2]/100) ic2<- ic2 + 1


  ##### (4) cdf weighted
      resp <- df2[, isel] > 0
      sel <- df2[resp,]
#      eout <- ecdf.w(sel[,xvar], sel$wt)
#      ic3 <- 1
#      while(eout(stress.u[ic3]) < taus[2]/100) ic3 <- ic3 + 1
      eout2 <- wtd.quantile(sel[,xvar], sel$wt, normwt = TRUE, prob = taus[2]/100)
  ######### (5)(6)(7) logistic regression model
      resp <- df2[, isel]
      dose <- df2[, xvar]
      lrm1 <- glm(resp ~ dose, family="binomial")
      lrm2 <- glm(resp ~ dose + I(dose^2), family="binomial")
      lrm3 <- gam(resp ~ s(dose,k=3), family= "binomial")
      if(mtype == 1) model <- lrm1
      if(mtype == 2) model <- lrm2
      if(mtype == 3) model <- lrm3

      predresp<- predict(model, newdata = data.frame(dose=xnew), type="link", se.fit=T)
      # Compute upper and lower 90% confidence limits
      up.bound.link <- predresp$fit + qnorm(0.95)*predresp$se.fit
      low.bound.link <- predresp$fit - qnorm(0.95)*predresp$se.fit
      mean.resp.link <- predresp$fit

      # Convert from logit transformed values to probability.
      up.bound <- exp(up.bound.link)/(1+exp(up.bound.link))
      low.bound <- exp(low.bound.link)/(1+exp(low.bound.link))
      mean.resp <- exp(mean.resp.link)/(1+exp(mean.resp.link))

  #### 5,6,7 full range
      ######## find modeled xc95 full data range
      lrm1.95f <- auc(xrange = xrange, mod = lrm1, dense.N = dense.N)
      lrm2.95f <- auc(xrange = xrange, mod = lrm2, dense.N = dense.N)
      lrm3.95f <- auc(xrange = xrange, mod = lrm3, dense.N = dense.N)
  #### 8,9,10 observed range
      ######## find modeled xc95 observed data range
      shortrange <- range(df2[df2[,isel]>0,xvar])
      lrm1.95p <- auc(xrange = shortrange, mod = lrm1, dense.N = dense.N)
      lrm2.95p <- auc(xrange = shortrange, mod = lrm2, dense.N = dense.N)
      lrm3.95p <- auc(xrange = shortrange, mod = lrm3, dense.N = dense.N)
  ##### 11  n values  12 quantiles
      samplen <- length(df2[df2[,isel]>0,xvar])
      limits<- quantile(df2[df2[,isel]>0,xvar],probs=taus/100)   # quantile calculation

  ###### Save results
      if(log.x) {##IF.log.x.START
      optsave[i,1:(length(taus)+11)] <- round(10^(c(WA, tol, df1[ic1,xvar],df1[ic2,xvar], eout2,
         lrm1.95f[1],lrm2.95f[1],lrm3.95f[1],lrm1.95p[1],lrm2.95p[1],lrm3.95p[1], limits)*
         diff(xlims) + xlims[1]), rounder)             # WA optimum
      } else {
      optsave[i,1:(length(taus)+11)] <- round(c(WA, tol,df1[ic1,xvar],df1[ic2,xvar],eout2,
         lrm1.95f[1],lrm2.95f[1],lrm3.95f[1],lrm1.95p[1],lrm2.95p[1],lrm3.95p[1], limits)*
         diff(xlims) + xlims[1], rounder)
      }##IF.log.x.END

      optsave[i,(length(taus)+12)] <- samplen
      #### The following code calcualte probabilities of occurrence

      binf <- cut(df2[,xvar], cutp, include.lowest = T)
      bvals <- tapply(df2[, isel] > 0, binf, mean, na.rm =T)

  # optiion for plot single species env graph
     if(plot.pdf) {##IF.plot.pdf.START
       dev.set(3)
       plot(xnew, mean.resp, axes = FALSE, type = "l", xlab = "", ylab = "", xlim=c(0,1),
                  ylim = range(low.bound, up.bound, bvals, na.rm =T), col=1)
       points(xnew, low.bound, type="l", col=1, lty="dashed")
       points(xnew, up.bound, type="l", col=1,lty="dashed")
       points(cutm, bvals, pch=21, col= 1, cex = 0.7, bg ="gray")                      # probability

  #       full <- get(paste("lrm",mtype,".95f", sep=""))
  #       part <- get(paste("lrm",mtype,".95p", sep=""))

  #       segments(part[1], part[2],part[1],0, col="red", lty="dashed")     # partial range
  #       segments(full[1], full[2],full[1],0, col="blue", lty= "dashed")   # full range

       if(lim =="GAM") { abline(v= lrm3.95f[1], lty= "dashed", col="red") ##IF.lim.START
       } else if(lim =="CDF"){
       abline(v= eout2, lty= "dashed", col="blue")
       }##IF.lim.END

       max.pow <- ceiling(max(xlims)); min.pow <- floor(min(xlims))  ## add x range for log formation

       if(log.x) {##IF.log.x.START
            at0 <- (min.pow:max.pow - xlims[1])/diff(xlims)             # add ticks at log10 = even
            lab1 <- 10^(min.pow:max.pow)
            axis(1, at = at0, labels = lab1, lwd.ticks = 1.2)    # major ticks with labels
            axis(1, at = (log10(1:10 * rep(lab1[-1]/10, each = 10)) -xlims[1])/diff(xlims), tcl = -0.3,
                labels = FALSE, lwd.ticks = 0.8)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
            xtick <- at0 <= max(df1[,xvar]) & at0 >= min(df1[,xvar])          # major labels
                       ### if only two or fewer major lables, then add tick labels at 2 and 5 log
            if(sum(xtick)<=2) {##IF.xtick.START
              axis(1, at = (log10(c(2, 5) * rep(lab1[-1]/10, each = 2)) -xlims[1])/diff(xlims), tcl = -0.4,lwd.ticks= 1.2,
                 labels = (c(2, 5) * rep(lab1[-1]/10, each = 2)))
            }##IF.xtick.END
          } else  {
            at0 <- (pretty(xlims) - xlims[1])/diff(xlims)             # add ticks
    #          lab1 <- round(at0*diff(xlims) + xlims[1], digits = rounder)    # scale back to orginial value
            axis(1, at = at0, labels = pretty(xlims), lwd.ticks = 1.5, las = 1)    # major ticks with labels
          }##IF.log.x.END

       axis(2)
  #       axis(4)
       #  mtext("Relative Abundance", side = 4, line = 2.5, col=1)
       mtext(xlabs, side = 1, line = 2.3)
       mtext("Capture Probability", side = 2, line = 2.3,col=1)
       mtext(bquote(italic(.(tnames[i]))), side = 3, line = 0.5)
       box()

       if((i+5)/6-round((i+5)/6)==0) {##IF.i+5.START
          title(main= main, outer=TRUE, cex.main= 1.8)
          mtext(paste("Page", (i+5)/6, "of", totpage), outer = TRUE, side =1)
        }##IF.i+5.END
  #### plot cdf
       dev.set(2)
       Ecdf(sel[,xvar], weights = sel$wt, xlim = c(0,1), col = "blue", pch = 1,
            axes = F, main = tnames[i], xlab = xlabs)
  #       plot.ecdf(sel[,xvar], xlim = c(0,1), col = "green", pch = 1,add = T,
  #            axes = F, main = tnames[i], xlab = xlabs)
            lab1 <- round(axTicks(1)*diff(xlims) + xlims[1], digits = rounder)    # scale back to orginial value
            axis(1, at = axTicks(1), labels = lab1, lwd.ticks = 1.5, las = 1)    # major ticks with labels
       axis(2)
       abline(h = c(taus[2]/100), col = "red", lty = 2)
       abline(v = eout2, col = "red", lty = 2)
       box(bty = "l")
    }##IF.plot.pdf.END
   }##FOR.i.END
    optout <- data.frame(tnames, optsave)
    names(optout) <- c("taxaname", "WAopt","WAtol", "CDF_ABD","CDF_PA", "CDF_WT_QTILE", "LRM_Full",
       "QLRM_Full","GAM_Full", "LRM_Ob", "QLRM_Ob", "GAM_Ob",
        paste("Observ",taus, sep=""),"N")

    if(plot.pdf) graphics.off()

    return(optout)

}##FUNCTION.taxon.response.END

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## inverse deshrinking function
## `inv.deshrink` <- function(env, wa.env) {
##     X <- cbind(rep(1, length(wa.env)), wa.env)
##     QR <- qr(X)
##     coef <- qr.coef(QR, env)
##     pred <- qr.fitted(QR, env)
##     return(list(coefficients = coef, env = pred))
## }
##
## classical deshrinking
## `class.deshrink` <- function(env, wa.env) {
##     X <- cbind(rep(1, length(env)), env)
##     QR <- qr(X)
##     coef <- drop(qr.coef(QR, wa.env))
##     coef <- c(-coef[1], 1)/coef[2]
##     pred <- deshrink.pred(wa.env, coef)
##     return(list(coefficients = coef, env = pred))
## }
##
## deshrinking to equal sd
## A bit like in vegan:::wascores, but wascores uses weighted sd which
## would need row and column sums in the function call, and this would
## make the function API incompatible with other *.deshrink functions.
## `expand.deshrink` <- function(env, wa.env) {
##     b1 <- sd(env)/sd(wa.env)
##     b0 <- mean(env) - b1 * mean(wa.env)
##     pred <- b0 + b1 * wa.env
##     return(list(coefficients = c(b0, b1), env = pred))
## }
##
# Do not deshrink: for those who think they know what they do
## `no.deshrink` <- function(env, wa.env) {
##     return(list(coefficients = c(0, 1), env = wa.env))
## }
