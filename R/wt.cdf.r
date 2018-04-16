  ##~~~~~~~~~~~~~~  4)  fish condutivity weighted method for appalanchian region
  # Chem.eco3 ; fish
##
##  weighted.cdf function to calculate XC95 values
##  data.env  environmental data
##  data.sp.wide, species crosstabed data
##  plot, a boolean to choose if plot cdf and gam plots
##  dogam, a booleen to choose if a gam fit is calculated
##  sortvect: to provide a vect of species list so plots will be sorted according to the list
##  nt: minimum number of occurence
##  addtrend: if a trend should be added ( = ">" etc) in the output
##  np: number of bins
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#' @title A function to calculate a weighted CDF for XC95 values.
#'
#' @description Weighted cdf function to calculate XC95 values on log10 transformed data.
#' 
#' @details Grouping variable requires three inputs.  The first two are 
#' categorical and the third is continuous.  For example, HUC, Ecoregion, and 
#' Watershed Area.  "groups" specifies the column names and "groupNamesOut" species
#' the output names.  The function will stop if groups is left as null.  
#' The function assumes log10 transformed data.
#' 
#' This function replaces fish.wt.cdf.  This function is more generic than fish.wt.cdf.
#'
#' @param data.env environmental data, default = "datafile" from global environment.
#' @param data.sp.wide Species wide (crosstabed) data; default = "ss" from global environment.
#' @param plot A boolean to choose if plot cdf and gam plots; default = F.
#' @param dogam A booleen to choose if a gam fit is calculated; default = F.
#' @param SampleID Site/sample id column used in both data files; default = "Station_Date"
#' @param tag Default = "".
#' @param sortvect to provide a vector of species list so plots will be sorted according to the list; default = NULL.
#' @param np Number of bins; default = 61.
#' @param nt Minimum number of occurence; default = 25.
#' @param addtrend A booleen if a trend should be added ( = ">" etc) in the output (T or F); default = F.
#' @param wd Working directory for saving files.
#' @param groups Three grouping variables (first 2 categorical and third is continuous).  Column names in data.env, e.g., HUC (BigHUC), Ecoregion (ECOREGL3), and Watershed Area (WS_Area);  default = NULL.
#' @param xvar variable on which to base calculations; default  = "cond"
#' @param groupNamesOut Group names in output.  Keep short. 
#' @param xvar.PlotName Plot name for xvar; default = paste0("Conductivity ( ", mu, "S/cm)").
#'
#' @return A dataframe of XC95 values and "gam" and "cdf" subfolders of "wd" with TIFF files of plots.
#'  
#' @examples
#' # data
#' data(dta.do)
#' data(ss.sites) 
#' # function inputs
#' data.env     <- dta.do
#' data.sp.wide <- ss.sites
#' plot         <- T
#' dogam        <- F
#' SampleID     <- "Station_Date"
#' tag          <- "wt"
#' sortvect     <- NULL
#' np           <- 61
#' nt           <- 25
#' addtrend     <- F
#' wd           <- getwd()
#' groups       <- c("BigHUC","ECOREGL3","WS_AREA")
#' xvar         <- "cond"
#' groupNames   <- c("HUC04", "EcoL3", "Area")
#' xvar.PlotName <- "Conductivity (uS/cm)"
#' # run function (~20 seconds)
#' df.xc95 <- wt.cdf (data.env, data.sp.wide, plot = T, dogam = T
#'                  , SampleID = SampleID, tag = tag, sortvect = NULL
#'                  , np = 61, nt = 25, addtrend = T, wd = getwd()
#'                  , groups = groups, xvar = xvar, groupNames = groupNames
#'                  , xvar.PlotName = xvar.PlotName) 
#' View(df.xc95)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC
#
# library(XC95)
# data.env     <- dta.do
# data.sp.wide <- ss.sites
# plot         <- T
# dogam        <- F
# SampleID     <- "Station_Date"
# tag          <- "wt"
# sortvect     <- NULL
# np           <- 61
# nt           <- 25
# addtrend     <- F
# wd           <- getwd()
# groups       <- c("BigHUC","ECOREGL3","WS_AREA")
# xvar         <- "cond"
# groupNames   <- c("HUC04", "EcoL3", "Area")
# xvar.PlotName <- "Conductivity (uS/cm)"
# #
# #names(data.sp.wide)[1] <- SampleID
# df.xc95 <- wt.cdf (data.env, data.sp.wide, plot = T, dogam = T
#                  , SampleID = SampleID, tag = tag, sortvect = NULL
#                  , np = 61, nt = 25, addtrend = T, wd = getwd()
#                  , groups = groups, xvar = xvar, groupNames = groupNames
#                  , xvar.PlotName = xvar.PlotName) 
# View(df.xc95)
# 
# # new data
# data.env     <- read.csv("env.bio69.csv")
# data.sp.wide <- read.csv("ss.csv")
# SampleID      <- "Sample.ID"
# groups       <- c("Level4", "Level3", "ACRES")
# xvar         <- "Conductivity"
# #
# data.env[,xvar] <- log10(data.env[,xvar])
# 
# #
# index <- 1
#
#' @export
wt.cdf <- function(data.env, data.sp.wide, plot = T, dogam = F
                   , SampleID = "Station_Date", tag = "", sortvect = NULL
                   , np = 61, nt = 25, addtrend = F, wd = getwd()
                   , groups = NULL, xvar = "cond", groupNames=NULL
                   , xvar.PlotName = "Specific Conductivity ( uS/cm)") 
{##FUNCTION.fish.wt.cdf.START
  #
  # 20170424, change df1$cond to df1[,xvar]
  #
  # 20170424, add directory check
  # 20180411, Remove region, alreadying using wd in function call.
  #region <- "Results"
  #dir.check.add(wd, region)
  # dir.check.add(file.path(wd,region),"cdf")
  # dir.check.add(file.path(wd,region),"gam")
  dir.check.add(wd, "cdf")
  dir.check.add(wd, "gam")
  #wd <- file.path(wd,region)
  
  # QC, Fields ####
  # 20180412
  if((SampleID %in% names(data.env))==FALSE){
    stop("SampleID missing from data.env.")
  }
  if((SampleID %in% names(data.sp.wide))==FALSE){
    stop("SampleID missing from data.sp.wide.")
  }
  if((xvar %in% names(data.env))==FALSE){
    stop("xvar missing from data.env.")
  }
  if(sum(groups %in% names(data.env))!=length(groups)){
    stop("One or more groups is missing from data.env.")
  }
  if(length(groups)!=length(groupNames)){
    stop("Number of groups and groupNames does not match.")
  }
  if(is.null(groups)){
    stop("At least one group *must* be specified.")
  }
  
  # Merge ####
  #df.sp.merge <- merge(data.env[c(SampleID, "HUC", "cond", "BigHUC","ECOREGL3","WS_AREA")], data.sp.wide); dim(df.sp.merge)   # merge env data with species data
  # 20170418
  #df.sp.merge <- merge(data.env, data.sp.wide); dim(df.sp.merge)
  # 20180412
  col.Keep <- c(SampleID, xvar, groups)
  df.sp.merge <- merge(data.env[, col.Keep], data.sp.wide, by=SampleID)

  # 20170418, size of data.env
  #tcol.start <- ncol(data.env)+1
  # 20180412
  tcol.start <- length(col.Keep) + 1
  
  # preliminary count of number of occurences
  tnames.count <- apply(df.sp.merge[tcol.start:ncol(df.sp.merge)] > 0, 2, sum)  
  tnames.sav <- tnames.count[tnames.count >= nt ]; length(tnames.sav)
  if( is.null(sortvect)) {##IF.sortvect.START
    mod.tnames <- names(tnames.sav)
  } else mod.tnames <- sortvect
  ##IF.sortvect.END

  tolval.cdf <- rep(NA, times = length(mod.tnames))  # xc95 from cdf
  tolval.gam <- rep(NA, times = length(mod.tnames))  # xc95 from gam
  total.n    <- rep(NA, times = length(mod.tnames))  # Xc95 for number of occurence sites
  total.N    <- rep(NA, times = length(mod.tnames))  # xc95 for total number of sites
  trend      <- rep(NA, times = length(mod.tnames))  ##  > or ~ or <
  # wtshd.df   <- rep(NA, times = length(mod.tnames))
  # eco3.df    <- rep(NA, times = length(mod.tnames))
  # huc.df     <- rep(NA, times = length(mod.tnames))
  Grp1.df   <- rep(NA, times = length(mod.tnames))  # HUC
  Grp2.df   <- rep(NA, times = length(mod.tnames))  # Ecoregion
  Grp3.df   <- rep(NA, times = length(mod.tnames))  # Watershed Area
  
  # Groups, 20180412
  # for (i in groups){##FOR.i.START
  #   assign(paste0(i,".df"), rep(NA, times = length(mod.tnames)))
  #   i.num <- match(i, groups)
  #   assign(paste0("groups",i.num), i)
  # }##FOR.i.END
  
  


#  row.names(wtshd.df) <- mod.tnames
#  if(plotcdf) {        # device setup



#    pdf(file = paste(wd, "/Results/",tag,"cdf.fish.wtshed.pdf", sep = ""), width = 9, height = 6, pointsize = 11)
#    par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#  } else if (plotgam) {
#    pdf(file = paste(wd, "/Results/",tag,"gam.fish.wtshed.pdf", sep = ""), width = 9, height = 6, pointsize = 11)
#    par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#  }

  # # parameter-ize (EWL, 20170424)
  # groups1 <- groups[1] # "BigHUC"
  # groups2 <- groups[2] # "ECOREGL3"
  # groups3 <- groups[3] # "WS_AREA"
  # # 20180412, plan for any number of groups (see code above)


  for (index in 1:length(mod.tnames)) {##IF.index.START
    ## for mod.taxa list
    imatch <- match(mod.tnames[index], names(df.sp.merge)); imatch      # find colum index
    #HUCs <- unique(df.sp.merge[df.sp.merge[,imatch] > 0 ,"BigHUC"])   # find all Site HUCs
    Grp1 <- unique(df.sp.merge[df.sp.merge[,imatch] > 0 ,groups[1]])   # find all Site Grp1
    Grp1.list <- paste(sort(Grp1), collapse ="_")
    #eco3 <- paste(unique(df.sp.merge[df.sp.merge[,imatch]>0,"ECOREGL3"]), collapse="_")
    Grp2 <- paste(sort(unique(df.sp.merge[df.sp.merge[,imatch]>0, groups[2]])), collapse="_")
    #wtshed <- paste(range(df.sp.merge[df.sp.merge[,imatch]>0,"WS_AREA"],na.rm =T), collapse="_")
    Grp3 <- paste(range(df.sp.merge[df.sp.merge[,imatch]>0,groups[3]],na.rm =T), collapse="_")
#    df1 <- subset(df.sp.merge, BigHUC %in% HUCs); dim(df1)
    #20180412, subset using Group1
    df1 <- df.sp.merge[df.sp.merge[,groups[1]] %in% Grp1, ]
    #print(paste(nrow(df.sp.merge), nrow(df1)))
    
    # groups.1.items <- unique(df.sp.merge[df.sp.merge[,imatch] > 0 ,groups[1]])
    # groups.1.sort  <- paste(sort(groups.1.items), collapse ="_")
    # groups.2.items <- unique(df.sp.merge[df.sp.merge[,imatch] > 0 ,groups[2]])
    # groups.2.sort  <- paste(sort(groups.1.items), collapse ="_")

    
    # cnew <- seq(from = min(df1$cond), to = max(df1$cond), length = 100)      # for plot
    # new.data <- data.frame(cond = cnew)
    # cond.u <- sort(unique(df1$cond))
    cnew <- seq(from = min(df1[,xvar]), to = max(df1[,xvar]), length = 100)      # for plot
    new.data <- data.frame(xvar = cnew)
    # 20180413
    names(new.data) <- xvar
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
    # eco3.df[index] <- eco3
    # wtshd.df[index] <- wtshed
    # huc.df[index] <- huc.list
    Grp1.df[index] <- Grp1.list
    Grp2.df[index] <- Grp2
    Grp3.df[index] <- Grp3
    
    #wtshd.df[index,] <- c(mywtshed[1], mywtshed[2], huc.list)
    # Do GAM ####
    if (dogam) {##IF.dogam.START
      # 20180413, "fix" code by changing xvar to "cond"
      names(df1)[match(xvar, names(df1))] <- "xvar"
      names(new.data) <- "xvar"
      ### only if dogam is selected, calculate gam based xc95 and plot
      mod <- mgcv::gam(resp ~ s(xvar, k = 3), data = df1, family = "binomial")
      predresp <- predict(mod, new.data, type = "link", se.fit = T)
      # Compute upper and lower 90% confidence limits
      up.bound.link <- predresp$fit + qnorm(0.95) * predresp$se.fit
      low.bound.link <- predresp$fit - qnorm(0.95) * predresp$se.fit

      # Convert from logit transformed values to probability.
      up.bound <- exp(up.bound.link)/(1 + exp(up.bound.link))
      low.bound <- exp(low.bound.link)/(1 + exp(low.bound.link))

      predresp.vline <- predict(mod, newdata = data.frame(xvar = max(df1[, "xvar"], na.rm=T)), type = "link", se.fit = T)
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
      # 20180413, change xvar back
      names(df1)[match("xvar", names(df1))] <- xvar
      names(new.data) <- xvar
    }##IF.dogam.END
    
    # Plot ####
    if(plot) {##IF.plot.START
       if(index%%6==1 ) {##IF.index.START
       pgs <- (index-1)%/%6 + 1
       tiff(file = paste(wd,"/cdf/", pgs, ".taxon.cdf.tiff",sep=""),
           width = 650, height = 450, pointsize = 13)
       par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))

       tiff(file = paste(wd,"/gam/", pgs, ".taxon.gam.tiff",sep=""),
           width = 650, height = 450, pointsize = 13)
        par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
        }##IF.index.END
     dev.set(2)
     Hmisc::Ecdf(df2[,xvar], weights = df2$wt, ylim = c(0,1), col = "blue", pch = 1, axes = F,
           main = bquote(italic(.(mod.tnames[index]))), ylab = "Proportion of occurrence ",
           #xlab = expression(paste("Conductivity ( ", mu, "S/cm)")))
           xlab = xvar.PlotName)
      abline(h = 0.95, col = "red", lty = 2)
      abline(v = tolval.cdf[index], col = "red", lty = 2)
      box(bty = "l")
      min.pow <- floor(min(df1[,xvar])); max.pow <- ceiling(max(df1[,xvar]))
      at0 <- min.pow:max.pow              # add ticks at log10 = even
      lab0 <- 10^(min.pow:max.pow)
      axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
      axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
           labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
      xtick <- at0 <= max(df1[,xvar]) & at0 >= min(df1[,xvar])          # major labels
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
      lines(new.data[,xvar], plogis(up.bound.link), lty = 2)
      lines(new.data[,xvar], plogis(low.bound.link), lty = 2)
      lines(new.data[,xvar], mn)
      min.pow <- floor(min(df1[,xvar])); max.pow <- ceiling(max(df1[,xvar]))
      at0 <- min.pow:max.pow              # add ticks at log10 = even
      lab0 <- 10^(min.pow:max.pow)
      axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
      axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
           labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
      xtick <- at0 <= max(df1[,xvar]) & at0 >= min(df1[,xvar])          # major labels
      ### if only two or fewer major lables, then add tick labels at 2 and 5 log
      if(sum(xtick)<=2) {
        axis(1, at = log10(c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 1,
             labels = (c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))))
      }
      axis(2)
      #mtext(expression(paste("Specific conductivity ( ", mu, "S/cm)")), side = 1, line = 2.3, cex = 0.8)
      mtext(xvar.PlotName, side = 1, line = 2.3, cex = 0.8)
      mtext("Probability of observing", side = 2, line = 2.3, cex =0.8)
      mtext(bquote(italic(.(mod.tnames[index]))), side = 3.5)
      box()
      abline(v = tolval.cdf[index], col = "red", lty = 2)

      if (index%%6==0 |index== length(mod.tnames)) graphics.off()
    }##IF.plot.END
  }##IF.index.START
  #
  # Add Trend ####
  # dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf), trend = trend,
  #                    XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
  #                    HUC = Grp1.df, eco3 = Grp2.df, wtshed = Grp3.df, stringsAsFactors = F)
  # 20180413, remove 10^ for variable
  dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf), trend = trend,
                     XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
                     HUC = Grp1.df, eco3 = Grp2.df, wtshed = Grp3.df, stringsAsFactors = F)
  # rename Group columns to thos provided by user
  names(dftv)[7:9] <- groupNames
  if(addtrend) {##IF.addtrend.START
    # dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf), trend = trend,
    #                    XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
    #                    HUC = huc.df, eco3 = eco3.df, wtshed = wtshd.df, stringsAsFactors = F)
  } else {
    # dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf),
    #                    XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
    #                    HUC = huc.df, eco3 = eco3.df, wtshed = wtshd.df, stringsAsFactors = F)
    dftv <- dftv[,-3]  # drop "trend"
  }##IF.addtrend.END
  dftv <- dftv[order(dftv$XC95.cdf),]
#  write.table(wtshd.df,paste(wd, "/Results/wtshed.size.csv", sep= ""),sep=",", row.names = T)
  return(dftv)                     # return a data.frame with xc95s
}##FUNCTION.fish.wt.cdf.END
