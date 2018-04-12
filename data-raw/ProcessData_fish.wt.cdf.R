# Prepare data for example for fish.wt.cdf.R; fish.weight.cdf()
#
# From _Figure_G12_G13_G14_G15.R
# named weighted.fish.cdf() in this file
#
# Erik.Leppo@tetratech.com
# 20170424
##~~~~~~~~~~~~~~~~~~~~

# 0. Prep####
# library(devtools)
# library(reshape)
wd <- getwd()

# 1. Get data and process#####

# Figure G-12, G-13, G-14, and G-15 code only
#

# clear the workspace
rm(list=ls())

# libraries
# library(reshape)
# library(Hmisc)
# library(mgcv)

# get working directory, using TINN-R or R-Studio so opening script sets the wd
wd <- getwd()

# Source helper scripts
#source(paste(wd,"/R_Function/weightcdf2.r", sep =""))
#source(paste(wd,"/R_Function/curve.shape.r", sep =""))




#Fish <- read.csv(paste(wd, "/Data/Combined_Less.csv", sep = ""), header = T); dim(Fish)   # 42336
#Fish <- read.csv(paste(wd, "/AppendixG/Data/Combined_Less.csv", sep = ""), header = T); dim(Fish)   # 42336
Fish <- read.csv(file.path(wd,"data-raw","Combined_Less.csv"), header = T); dim(Fish)   # 42336
Fish$RA[is.na(Fish$RA)] <- 1
#  length(unique(Fish$SITE_ID))
#   3278
#  length(unique(Fish$Station_Date))
#   3789
#  length(unique(Fish$Station_Year))
#   3559
# station year june 609230_40401
Chem <- subset(Fish, !duplicated(Station_Date));dim(Chem)  # 3789
Chem$cond <- Chem$lcond
Chem <- subset(Chem, (cond>log10(5))) ; dim(Chem) # 3788
Chem <- Chem[order(Chem$YEAR, Chem$month),]
Chem.ecosites <- subset(Chem, !duplicated(SITE_ID, fromLast=T)) ; dim(Chem.ecosites) # 3277
#write.table(Chem.ecosites, paste(wd, "/Results/chem.csv", sep =""), sep =",", row.names =F)

summary(Chem$COND  )
species <- subset(Fish, Station_Date %in% Chem.ecosites$Station_Date); dim(species)
ss.sites <-cast(species, Station_Date ~ species, sum, value = "RA"); dim(ss.sites) # 3277 213
ss.genera <- cast(species, Station_Date ~ GENUS, sum, value = "RA"); dim(ss.genera) # 3277 70
#write.table(species, paste(wd, "/Results/species_all.csv", sep =""),sep=",", row.names =F)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SSD ~~~~~~~~~~~~~~~~~~~~~~



# allsites <-   weightcdf(df1 = Chem.ecosites, ss = ss.sites, SampID = "Station_Date", xvar = "cond",
#                         nt = 25)  ; dim(allsites) # 101
# allsites <- allsites[order(allsites$XC95.all),]
# hc05.site <- apply(allsites[-1], 2, quantile, prob =0.05, type =6, na.rm = T)[1]      # 512
#
# allgenera <-   weightcdf(df1 = Chem.ecosites, ss = ss.genera, SampID = "Station_Date", xvar = "cond",
#                          nt = 25)  ; dim(allgenera)
# allgenera <- allgenera[order(allgenera$XC95.all),]
# hc05.genera <- apply(allgenera[-1], 2, quantile, prob =0.05, type =6, na.rm = T)[1]  # 546

#~~~~~~~~~~~~~~~~~~~~~~~~~~ 4)  fish condutivity weighted method ~~~~~~~~~~~~~~~~~~~~~~~~~~
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
# weighted.cdf <- function(datafile = datafile, ss = ss, plot=F, dogam = F,
#                          SampleID = "Station_Date",tag="", sortvect = NULL, np = 61, nt = 25, addtrend =F) {
#
#   my.ss <- merge(datafile[c(SampleID, "HUC", "cond", "BigHUC","ECOREGL3","WS_AREA")], ss); dim(my.ss)   # merge env data with species data
#
#   tnames.count <- apply(my.ss[7:ncol(my.ss)] > 0, 2, sum)            # preliminary count of number of occurences
#   tnames.sav <- tnames.count[tnames.count >= nt ]; length(tnames.sav)
#   if( is.null(sortvect)) {
#     mod.tnames <- names(tnames.sav)
#   } else mod.tnames <- sortvect
#
#   tolval.cdf <- rep(NA, times = length(mod.tnames))  # xc95 from cdf
#   tolval.gam<- rep(NA, times = length(mod.tnames))   # xc95 from gam
#   total.n <- rep(NA, times = length(mod.tnames))    # Xc95 for number of occurence sites
#   total.N <- rep(NA, times = length(mod.tnames))    # xc95 for total number of sites
#   trend <- rep(NA, times = length(mod.tnames))      ##  > or ~ or <
#   wtshd.df <- rep(NA, times = length(mod.tnames))
#   eco3.df <- rep(NA,times = length(mod.tnames))
#   huc.df <- rep(NA,times = length(mod.tnames))
#
#   #  row.names(wtshd.df) <- mod.tnames
#   #  if(plotcdf) {        # device setup
#
#
#
#   #    pdf(file = paste(wd, "/Results/",tag,"cdf.fish.wtshed.pdf", sep = ""), width = 9, height = 6, pointsize = 11)
#   #    par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#   #  } else if (plotgam) {
#   #    pdf(file = paste(wd, "/Results/",tag,"gam.fish.wtshed.pdf", sep = ""), width = 9, height = 6, pointsize = 11)
#   #    par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#   #  }
#
#   for (index in 1:length(mod.tnames)) {  ## for mod.taxa list
#     imatch <- match(mod.tnames[index], names(my.ss)); imatch      # find colum index
#     HUCs <- unique(my.ss[my.ss[,imatch] > 0 ,"BigHUC"])   # find all Site HUCs
#     huc.list <- paste(sort(HUCs), collapse ="_")
#     eco3 <- paste(unique(my.ss[my.ss[,imatch]>0,"ECOREGL3"]), collapse="_")
#
#     wtshed <- paste(range(my.ss[my.ss[,imatch]>0,"WS_AREA"],na.rm =T), collapse="_")
#     df1 <- subset(my.ss, BigHUC %in% HUCs); dim(df1)
#
#
#     cnew <- seq(from = min(df1$cond), to = max(df1$cond), length = 100)      # for plot
#     new.data <- data.frame(cond = cnew)
#     cond.u <- sort(unique(df1$cond))
#
#     cutp <- seq(from = min(df1$cond), to = max(df1$cond), length = np)  # cut entire gradient into np bins
#     cutm <- 0.5*(cutp[-1] + cutp[-np])                   # find middle point
#     df1$cutf <- cut(df1$cond, cutp, include.lowest = T)   # define bins
#
#     wt <- 1/table(df1$cutf)                 # weight for each bin
#     wt[is.infinite(wt)] <- NA
#     wt <- wt/sum(wt, na.rm = T)
#
#     dftemp <- data.frame(cutf = names(wt), wt = as.vector(wt))
#     df1 <- merge(df1, dftemp, by = "cutf"); dim(df1); names(df1)  # final data
#
#     #  ECDF of observed conductivity
#     resp <- df1[, mod.tnames[index]] > 0              # response variable
#     df2 <- df1[resp,]
#
#     tolval.cdf[index] <- wtd.quantile(df2$cond, df2$wt, normwt = TRUE, prob = 0.95) # xc95 calculation
#     total.n[index] <- nrow(df2)
#     total.N[index] <- nrow(df1)
#     eco3.df[index] <- eco3
#     wtshd.df[index] <- wtshed
#     huc.df[index] <- huc.list
#     #   wtshd.df[index,] <- c(mywtshed[1], mywtshed[2], huc.list)
#     if (dogam) {  ### only if dogam is selected, calculate gam based xc95 and plot
#       mod <- gam(resp ~ s(cond, k = 3), data = df1, family = "binomial")
#       predresp <- predict(mod, new.data, type = "link", se.fit = T)
#       # Compute upper and lower 90% confidence limits
#       up.bound.link <- predresp$fit + qnorm(0.95) * predresp$se.fit
#       low.bound.link <- predresp$fit - qnorm(0.95) * predresp$se.fit
#
#       # Convert from logit transformed values to probability.
#       up.bound <- exp(up.bound.link)/(1 + exp(up.bound.link))
#       low.bound <- exp(low.bound.link)/(1 + exp(low.bound.link))
#
#       predresp.vline <- predict(mod, newdata = data.frame(cond = max(df1$cond, na.rm=T)), type = "link", se.fit = T)
#       low500.lk <- predresp.vline$fit - qnorm(0.95) * predresp.vline$se.fit
#       low500.resp <- exp(low500.lk)/(1 + exp(low500.lk))
#       mean500.lk <- predresp.vline$fit
#       mean500.resp <- exp(mean500.lk)/(1 + exp(mean500.lk))
#       mn <- plogis(predresp$fit)
#       crit <- 0.05 * max(mn)           # arbitrarily define a criterion to assign trend
#       if(mean500.resp <= crit) {
#         trend[index] <- "="
#       }  else  if(low500.resp <= crit & mean500.resp> crit ) {
#         trend[index] <- "~"
#       }  else {
#         trend[index] <- ">"
#       }
#
#
#       bararea <- (cnew[2] - cnew[1])* mn         ### determine gam xc95 values
#       tot <- sum(bararea)
#       ii <- 1
#       while(sum(bararea[1:ii]) < 0.95 * tot) ii <- ii + 1
#       tolval.gam[index] <- cnew[ii]
#     }
#
#     if(plot) {
#       if(index%%6==1 ) {
#         pgs <- (index-1)%/%6 + 1
#         tiff(file = paste(wd,"/Results/cdf/", pgs, ".taxon.cdf.tiff",sep=""),
#              width = 650, height = 450, pointsize = 13)
#         par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#
#         tiff(file = paste(wd,"/Results/gam/", pgs, ".taxon.gam.tiff",sep=""),
#              width = 650, height = 450, pointsize = 13)
#         par(mfrow = c(2, 3), pty = "m", mar = c(4, 5, 3, 2))
#       }
#       dev.set(2)
#       Ecdf(df2[,"cond"], weights = df2$wt, ylim = c(0,1), col = "blue", pch = 1, axes = F,
#            main = bquote(italic(.(mod.tnames[index]))), ylab = "Proportion of occurrence ",
#            xlab = expression(paste("Conductivity ( ", mu, "S/cm)")))
#       abline(h = 0.95, col = "red", lty = 2)
#       abline(v = tolval.cdf[index], col = "red", lty = 2)
#       box(bty = "l")
#       min.pow <- floor(min(df1$cond)); max.pow <- ceiling(max(df1$cond))
#       at0 <- min.pow:max.pow              # add ticks at log10 = even
#       lab0 <- 10^(min.pow:max.pow)
#       axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
#       axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
#            labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
#       xtick <- at0 <= max(df1[,"cond"]) & at0 >= min(df1[,"cond"])          # major labels
#       ### if only two or fewer major lables, then add tick labels at 2 and 5 log
#       if(sum(xtick)<=2) {
#         axis(1, at = log10(c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 1,
#              labels = (c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))))
#       }
#       axis(2)
#
#       dev.set(3)
#       val <- tapply(resp, df1$cutf, mean)
#       plot(cutm, val, type = "n", ylim = range(c(val, plogis(up.bound.link), plogis(low.bound.link)),na.rm = T), axes= F,
#            xlab = "" , ylab = "")
#       points(cutm, val)
#       lines(new.data$cond, plogis(up.bound.link), lty = 2)
#       lines(new.data$cond, plogis(low.bound.link), lty = 2)
#       lines(new.data$cond, mn)
#       min.pow <- floor(min(df1$cond)); max.pow <- ceiling(max(df1$cond))
#       at0 <- min.pow:max.pow              # add ticks at log10 = even
#       lab0 <- 10^(min.pow:max.pow)
#       axis(1, at = at0, labels = lab0, lwd.ticks = 1)    # major ticks with labels
#       axis(1, at = log10(1:10 * rep(lab0[-1]/10, each = 10)) , tcl = -0.3,
#            labels = FALSE, lwd.ticks = 0.9)  # minor ticks            xtick <- at0 <= max(x) & axis.at >= xmin          # major labels
#       xtick <- at0 <= max(df1[,"cond"]) & at0 >= min(df1[,"cond"])          # major labels
#       ### if only two or fewer major lables, then add tick labels at 2 and 5 log
#       if(sum(xtick)<=2) {
#         axis(1, at = log10(c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))), tcl = -0.4,lwd.ticks= 1,
#              labels = (c(2, 5) * rep(lab0[-1]/10, each = sum(xtick))))
#       }
#       axis(2)
#       mtext(expression(paste("Specific conductivity ( ", mu, "S/cm)")), side = 1, line = 2.3, cex = 0.8)
#       mtext("Probability of observing", side = 2, line = 2.3, cex =0.8)
#       mtext(bquote(italic(.(mod.tnames[index]))), side = 3.5)
#       box()
#       abline(v = tolval.cdf[index], col = "red", lty = 2)
#
#       if (index%%6==0 |index== length(mod.tnames)) graphics.off()
#     }
#   }
#
#   if(addtrend) {
#     dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf), trend = trend,
#                        XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
#                        HUC = huc.df, eco3 = eco3.df, wtshed = wtshd.df, stringsAsFactors = F)
#   } else {
#     dftv <- data.frame(GENUS = mod.tnames, XC95.cdf = round(10^tolval.cdf),
#                        XC95.gam = round(10^tolval.gam), total.n = total.n, total.N = total.N,
#                        HUC = huc.df, eco3 = eco3.df, wtshed = wtshd.df, stringsAsFactors = F)
#   }
#   dftv <- dftv[order(dftv$XC95.cdf),]
#   #  write.table(wtshd.df,paste(wd, "/Results/wtshed.size.csv", sep= ""),sep=",", row.names = T)
#   return(dftv)                     # return a data.frame with xc95s
# }
# ###################################  Chapter 2  calcualte hc05 for different conditions
# t1 <- proc.time()
# taxalist2 <- as.vector(allsites$taxaname)
# dftv2 <- weighted.cdf(datafile = Chem.ecosites, ss=ss.sites, SampleID="Station_Date", tag="sites",
#                       plot=F, dogam = F,  addtrend =F, sortvect =taxalist2)
#
# taxalist2 <- as.vector(dftv2$GENUS[order(dftv2$XC95.cdf)])
#
# dftv2 <- weighted.cdf(datafile = Chem.ecosites, ss=ss.sites, SampleID="Station_Date", tag="sites",
#                       #     plot=T, dogam = T,  addtrend =T, sortvect =taxalist2)
#                       plot=F, dogam = T,  addtrend =T, sortvect =taxalist2)
# (hcs0 <- apply(dftv2[, c(2,4)], 2, quantile, prob = 0.05, type=6, na.rm = T))
# #write.table(dftv2, paste(wd, "/Results/xc95_allfish.csv", sep =""), sep =",", row.names =F)
#
#
# taxalist3 <- as.vector(allgenera$taxaname)
#
# dftv3 <- weighted.cdf(datafile = Chem.ecosites, ss=ss.genera, SampleID="Station_Date", tag="genera",
#                       plot =F, dogam = T,
#                       addtrend =T, sortvect =taxalist3)
#
# dftv3 <- merge(allgenera, dftv3,  by.x ="taxaname", by.y = "GENUS", all=T)
#
# names(dftv3)[c(2,4)] <-c("XC95.nowt","XC95.wt")
#
# (hcs1 <- apply(dftv3[, c(2,4)], 2, quantile, prob = 0.05, type=6, na.rm = T))
#
# #write.table(dftv3, paste(wd, "/Results/xc95_allgenera.csv", sep =""), sep =",", row.names =F)


###############################confounder analysis

# dta.rbp <- subset(Chem.ecosites, rbpscore>135|is.na(rbpscore)); dim(dta.rbp)
# dftv.rbp <- weighted.cdf(datafile = dta.rbp, ss=ss.sites,plot = F,
#                          dogam = F) ; dim(dftv.rbp)
# (hcs2 <- apply(dftv.rbp[, 2:5], 2, quantile, prob = 0.05, type=6, na.rm = T)) #  432

dta.do <- subset(Chem.ecosites, do > 4|is.na(do)); dim(dta.do)
dftv.do <- weighted.cdf(datafile = dta.do, ss=ss.sites,plot= F,
                        dogam = F)  ; dim(dftv.do)
(hcs3 <- apply(dftv.do[, 2:5], 2, quantile, prob = 0.05, type=6, na.rm = T))     #392

# dta.tmp <- subset(Chem.ecosites, Temp <22 | is.na(Temp)); dim(dta.tmp)
# dftv.tmp <- weighted.cdf(datafile = dta.tmp, ss=ss.sites, plot = F,
#                          dogam = F) ;  dim(dftv.tmp)
# (hcs4 <- apply(dftv.tmp[, 2:5], 2, quantile, prob = 0.05, type=6, na.rm = T))  #419
#
# dta.wtshd <- subset(Chem.ecosites, WS_AREA>10 | is.na(WS_AREA)); dim(dta.wtshd)
# dftv.wtshd <- weighted.cdf(datafile = dta.wtshd, ss=ss.sites, plot= F,
#                            dogam = F)  ; dim(dftv.wtshd)
# (hcs5 <- apply(dftv.wtshd[, 2:5], 2, quantile, prob = 0.05, type=6, na.rm = T))   #274
#
# #~~~~~~~~~~~~~~~~~~~~~~~~~~
# ## Figure G12
# #jpeg(filename = paste(wd, "/Results/Figure_G_confounder_rbp.jpg",sep=""), width = 400, height = 400,
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/AppendixG/Figure_G12_confounder_rbp.jpg",sep=""), width = 400, height = 400,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# dftv.1 <- merge(dftv2[1:2], dftv.rbp[1:2], by = "GENUS", all = T)
# dftv.2 <- merge(dftv2[1:2], dftv.tmp[1:2], by = "GENUS", all = T)
# dftv.3 <- merge(dftv2[1:2], dftv.do[1:2], by = "GENUS", all = T)
# dftv.4 <- merge(dftv2[1:2], dftv.wtshd[1:2], by = "GENUS", all = T)
# #  par(mfrow = c(2,2))
# par(mar = c(5,4,1,1) )
# plot.ecdf(dftv.1$XC95.cdf.x,pch=1, ylab="Proportion of species", col = 1, y = c(0,0.5),
#           xlim=range(100, 3200), xlab="Specific conductivity (?S/cm)",
#           axes=F, log="x", main="") #Excluding rbpscore < 135", font.main = 1)                              # quantile xc95
# plot.ecdf(dftv.1$XC95.cdf.y, pch=19,add=T, col=1)
# box(bty = "l")
# axis(1)
# axis(2)
# grid(equilogs = F)
# abline(h=0.05, lty=3, col="black")
# box()
# legend("topleft", pch = c(1, 19), col = 1, legend=c("All Stations", "Subset"))
# graphics.off()
# ## Figure G15
# #jpeg(filename = paste(wd, "/Results/Figure_G_confounder_tmp.jpg",sep=""), width = 400, height = 400,
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/AppendixG/Figure_G15_confounder_tmp.jpg",sep=""), width = 400, height = 400,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mar = c(5,4,1,1) )
# plot.ecdf(dftv.2$XC95.cdf.x, pch = 1, ylab = "Proportion of species", col = 1, y = c(0,0.5),
#           xlim = range(100, 3200), xlab = "Specific conductivity (?S/cm)",
#           axes = F, log = "x", main= "") #expression(paste("Excluding Temperature > 22", ~degree~C, sep = "")), font.main = 2)                              # quantile xc95
# #            mtext(expression(paste("Excluding Temperature > 22", ~degree~C, sep = "")), font = 2)
# plot.ecdf(dftv.2$XC95.cdf.y, pch=19,add=T, col=1)
# box(bty = "l")
# axis(1)
# axis(2)
# grid(equilogs = F)
# abline(h=0.05, lty=3, col="black")
# box()
# legend("topleft", pch = c(1, 19), col = 1, legend=c("All Stations", "Subset"))
# graphics.off()
# ## Figure G14
# #jpeg(filename = paste(wd, "/Results/Figure_G_confounder_DO.jpg",sep=""), width = 400, height = 400,
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/AppendixG/Figure_G14_confounder_DO.jpg",sep=""), width = 400, height = 400,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mar = c(5,4,1,1) )
# plot.ecdf(dftv.3$XC95.cdf.x,pch=1, ylab="Proportion of species", col=1, y=c(0,0.5),
#           xlim=range(100, 3200), xlab="Specific conductivity (?S/cm)", font.main = 1,
#           axes=F, log="x", main="") # "Excluding DO < 4 mg/L")                              # quantile xc95
# plot.ecdf(dftv.3$XC95.cdf.y, pch=19,add=T, col=1)
# box(bty = "l")
# axis(1)
# axis(2)
# grid(equilogs = F)
# abline(h=0.05, lty=3, col="black")
# box()
# legend("topleft", pch = c(1, 19), col = 1, legend=c("All Stations", "Subset"))
# graphics.off()
# ## Figure G13
# #jpeg(filename = paste(wd, "/Results/Figure_G_confounder_wtshd.jpg",sep=""), width = 400, height = 400,
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/AppendixG/Figure_G13_confounder_wtshd.jpg",sep=""), width = 400, height = 400,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mar = c(5,4,1,1) )
# plot.ecdf(dftv.4$XC95.cdf.x,pch=1, ylab="Proportion of species", col=1, y=c(0,0.5),
#           xlim=range(100, 3200), xlab="Specific conductivity (?S/cm)", font.main = 1,
#           axes=F, log="x", main="") # "Excluding Embeddedness > 75%")                              # quantile xc95
# plot.ecdf(dftv.4$XC95.cdf.y, pch=19,add=T, col=1)
# box(bty = "l")
# axis(1)
# axis(2)
# grid(equilogs = F)
# abline(h=0.05, lty=3, col="black")
# box()
# legend("topleft", pch = c(1, 19), col = 1, legend=c("All Stations", "Subset"))
# dev.off()




###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Save as RDA for use in package#####
# ss and bio.sample
devtools::use_data(dta.do)
devtools::use_data(ss.sites)

