# Prepare data for example for taxon.response.R; taxon.response()
#
# From:
# file:///C:/Users/Erik.Leppo/Documents/GitHub/XC95/data-raw/ProcessData_taxon.response.sort.R
#
# Erik.Leppo@tetratech.com
# 20170417
##~~~~~~~~~~~~~~~~~~~~~

# 0. Prep####
# library(devtools)
# library(reshape)
wd <- getwd()

# 1. Get data and process#####

#  Conductivity document files
#  Lei.Zheng@tetratech.com
#rm(list = ls())
# library(mgcv)
# library(reshape)           # reshape list to dataframe
# library("maps")
# library(Hmisc)
# library(lattice)
# library(DAAG)         # pause
# wd<-getwd()                                 # get working directory
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/indicators.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/curve.shape.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/scatter.plot.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/taxon.response.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/multiplots.r")

############ Metrics and environmental variables
envdata <- read.delim("Env_metrics.txt", header =T); dim(envdata)          # 10455 samples
envdata$lgX1daymax <- log10(envdata$X1daymax)
envdata$lgX3daymax <- log10(envdata$X3daymax)
envdata$lgMAF <- log10(envdata$MeanAnnualFlow)
envdata$lgFallrate <- log10(-envdata$Fallrate)
envdata$lgHigh1fall <- log10(-envdata$High1fall)
envdata$X3daymin[envdata$X3daymin==0] <- 0.01
envdata$lgX3daymin <- log10(envdata$X3daymin)
envdata$Baseflow[envdata$Baseflow==0] <- NA
envdata$lgBaseflow <- log10(envdata$Baseflow)
envdata$lgpro_Spring <- log10(envdata$pro_Spring)
# envdata$pro_Index[envdata$pro_Index>20] <- NA
envdata$lgpro_Index <- log10(envdata$pro_Index)
envdata$lgpro_X3daymin <- log10(envdata$pro_X3daymin)
envdata$lgpro_X3daymax <- log10(envdata$pro_X3daymax)
envdata$lgpro_FlowPredict <- log10(envdata$pro_FlowPredict)

# # regref <- read.delim("RegionalRefSites.txt", header = T); dim(regref)           # 953 reference samples 1416 samples
# #cors <- cor(envdata[varlist],   use = "pairwise.complete.obs", method ="spearman")
# cor1 <- cor(envdata[c("EPTRich120", "pWarm")], envdata[17:125], "pairwise.complete.obs",  "spearman")
# cor1 <- t(data.frame(cor1))
# cor1 <- cor1[order(cor1[,1]),]
# write.table(cor1, "cor1.txt")
# # sort(apply(cors, 2, function(x) mean(abs(x))))

varlist <- c("lgX3daymax", "lgMAF", "lgFallrate", "lgHigh1fall", "RBI",	"lgX3daymin",
             "lgBaseflow",	"lgpro_Spring",	"lgpro_Index", "lgpro_X3daymin",	"lgpro_X3daymax",	"lgpro_FlowPredict")
varnames <- c("3-day Max", "Mean Annual Flow", "Fall Rate", "High 1 Fall", "RBI", "3-Day Min",
              "Baseflow", "Spring Proportion", "Index Proportion", "3-day Min Proportion", "3day Max Proportion",
              "Flow Predicted Proprotion")

my.data <- envdata[c("BenSampID", "Entity", "BioSta_Long", "BioSta_Lat", varlist)]


# for(i in 1:length(varlist)) {
#   hist(envdata[varlist[i]]) #(envdata[varlist[i]])
#   #   print(summary(envdata[varlist[i]]))
#   #   pause()
# }

############ raw taxonomical composition
species.all <- read.delim(file = "QC_Benthic.txt", header = T); dim(species.all)
# ITIS <- read.delim(file = paste(wd,"/data/ITIS.txt",sep=""), header = T); dim(ITIS)
species.all$RA <- species.all$RelAbund/100



###### select samples and species/OTU
species.lim <- merge(species.all, my.data[c("BenSampID", varlist)]); dim(species.lim)

tcount <- tapply(species.lim$BenSampID, species.lim$OTU, length)
tnames <- names(tcount)[tcount>=10]; length(tnames)

species <- subset(species.lim, OTU %in% tnames );  dim(species)
states <- c("Vermont", "New York", "New Hampshire","Rhode Island","Maine","Connecticut", "Massachusetts")


# (entity.n <- table(my.data$Entity))
# (entity <- names(entity.n)[entity.n>100])
# len <- length(entity)
# tile <- c(50, 5)   # either optima or tolerance
# final.perf <- matrix(NA, ncol = 14, nrow = 50)
# colnames(final.perf) <- c("Test_set", "Test_N", "Cal_N", "WAopt","Weighted_WA",
#                           "CDF_ABD","CDF_PA", "CDF_WT_QTILE", "LRM_Full",
#                           "QLRM_Full","GAM_Full", "LRM_Ob", "QLRM_Ob", "GAM_Ob")

# all.i <- 1                             # for row index of r squre file final.perf
# #for(vari in 1:length(varlist))  {      ### for all variables listed
# # 20170417
# vari<-1
#   mydata <- my.data[!is.na(my.data[,varlist[vari]]),]
#   # #for(i in 1:2) {         # optima or tolerance
#   #   for(ind in 1:len) {       # length of test site datasets, only two states qualified as test data
#   #     test.envdata <- subset(mydata, Entity == entity[ind]);dim(test.envdata)
#   #     ind.values <- taxon.response(spdata = species, envdata = test.envdata,  sp.siteid="BenSampID", species="OTU",
#   #                                  sp.abndid="RA", env.siteid="BenSampID", xvar = varlist[vari], cutoff = 20,
#   #                                  region = paste("tolerance_",tile[i],"/", varlist[vari],"_",entity[ind], sep =""),
#   #                                  mtype = 3, dense.N = 201, plot.pdf = F,
#   #                                  add.map = FALSE, log.x = F, rounder=3, taus=c(0,tile[i],100), nbin = 51)
#   #     write.table(ind.values, paste(wd, "/tolerance_",tile[i],"/",varlist[vari],"_", entity[ind], ".indicators.txt",sep=""))
#   #
#   #     if(all.i == 1) {     # save indicate values to all.indicator
#   #       all.indicator <- ind.values[c("taxaname", "WAopt", "CDF_WT_QTILE", "GAM_Full")]
#   #       names(all.indicator)[-1] <- paste(entity[ind], names(all.indicator)[-1],tile[i], sep ="_")
#   #     } else {
#   #       tmp <- ind.values[c("taxaname", "WAopt", "CDF_WT_QTILE", "GAM_Full")]
#   #       names(tmp)[-1] <- paste(entity[ind], names(tmp)[-1], tile[i], sep = "_")
#   #       all.indicator <- merge(all.indicator, tmp, by ="taxaname", all = T)
#   #     }
#   #     if(tile[i]==50)  { # only median do calibration
#   #       pdf(file =  paste(wd, "/tolerance_",tile[i],"/",varlist[vari],"_", entity[ind], "_calibration.pdf",sep=""), width = 9, height = 6.5, pointsize = 12)
#   #       par(mfrow = c(2, 3), pty = "m", mar = c(4, 4, 3, 2),oma = c(1.5, 0.5, 3,0.5) )
#   #       ################## calibration
#   #       cali.set <- subset(species, !(BenSampID %in% test.envdata$BenSampID) & BenSampID %in% mydata$BenSampID); (dim(cali.set))
#   #       cali.set <- merge(cali.set, ind.values, by.x = "OTU", by.y = "taxaname"); dim(cali.set)
#   #       cal.weight <- aggregate(cali.set[c(varlist[vari],colnames(final.perf)[-(1:3)])], list(BenSampID=cali.set$BenSampID), mean); dim(cal.weight)
#   #       cal.weight <- cal.weight[order(cal.weight[,varlist[vari]]),]
#   #       xrange <- range(cal.weight[, varlist[vari]], na.rm =T)
#   #       final.perf[all.i, 1:3] <- c(entity[ind],nrow(test.envdata), nrow(cal.weight))
#   #
#   #       for(index in 3:13) {  # models index in data frame
#   #         yrange <- range(cal.weight[,index])
#   #         y <- (cal.weight[,index] - min(yrange)) * diff(xrange)/diff(yrange) + min(xrange)
#   #         x <- cal.weight[, varlist[vari]]
#   #         plot(y ~ x, xlab = "Observed Value", ylab = "Inferred Value",
#   #              main = names(cal.weight)[index], xlim = c(xrange[1], xrange[2]), pch =21, bg ="gray",
#   #              cex = 0.5, ylim = c(xrange[1], xrange[2]))
#   #         lm1 <- lm(y ~ x)
#   #         legend("bottomright", bty ="n", legend = bquote(italic(r^2) ==.(round(summary(lm1)$r.squared, digit = 2))))
#   #         abline(a=0, b =1, lty = 3, col = 4)
#   #         lines(x, fitted(lm1), lty = 2, col = "red")
#   #         grid()
#   #         final.perf[all.i, 1 + index] <- round(summary(lm1)$r.squared, digit = 3)
#   #       }
#   #       all.i <- all.i + 1
#   #       #     graphics.off()
#   #
#   #     }   # end calibration
#   #   }    # end two states
#     ######## random samples
#     ### random sample indicators
#     test.samp <- sample(mydata$BenSampID, size = round(2/3 * nrow(mydata)), replace=F); length(test.samp)
#     test.set1 <- subset(mydata, BenSampID %in% test.samp); dim(test.set1)
#
#     ind.values <- taxon.response(spdata = species, envdata = test.set1,  sp.siteid="BenSampID", species="OTU",
#                                  sp.abndid="RA", env.siteid="BenSampID", xvar = varlist[vari], cutoff = 20,
#                                  region = paste("tolerance_",tile[i],"/", varlist[vari],"_all", sep =""),
#                                  mtype = 3, dense.N = 201, plot.pdf = F, add.map = FALSE,
#                                  log.x = F, rounder=3, taus=c(0,tile[i],100), nbin = 51)
#     tmp <- ind.values[c("taxaname", "WAopt", "CDF_WT_QTILE", "GAM_Full")]
#     names(tmp)[-1] <- paste("random", names(tmp)[-1], tile[i], sep = "_")
#     all.indicator <- merge(all.indicator, tmp, by ="taxaname", all = T)
#     write.table(ind.values, paste(wd, "/tolerance_",tile[i],"/",varlist[vari],"_random.indicators.txt",sep = ""))
#
#     if(tile[i]==50)  { # only median do calibration
#       pdf(file =  paste(wd, "/tolerance_",tile[i],"/", varlist[vari],"_random.calibration.pdf",sep=""), width = 9, height = 6.5, pointsize = 12)
#       par(mfrow = c(2, 3), pty = "m", mar = c(4, 4, 3, 2),oma=c(1.5, 0.5, 3,0.5) )
#       cal.samp <-  mydata$BenSampID[!(mydata$BenSampID %in% test.samp)]
#       cal.set1 <- subset(species, BenSampID %in% cal.samp); dim(cal.set1)
#       cal.set1 <- merge(cal.set1, ind.values, by.x = "OTU", by.y = "taxaname"); dim(cal.set1)
#       #      cal.wt1 <- aggregate(cal.set1[colnames(final.perf)[4:5]]*cal.set1$RA, list(BenSampID= as.vector(cal.set1$BenSampID)), mean); dim(cal.wt1)
#       cal.wt2 <- aggregate(cal.set1[c(varlist[vari],colnames(final.perf)[-(1:3)])], list(BenSampID= as.vector(cal.set1$BenSampID)), mean); dim(cal.wt2)
#       cal.weight1 <- cal.wt2
#       #      cal.weight1 <- merge(cal.wt1, cal.wt2, by = "BenSampID"); dim(cal.weight1)
#       #      cal.weight1 <- cal.weight1[c("BenSampID", varlist[vari],colnames(final.perf)[-(1:3)])]
#
#       #      cal.weight1 <- aggregate(cal.set1[c(varlist[vari],colnames(final.perf)[-(1:3)])], list(BenSampID= as.vector(cal.set1$BenSampID)), mean); dim(cal.weight1)
#       #      cal.weight1 <- cal.weight1[order(cal.weight1[,varlist[vari]]),]
#       xrange <- range(cal.weight1[, varlist[vari]], na.rm =T)
#       final.perf[all.i, 1:3] <- c(varlist[vari], nrow(test.set1), length(cal.samp))
#       for(index in 3:13) {
#         yrange <- range(cal.weight1[,index], na.rm =T)
#         x <- cal.weight1[, varlist[vari]]
#         y <- (cal.weight1[,index] - min(yrange)) * diff(xrange)/diff(yrange) + min(xrange)
#         plot(y ~ x, data = cal.weight1, xlab = paste("Observed", varlist[vari]), ylab = paste("Inferred", varlist[vari]),
#              main = names(cal.weight1)[index], xlim = c(xrange[1], xrange[2]), pch =21, bg ="gray",
#              cex = 0.3, ylim = c(xrange[1], xrange[2]))
#         lm2 <- lm(y ~x, data = cal.weight1)
#         legend("bottomright", bty ="n", legend = bquote(italic(r^2) ==.(round(summary(lm2)$r.squared, digit = 2))))
#         abline(a=0, b =1, lty = 3, col = 4)
#         lines(x, fitted(lm2), lty = 2, col = "red")
#         grid()
#         final.perf[all.i, index + 1] <- round(summary(lm2)$r.squared, digit = 3)
#       }
#       all.i <- all.i + 1
#       dev.off()
#     }    # end calibration only for median
#
#     whole.values <- taxon.response(spdata = species, envdata = mydata,  sp.siteid="BenSampID", species="OTU",
#                                    sp.abndid="RA", env.siteid="BenSampID", xvar =varlist[vari], cutoff = 20,
#                                    region = paste("tol_all/",varlist[vari], sep =""), xlabs= varnames[vari],
#                                    mtype = 1, dense.N = 201, plot.pdf = T,
#                                    add.map = FALSE,  main = paste("Macroinvertebrates Response to", varnames[vari]),
#                                    log.x = T, rounder=3, taus=c(0,tile[i],100), nbin = 51)
#     tmp <- whole.values[c("taxaname", "CDF_PA", "LRM_Full","N")]
#     names(tmp)[-1] <- paste(names(tmp)[-1], tile[i], sep = "_")
#     all.indicator <- merge(all.indicator, tmp, by ="taxaname", all = T)
#   #}     # end optima and tolerance
# #  graphics.off()
#   #write.table(all.indicator, paste(wd, "/tol_all/",varlist[vari],"all.indicators.txt",sep = ""))
# #}

#write.table(final.perf, paste(wd, "/tol_all/_rsquared.txt",sep = ""))

### final table writing
varlist <- c("lgX3daymax", "lgMAF", "lgFallrate", "lgHigh1fall", "RBI",	"lgX3daymin")
varnames <- c("3-day Max", "Mean Annual Flow", "Fall Rate", "High 1 Fall", "RBI", "3-Day Min")
#all.i <- 1                             # for row index of r squre file final.perf
vari <- 1
i <- 1
#vari <- vari+1
#for(vari in 1:length(varlist))  {      ### for all variables listed
  mydata <- my.data[!is.na(my.data[,varlist[vari]]),]
  #for(i in 1:2) {         # optima or tolerance
    print(varlist[vari])
    whole.values <- taxon.response(spdata = species, envdata = mydata,  sp.siteid="BenSampID", species="OTU",
                                   sp.abndid="RA", env.siteid="BenSampID", xvar =varlist[vari], cutoff = 20,
                                   #region = paste("tol_all/",varlist[vari], sep =""), xlabs= varnames[vari],
                                   region = "tol_all", lim ="GAM", coord = c("BioSta_Long", "BioSta_Lat"), mtype = 3,
                                   dense.N = 201, plot.pdf = T, add.map = F,
                                   # statename and add.lab
                                   statename=NULL, add.lab=F,
                                   main = paste("Macroinvertebrates Response to", varnames[vari]),
                                   # add mar
                                   mar = c(5,4,3, 2),
                                   # modified taus.  #taus=c(0,tile[i],100),
                                   xlabs = varnames[vari], log.x = F, rounder=3, taus=c(0,50,100), nbin = 51,
                                   # added
                                   wd=wd
                                   )

  #   tmp <- whole.values[c("taxaname", "CDF_PA", "LRM_Full","N")]
  #   names(tmp)[-1] <- paste(names(tmp)[-1], tile[i], varlist[vari], sep = "_")
  #   if(all.i == 1 & i ==1) {     # save indicate values to all.indicator
  #     all.indicator <- tmp
  #   } else {
  #     all.indicator <- merge(all.indicator, tmp, by ="taxaname", all = T)
  #   }
  #   all.i <- all.i + 1
  # }
#}
#write.table(all.indicator, paste(wd, "/all_important_indicators.txt",sep = ""))

## 1.5. QC function input ####
    spdata = species
    envdata = mydata
    sp.siteid="BenSampID"
    species="OTU"
    sp.abndid="RA"
    env.siteid="BenSampID"
    xvar =varlist[vari]
    cutoff = 20
    region = "tol_all"
    lim ="GAM"
    coord = c("BioSta_Long", "BioSta_Lat")
    mtype = 3
    dense.N = 201
    plot.pdf = T
    add.map = F
    statename=NULL
    add.lab=F
    main = paste("Macroinvertebrates Response to", varnames[vari])
    mar = c(5,4,3, 2)
    xlabs= varnames[vari]
    log.x = F
    rounder=3
    taus=c(0,50,100)
    nbin = 51
    wd=wd


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Save as RDA for use in package####
# use names in function; spdata and envdata
#species
envdata.all <- my.data
devtools::use_data(species)
devtools::use_data(envdata.all)

