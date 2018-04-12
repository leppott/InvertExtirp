# Prepare data for example for tolerance.R; tolerance()
#
# From SD20150413.R
# V:\Lei.Zheng\20170410_AllFiles\Current_Project\SD_CausalAnalysis\SanDiego2015\Data
#
# Erik.Leppo@tetratech.com
# 20170413
##############################

# 0. Prep####
# library(devtools)
# library(reshape)
wd <- getwd()
#library(XC95)

# 1. Get data and process#####

# library(MatchIt)
# library(mgcv)
 #require(Hmisc)
# library(maps)
# library(rgdal)
# library(maptools)
# library(randomForest)
 #library(reshape)
# #wd <- getwd()
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/scatter.plot2.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/weightcdf.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/curve.shape.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/taxon.response.sort.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/tolerance.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/north.arrow.r")
# source("C:/Documents and Settings/Lei.Zheng/Desktop/R_Function/mapscale.r")

condunit <- expression(paste("Conductivity ( ", mu, "S/cm)", sep = ""))

# Biochem <- read.csv(paste(wd, "/Biochem20150413.csv", sep = ""), header =T); dim(Biochem)
# model <- lm(log10(Cond) ~ log10(TDS), data = Biochem)
#
# #    outliers <- pred[abs(BugEnv[names(pred),"cond"] - pred)>0.4 ]         # outliers
# #   with(Biochem, plot(tds, cond))
# #  with(Biochem, identify(tds, cond, labels = Cond))
#
# Biochem$Cond_cor <- Biochem$Cond
# tds.x <- log10(Biochem$TDS[is.na(Biochem$Cond_cor)&!is.na(Biochem$TDS)])
# Biochem$Cond_cor[is.na(Biochem$Cond_cor)&!is.na(Biochem$TDS)] <- 10^(model$coef[1]+model$coef[2]*tds.x)
#
# Biochem$cond <- log10(Biochem$Cond_cor)
# Landuse <- read.csv(paste(wd, "/07_StationRefStatus_LULC_rm_noland.csv",sep=""),
#                     header =T); dim(Landuse)   # Land use data only, add 4 new target stations with land use information
# Habitat <- read.csv(paste(wd, "/GeoHabitat.csv",sep=""), header =T); dim(Habitat)   #
#
# env <- subset(Biochem, !is.na(cond)); dim(env)  # 542
# env <- merge(env, Landuse[c(1:7,10,13,23,26)],
#              by ="CanonicalStationID", all.x=T); dim(env)
# env <- merge(env, Habitat[c("NewSampID", "W1_HALL", "Slope", "Embed","Sed","PCT_SAFN")], all.x =T); dim(env)
# #write.table(env, paste(wd, "/Results/env.csv", sep=""), sep=",", row.names =F)

 env <- read.csv(paste(wd, "/env.csv", sep=""), header =T); dim(env)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~ propensity score matching ~~~~~~~~~~~~~~~
env.site <- aggregate(env[,c(4:113,119:127)], list(CanonicalStationID= env$CanonicalStationID), mean, na.rm =T)
env.sites <- subset(env, !duplicated(CanonicalStationID))[c("CanonicalStationID","InCaseYN",
                                                            "GIS_InSDR")]; dim(env.sites)  # 367
env.site <- merge(env.site, env.sites); dim(env.site) # 367
env.cond <- subset(env, CanonicalStationID!="SDR-MLS"); dim(env.cond)  # remove cases 497

# selvar <- c("InCaseYN", "CanonicalStationID", "URBAN_2000_WS","Sed", "Embed","RoadDens_WS","Ag_2000_WS","Cond") #"W1_HALL",
#
# #write.table(env.site[selvar], paste(wd, "/Results/env.site.csv", sep =""),sep=",", row.names=F)
#
#
# apply(env.site[env.site$InCaseYN&!is.na(env.site$InCaseYN),],2, function(x) sum(!is.na(x)))
#
# seldata <- na.omit(env.site[selvar])
# row.names(seldata) <- seldata$CanonicalStationID
# imatch <- matchit(InCaseYN ~ URBAN_2000_WS + Sed + Embed, data = seldata, method="nearest" )
#
# summary(imatch)
# m.data <- match.data(imatch)
# env.site[rownames(m.data)[m.data$InCaseYN],"cond"]
# env.site[rownames(m.data)[!m.data$InCaseYN],"cond"]
#
# match.df <- subset(seldata, CanonicalStationID%in%as.vector(imatch$match.matrix))#["Cond"]; dim(match.df)
# org.df <- subset(seldata, CanonicalStationID%in% row.names(imatch$match.matrix))#["Cond"]; dim(match.df)


# ##########################  summary statistics
# sumstats <- function(x) {
#   len <- length(x[!is.na(x)])
#   avg <- as.vector(mean(x, na.rm=T))
#   #   geomean <- exp(mean(log(x), na.rm=T))
#   pct <- c(0, 5, 10, 25, 50, 75, 90, 95, 100)
#   tiles <- as.vector(quantile(x, prob=pct/100, na.rm=T))
#   final <- c(len, avg, tiles)
#   names(final) <- c("Length", "Mean", paste(pct,"th", sep="_"))
#   return(final)
# }
# varlist <- names(env)[c(4:53, 112,119:128)]
# varlist2 <- names(apply(env[varlist], 2, function(x) sum(!is.na(x)))[apply(env[varlist], 2, function(x) sum(!is.na(x)))>200])
#
# listout <- t(apply(env[varlist2], 2, sumstats))
# #write.table(listout,paste(wd, "/Results/Table_summary.csv", sep =""),sep=",")



##~~~~~~~~~~~~~~~~~~~~~~~ tolerance values ~~~~~~~~~~~~~~~~~~~~~

Bugs <- read.csv(paste(wd, "/Invert_Count.csv", sep = ""), header =T); dim(Bugs) # 26277
Bentot <- aggregate(Bugs$Result, list(NewSampID = Bugs$NewSampID), sum, na.rm =T)  # 1107
Bugs <- merge(Bugs, Bentot, by = "NewSampID"); dim(Bugs)
Bugs$RA <- Bugs$Result/Bugs$x

species <- subset(Bugs, (Genus!="")); dim(species)     # 16283
ss <-  reshape::cast(species, NewSampID ~ Genus, sum, value = "RA"); dim(ss)      # 1095 340

# xc.cond <- weightcdf(df1 = env.cond, ss = ss, SampID = "NewSampID", xvar = "cond", nt = 25)
#
# taxalist <- as.vector(xc.cond[order(xc.cond[,"XC95.2"]), "taxaname"])
# df1 <- merge(ss, env.cond[c("NewSampID","cond")], by = "NewSampID"); dim(df1)
# graphics.off()
# full.results <- taxon.response.sort(df1 = df1, xvar="cond", cutoff=25,
#                                     mtype = 3, dense.N = 201, plot.pdf = T, xlabs= condunit,
#                                     add.map = F, maintext = "", region="results/all",
#                                     log.x=TRUE, rounder=0, taus=c(0,95,100), nbin =61,
#                                     sort.vect = taxalist)

# 1.4 Check the function ####
# cond.opt <- tolerance(spdata= ss, envdata= env.cond,  sp.siteid="NewSampID", species="Genus", covar=NULL,
#                       sp.abndid="RA", env.siteid="NewSampID", xvar="cond", cutoff = 25, cutoff2 = 10,
#                       region = "results", lim ="GAM", coord = NULL, mtype = 3, dense.N = 201, cast = FALSE,
#                       plot.pdf = T, add.map = F, statename = NULL,  add.lab = F, add.abund=T,
#                       main = "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient",
#                       mar = c(5,4,3,4), xlabs = condunit, log.x = T,
#                       plus = F, rounder = 3, taus = c(50,95), nbin = 61, wd=getwd()
#                       )
#base = 10,
#write.table(cond.opt, paste(wd, "/Results/",  "opt.cond.csv",sep=""), sep=",", row.names =F)

# #~~~~~~~~~~~~~~
# # 1.5. QC function ####
# spdata <- ss
# envdata<- env.cond
# sp.siteid<-"NewSampID"
# species<-"Genus"
# covar<-NULL
# sp.abndid<-"RA"
# env.siteid<-"NewSampID"
# xvar<-"cond"
# cutoff <- 25
# cutoff2 <- 10
# region <- "results"
# lim <-"GAM"
# coord <- NULL
# mtype <- 3
# dense.N <- 201
# cast <- FALSE
# plot.pdf <- T
# add.map <- F
# statename <- NULL
# add.lab <- F
# add.abund<-T
# main <- "Capture Probability of Macroinvertebrate Taxon Along Conductivity Gradient"
# mar <- c(5,4,3,4)
# xlabs <- condunit
# log.x <- T
# #base <- 10
# plus <- F
# rounder <- 3
# taus <- c(50,95)
# nbin <- 61
# wd<-getwd()
# i<-22

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Save as RDA for use in package####
# ss and env.data
tol.ss <- ss
tol.env.cond <- env.cond
devtools::use_data(tol.ss)
devtools::use_data(tol.env.cond)

