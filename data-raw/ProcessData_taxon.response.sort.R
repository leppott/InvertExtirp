# Prepare data for example for taxon.response.sort.R; taxon.response.sort()
#
# From _Appendix_F_all.R
#
# Erik.Leppo@tetratech.com
# 20170406
##~~~~~~~~~~~~~~~~~~~~~

# 0. Prep####
# library(devtools)
# library(reshape)
wd <- getwd()

# 1. Get data and process#####


# #~~~~~~~~~~~~~~~~~~~~~~
# #
# rm(list = ls())
# # library(mgcv)
# library(lattice)
# library(reshape)           # reshape list to dataframe
# library("maps")
# library(Hmisc)
# library(mgcv)
# library(arm)
# wd<-getwd()                                 # get working directory  # QC, mod directory structure
# source(paste(wd,"/R_Function/threeaxes.r",sep=""))
# source(paste(wd,"/R_Function/scatter.plot2.r",sep=""))
# source(paste(wd,"/R_Function/xyplots.r",sep=""))
# source(paste(wd,"/R_Function/weightcdf2.r",sep=""))
# source(paste(wd,"/R_Function/multiplots.r",sep=""))
# source(paste(wd,"/R_Function/taxon.response.sort3.r",sep=""))
# source(paste(wd,"/R_Function/curve.shape.r",sep=""))
# t1 <- proc.time()
#
# ##~~~~~~~~~~~~~~~ do with the new dataset
# #species.all <- read.delim(file = paste(wd,"/data/BenthicCompositionGenera.txt",sep=""), header = T); dim(species.all)
# species.all <- read.delim(file = paste(wd,"/data-raw/BenthicCompositionGenera.txt",sep=""), header = T); dim(species.all)
# mydate <- strptime(species.all$SAMPLE_DATE, format="%m/%d/%Y")
# species.all$Year <- mydate$year + 1900
# species.all$Month <- mydate$mon + 1
# species.all$Julian <- as.numeric(format(mydate, "%j"))
# species.all$NumMonth <- species.all$Month
# species.all$myDate <- as.Date(species.all$SAMPLE_DATE, format="%m/%d/%Y")
# species.all$Sta_Year <- as.factor(paste(species.all$STATION_CODE, species.all$Year, sep="_"))
# newsum <- aggregate(species.all[,"Benthic.Count"], list(Sample.ID = species.all$Sample.ID), sum, na.rm =T)
# species.all <- merge(species.all, newsum)
# species.all$RA <- species.all$Benthic.Count/species.all$x
# ITIS.genera <- read.delim(paste(wd, "/data/ITIS_genera.txt", sep=""), header = T); dim(ITIS.genera)
#
# xc95 <- read.table(paste(wd, "/ecoregion/xc95_nowt.txt",sep=""), header = T)
# xc95 <- xc95[order(xc95$Eco69_xc95),]
# sens.taxa69 <- xc95$taxaname[xc95$Eco69_xc95<310]
# xc95 <- xc95[order(xc95$Eco70_xc95),]
# sens.taxa70 <- xc95$taxaname[xc95$Eco70_xc95<340]
# new.sens <- subset(species.all, GENUS %in% sens.taxa69 & Level3 == 69); dim(new.sens)
# new.sens$Sta_Year <- factor(new.sens$Sta_Year)
# new.sta69 <- aggregate(new.sens$GENUS, list(Sta_Year= new.sens$Sta_Year, Sample.ID= new.sens$Sample.ID), length); dim(new.sta69)  # 2379
# names(new.sta69)[3] <- "SensTaxa"
# new.sens <- subset(species.all, GENUS %in% sens.taxa70 & Level3 == 70); dim(new.sens)
# new.sens$Sta_Year <- factor(new.sens$Sta_Year)
# new.sta70 <- aggregate(new.sens$GENUS, list(Sta_Year= new.sens$Sta_Year, Sample.ID= new.sens$Sample.ID), length); dim(new.sta70)   # 2127
# names(new.sta70)[3] <- "SensTaxa"
# new.sta <- rbind(new.sta69, new.sta70) ; dim(new.sta)
#
# fullset <- read.delim(file = paste(wd,"/data/Data_combined.txt", sep =""), header = T); dim(fullset)
# fullset$RBP_7Sc <- with(fullset, Fish.Cover + Embeddedness + Channel.Alteration + Sediment.Deposition +
#                           Total.Bank.Stability + Total.Bank.Vegetation + Total.Undisturbed.Vegetation)
# RBP_avg <- aggregate(fullset["RBP_7Sc"], list(fullset$STATION_ID), mean, na.rm =T); dim(RBP_avg)
# names(RBP_avg) <- c("STATION_ID", "RBP_7avgSc")
# fullset <- merge(fullset, RBP_avg, by = "STATION_ID"); dim(fullset)
#
# fullset$Conductivity[fullset$Conductivity ==0] <- NA
# fullset$Alkalinity[fullset$Alkalinity < 1]  <- NA
# fullset$Na_Tot[fullset$Na_Tot < 0.1] <- NA
# fullset$Fecal[fullset$Fecal==0] <- 0.5               # set 0 to 0.5
# fullset$Ca_Tot[fullset$Ca_Tot < 0.5] <- NA
# fullset$Mg_Tot[fullset$Mg_Tot < 0.5] <- NA
# fullset$Fe_Dis[fullset$Fe_Dis < 0.01] <- NA
# fullset$Fe_Tot[fullset$Fe_Tot < 0.01] <- NA
# fullset$Al_Dis[fullset$Al_Dis < 0.01] <- NA
# fullset$Al_Tot[fullset$Al_Tot < 0.01] <- NA
# fullset$TP[fullset$TP==0] <- NA
# fullset$Hardness[fullset$Hardness < 1] <- NA
# fullset$Chloride[fullset$Chloride_Total_Q =="<"] <- fullset$Chloride[fullset$Chloride_Total_Q=="<"]/2
# fullset$Chloride_Total_Q[fullset$Chloride_Total_Q !="<"&fullset$Chloride_Total_Q !=""] <- NA
# fullset$Bicarbonate <- fullset$Alkalinity * 61 *2/100
# fullset$SO4HCO3 <- (fullset$Sulfate + fullset$Bicarbonate)
# fullset$lgSO4HCO3 <- log10(fullset$SO4HCO3)
# fullset$SO4CL <- (fullset$Sulfate + fullset$Chloride)
# fullset$HCO3_eq <- fullset$Alkalinity/100*2            # equivlent
# fullset$SO4_eq <- fullset$Sulfate/96*2                 # equivlent
# fullset$HCO3SO4_eq <- (fullset$HCO3_eq + fullset$SO4_eq)      # equivlent
# fullset$Cl_eq <- fullset$Chloride/35.45               # equivlent
# fullset$SO4_mol <- fullset$Sulfate/96
# fullset$HCO3SO4_mol <- (fullset$HCO3_eq + fullset$SO4_mol)
# fullset$lgSO4CL <- log10(fullset$SO4CL)
# fullset$Anion <-  fullset$SO4HCO3 + fullset$Chloride
# fullset$lgAnion <- log10(fullset$Anion)
# fullset$Cation <- fullset$Ca_Tot + fullset$Mg_Tot + fullset$Na_Tot + fullset$K_Tot
# fullset$CaMg <-  fullset$Ca_Tot + fullset$Mg_Tot
# fullset$cond <- log10(fullset$Conductivity)
# fullset$lgFecal <- log10(fullset$Fecal)
# fullset$Ratio <-  fullset$SO4HCO3/fullset$Chloride
# mydate <- strptime(fullset$SAMPLE_DATE, format="%m/%d/%Y")
# fullset$Year <- mydate$year + 1900
# fullset$Month <- mydate$mon + 1
# fullset$Julian <- as.numeric(format(mydate, "%j"))
# fullset$NumMonth <- fullset$Month
# fullset$Year2 <- fullset$Year
# fullset$Year2[fullset$Month >7 & !is.na(fullset$Month)] <-
#   fullset$Year2[fullset$Month >7 & !is.na(fullset$Month)] +1 # move last year summer samples to this year summer
# fullset$Sta_Year <- as.factor(paste(fullset$STATION_CODE, fullset$Year, sep="_"))
# fullset$Sta_Year2 <- as.factor(paste(fullset$STATION_CODE, fullset$Year2, sep="_"))
# fullset$Julian2 <- fullset$Julian
# fullset$Julian2[fullset$Julian2>182 & !is.na(fullset$Julian2)] <- fullset$Julian2[fullset$Julian2>182&!is.na(fullset$Julian2)]-365
# fullset$myDate <- as.Date(fullset$SAMPLE_DATE, format="%m/%d/%Y")
# fullset <- merge(fullset, new.sta[c("Sta_Year","SensTaxa")], all.x =T, by.x = "Sta_Year2", by.y = "Sta_Year"); dim(fullset)
# fullset$SensTaxa[is.na(fullset$SensTaxa)] <- 0
#
# ################## biological samples with conductivity
# envbio <- subset(fullset, !is.na(PCT_OF_THRESHOLD_GLIMPSS_CF) & !is.na(cond) & !is.na(Month) & !duplicated(Sample.ID)&
#                    (Level3 == 69 | Level3 == 70)); dim(envbio)
#
# biosta <- subset(envbio, !duplicated(STATION_ID)); dim(biosta)
# bio.sta <- aggregate(envbio[c("lgSO4HCO3", "WVSCI", "PCT_OF_THRESHOLD_GLIMPSS_CF")], #  "SensTaxa")],
#                      list(Sta_Year = envbio$Sta_Year, Sample.ID = envbio$Sample.ID), mean, na.rm =T); dim(bio.sta)
# bio.date <- aggregate(envbio["Month"], list(Sta_Year= envbio$Sta_Year), min, na.rm =T)  # only one biosample in a year
# bio.sta <- merge(bio.sta, bio.date, all.x =T)
#
# ################# chemistry samples in the two ecoregions
# chemistry <- subset(fullset, (Level3 == 69 | Level3 == 70)&!is.na(Conductivity)); dim(chemistry)
# chem.sta <- subset(chemistry, !duplicated(STATION_ID)); dim(chem.sta)
#
# ##################number of total stations and samples
# envbio69 <- subset(chemistry, !is.na(PCT_OF_THRESHOLD_GLIMPSS_CF) &Level3 == 69 & !duplicated(Sample.ID)&!is.na(Month)); dim(envbio69) #1911
# envbio69 <- subset(envbio69, pH>6);dim(envbio69)       # 1674
# envbio69 <- subset(envbio69, Ratio>1 ); dim(envbio69) #1661
# envbio69 <- subset(envbio69, !is.na(cond)); dim(envbio69) #1661
# length(table(envbio69$STATION_CODE)[table(envbio69$STATION_CODE)>1])    # number of stations with replicates    159
# length(table(envbio69$STATION_CODE)[table(envbio69$STATION_CODE)>0])     # total number of stations  1420
# write.table(envbio69, paste(wd, "/Ion_HC05/env.bio69.csv",sep=""), sep=",", row.names =F)
# envbio70 <- subset(chemistry, !is.na(PCT_OF_THRESHOLD_GLIMPSS_CF) &Level3 == 70 & !duplicated(Sample.ID)&!is.na(Month)); dim(envbio70) #2126
# envbio70 <- subset(envbio70, Ratio>1); dim(envbio70) #2123
# envbio70 <- subset(envbio70, pH>6);dim(envbio70)       # 2075
# envbio70 <- subset(envbio70, !is.na(cond)); dim(envbio70) #2075
# write.table(envbio70, paste(wd, "/Ion_HC05/env.bio70.csv",sep=""), sep=",", row.names =F)
# sum(!is.na(envbio70$Ratio)&!is.na(envbio70$Ca_Tot))/nrow(envbio70)
# sum((envbio70$Ratio>1)&!is.na(envbio70$Ca_Tot),na.rm =T)
#
#
# length(table(envbio70$STATION_CODE)[table(envbio70$STATION_CODE)>1])     # number of staions with replicates 317
# sum(table(envbio70$Sta_Year2)>0)    # total number of station year 2058, ouccred more than once 17
# length(table(envbio70$STATION_CODE)[table(envbio70$STATION_CODE)>0])     # total number of stations 1695
#
# length(unique(fullset$STATION_ID[fullset$Level3==69])); length(unique(fullset$STATION_ID[fullset$Level3==70]))   # 2299 2011 number of stations
# sum(table(fullset$STATION_ID[fullset$Level3==69])>0)     # 2299
# sum(table(fullset$STATION_ID[fullset$Level3==70])>0)    # 2011
# sum(fullset$Level3==69); sum(fullset$Level3==70)             # total number of samples   9824 12915
# sum(biosta$Level3==69) ;sum(biosta$Level3==70)                   ###  1620 and 1727     number of genus level biostation with Month
# sum(envbio$Level3==69); sum(envbio$Level3==70)          ###     1911 and 2126       number of genus level biossample with Month
# sum(table(envbio$STATION_ID[envbio$Level3==69])>1)/1620
# sum(table(envbio$STATION_ID[envbio$Level3==69])>1)/1727
#
# sum(envbio69$Ratio>1 & envbio69$Ca_Tot, na.rm =T)         # 825    + 13
# sum(envbio70$Ratio>1 & envbio70$Ca_Tot, na.rm =T)          # 867   + 3
#
# sta.mean <- aggregate(chemistry[c("lgFecal","pH", "RBP_7avgSc")], list(STATION_CODE = chemistry$STATION_CODE), mean, na.rm =T)
# sta.mean$Bad <-  sta.mean$lgFecal>3 | sta.mean$pH < 6 | sta.mean$RBP_7avgSc <91             # determine bad station impaired by other stressors
# chemistry <- merge(chemistry, sta.mean[c("STATION_CODE", "Bad")], by = "STATION_CODE");dim(chemistry)
# ####################### Define biological samples then weighting the conductivity values
# bio.sample <- subset(envbio, pH>6 & (Ratio > 1|is.na(Ratio))); dim(bio.sample)      # 3736
# sum(bio.sample$Level3==69); sum(bio.sample$Level3==70)          ###     1661 and 2075
#
# #   write.table(bio.sample, paste(wd,"/IntermediaFile/bio.sample.csv", sep =""), sep =",", row.names =T)
# bio.sample <- bio.sample[order(bio.sample$STATION_ID, bio.sample$Year),]
# # bio.sample <- subset(bio.sample, !duplicated(STATION_ID, fromLast=TRUE))    # 3115
# my70 <- subset(bio.sample, Level3==70 &Ratio>1); dim(my70)              # 1042
# my69 <- subset(bio.sample, Level3==69 & Ratio>1); dim(my69)              # 925
# sum(table(my69$STATION_ID)>1)
#
# ##~~~~~~~~~~~~
# ##~~~~~~~~~~~~
# # QC, 20160622, Table F-1, code from Table 401
#
# varlist <-c("Conductivity", "Hardness","Alkalinity","Sulfate", "Chloride","SO4HCO3", "Ca_Tot", "Mg_Tot", "Na_Tot",
#             "K_Tot","TSS","Fe_Tot","Fe_Dis","Al_Tot","Al_Dis", "Mn_Tot", "Se_Tot","DO","TP", "NOx","Fecal", "pH","SQKM", "Temperature",
#             "RBP_Sc","RBP_7Sc","Embeddedness","X..Fines..sand.silt." )
#
# ############### Table   summary statistics
# summ <- function(x) {
#   x <- x[!is.na(x)]
#   len <- round(length(x))
#   tiles <- round(quantile(x, prob = c(0, 0.25, 0.5, 0.75, 1)),3)
#   rmean <- round(mean(x),3)
#   gmean <- round(exp(mean(log(x[x>0]))),3)
#   y <- c(len, tiles, rmean, gmean)
#   names(y) <- c("N", "Min", "25th", "50th","75th", "Max", "Mean", "Geomean")
#   return(y)
# }#END.FUNCTION.summ
#
# my.data <- my69
# listout <- t(apply(my.data[varlist], 2, summ))
#
# #
# #   eco3 <- ifelse(switch0 == 1, 69, 70)
# #   if(switch0 == 1) myTableNum <- "401"
# #   if(switch0 == 2) myTableNum <- "501"
#
# myTableNum <- "F01"
#
# write.csv(listout,paste(wd, "/_Tables_Figures_Final/","Table_",myTableNum,"_summary.csv", sep =""))
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ##~~~~~~~~~~~~~~~~~~Ecoregion begin ~~~~~~~~~~~~~~~~~~~~~~~~
# ref.site <- subset(bio.sample, REFERENCE == "Level I"); dim(ref.site)      # 133
#
# ref.taxa <- unique((subset(species.all, GENUS !="" & Sample.ID %in% ref.site$Sample.ID))$GENUS) ; length(ref.taxa)
# ref.taxa69 <- unique((subset(species.all, GENUS !="" & Sample.ID %in% ref.site$Sample.ID & Level3 == 69 ))$GENUS) ; length(ref.taxa69) # 193
# ref.taxa70 <- unique((subset(species.all, GENUS !="" & Sample.ID %in% ref.site$Sample.ID & Level3 == 70 ))$GENUS) ; length(ref.taxa70) # 179
# ref.taxa.all <- ITIS.genera[ITIS.genera$Genus %in% ref.taxa,2:4]
# ref.fam <- unique((subset(species.all, FAMILY !="" & Sample.ID %in% ref.site$Sample.ID))$FAMILY) ; length(ref.fam)
# all.taxa <- unique((subset(species.all, GENUS !="" & Sample.ID %in% bio.sample$Sample.ID))$GENUS); length(all.taxa) #476
# ref.genera <- as.vector(ref.taxa[regexpr("/", ref.taxa)==-1]); length(ref.genera)       # 224
#
# species <- subset(species.all, GENUS %in% ref.genera & Sample.ID %in% envbio$Sample.ID &
#                     !is.na(RA)); dim(species)
# species <- merge(chem.sta[c("STATION_ID", "STATION_CODE")], species); dim(species)
# ss <- cast(species, Sample.ID ~ GENUS, sum, value = "RA");dim(ss)
#
# ##### ecoregion 70 and 69 have similar taxa
#
# taxa69 <- sort(unique(species[species$Level3==69,"GENUS"]))    # 227
#
# unitlab <- expression(paste("SO"[4]^{2-phantom()}," + HCO"[3]^{-phantom()}, " (mg/L)"))
#
# #### plot two ecoregion seperately
# my69.ref <- subset(ref.site, Level3==69); dim(my69.ref);sum(table(my69.ref$STATION_ID)>0)  # 87 / 64               # 87 samples,   64 stations
# my69.spr <- subset(my69, Month > 2 & Month < 7); dim(my69.spr)           # 627
# my69.sum <- subset(my69, Month > 6 & Month < 11); dim(my69.sum)          # 1016
# length(unique(my69.ref$STATION_ID))
#
# ##~~~~~~~~~~~~~~~~~~~~~  STEP 1 Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~
# ### Step 1A. use samples to calcualte HC05
# ### basic models for season, and ecoregion
# allyear <- weightcdf(df1 = bio.sample[!is.na(bio.sample$lgSO4HCO3),], ss = ss, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = 25)
# apply(allyear[-1], 2, quantile, prob =0.05, type =6, na.rm = T)
#
# xc70 <- weightcdf(df1 = my70, ss = ss, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = 25)
# names(xc70) <- c("taxaname","Eco70_xc95","Eco70_N")
#
# xc69 <- weightcdf(df1 = my69, ss = ss, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = 25)
# names(xc69) <- c("taxaname","Eco69_xc95",  "Eco69_N")
#
# xc69.spr <- weightcdf(df1 = my69.spr, ss = ss, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = 25)
# names(xc69.spr) <- c("taxaname","Eco69.spr_xc95",  "Eco69.spr_N")
# xc69.sum <- weightcdf(df1 = my69.sum, ss = ss, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = 25)
# names(xc69.sum) <- c("taxaname","Eco69.sum_xc95",  "Eco69.sum_N")
#
# #### merge them all togetther
# xc6970 <- merge(xc69, xc70, by = "taxaname", all = T)
#
# xc95 <- merge(allyear, xc6970, all = T, by = "taxaname")
# xc.season <- merge(xc69.spr, xc69.sum, all =T); dim(xc.season)
#
# xc95 <- merge(xc95, xc.season, all =T)
# xc95 <- merge(xc95, ITIS.genera[c("Genus","Order","Family")], all.x =T, by.x ="taxaname", by.y = "Genus"); dim(xc95)
# hc05.all <- c(round(apply(xc95[, c(2, 4,6, 8, 10)], 2, quantile, prob = 0.05, type = 6, na.rm = T)),
#               apply(xc95[, c(3, 5,7, 9, 11)],2, function(x) sum(!is.na(x))) )
#
# write.table(xc95, paste(wd, "/Ion_HC05/xc95_nowt.txt", sep=""))
# write.table(hc05.all, paste(wd, "/Ion_Hc05/hc05_nowt.txt", sep=""))
#
# ##~~~~~~~~~~~~~~~~~~~~~~~ now with single ecoregions ~~~~~~~~~~~~~~~~~~~~
# switch0 <- 1
#
# ecolab <- ifelse (switch0 ==1, "/Ion_HC05/eco69/", "/Ion_HC05/eco70/")
# eco3 <- ifelse(switch0 == 1, 69, 70)
# if(switch0 == 1) my.data <- my69
# if(switch0 == 2) my.data <- my70
# unitlab <- expression(paste("SO"[4]^{2-phantom()}," + HCO"[3]^{-phantom()}," (mg/L)"))
# tiff(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F01_Ion_vs_Cond.tiff",sep=""), width = 480, height = 480,
#      pointsize = 13,  restoreConsole = TRUE)
# #par(mfrow = c(3,2), mar = c(5,4.5,0.1,1))
# par(mfrow = c(2,2), mar = c(5,4.5,0.1,1))
# scatter.plot(my.data, xvar ="Conductivity", yvar ="Sulfate", add.fit = "linear", xlab ="Conductivity (?S/cm)",
#              ylab = expression(paste("SO"[4]^{2-phantom()},"(mg/L)")), log.flag ="xy", cex = 0.5, add.r2=T, pch =21)
# scatter.plot(my.data, xvar ="Conductivity", yvar ="Bicarbonate", add.fit = "linear",xlab ="Conductivity (?S/cm)",
#              ylab = expression(paste("HCO"[3]^{-phantom()}, " (mg/L)")), log.flag ="xy", add.r2=T,cex = 0.5,  pch =21)
# scatter.plot(my.data, xvar ="Conductivity", yvar ="SO4HCO3", xlab ="Conductivity (?S/cm)", add.reg =T, cex = 0.5,
#              add.fit = "linear", ylab = expression(paste("SO"[4]^{2-phantom()}," + HCO"[3]^{-phantom()}," (mg/L)")),
#              log.flag ="xy", pch = 21, add.r2=T)
# scatter.plot(my.data, xvar ="Conductivity", yvar ="Anion", xlab ="Conductivity (?S/cm)", add.reg =T,cex = 0.5,
#              add.fit = "linear", ylab = unitlab, log.flag ="xy", pch = 21, add.r2=T)
# dev.off()
#
#
#
#
# ###### Figure 3 box plot of monthly variation of random sites
#
# random <- subset(chemistry, Level3== eco3 & !is.na(SO4HCO3) & Med_pH > 6 &
#                    (SURVEY_TYPE == "WAP-Random")); dim(random)         # 585  681 samples respetive 544 and 617 sites,
# month.lab <-  month.abb[c(sort(unique(random$Month)))]
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F03_random_sites.jpg",sep=""), width = 480, height = 480,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mar = c(5,5,1,1))
# boxplot(SO4HCO3 ~ Month, data = random, ylab = unitlab,las =1,
#         xlab = "Month", log = "y", col = "lightgray", axes = F)
# axis(1, at = 1:length(month.lab), labels = month.lab)
# axis(2, las = 1)
# #   abline(h=200, lty =2)
# box()
# dev.off()##JPEG.END
#
# ####################### Figure 4  conductivity monthly variation
#
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F04_bxplot_ref.jpg",sep=""), width = 480, height = 480,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# ref.sub <- subset( chemistry, REFERENCE=="Level I" & Level3 == eco3 & !is.na(cond)) ; dim(ref.sub)     # 69: 112 samples, 82 ref sites  70: 30 stations, 48 samples,
# length(unique(ref.sub$STATION_ID))         # 82 ref sites, 30 stations
# month.lab <-  month.abb[c(sort(unique(ref.sub$Month)))]
# par(mar = c(5,5,1,1))
# boxplot(SO4HCO3 ~ Month, data = ref.sub, ylab = unitlab,las=1,
#         xlab = "Month", log = "y", col = "lightgray", axes = F)
# axis(1, at = 1:length(month.lab), labels = month.lab)
# axis(2, las = 1)
# #  abline(h=200, lty =2)
# box()
# dev.off()##JPEG.END
#
# #################  Figure 5 conductivity monthly variation in all sites
# month.lab <-  month.abb[c(sort(unique(my.data$Month)))]
#
# par(mar = c(5,5,1,1))
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F05_bxplot_allsites.jpg",sep=""), width = 480, height = 480,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mar = c(5,5,1,1))
# boxplot(SO4HCO3 ~ Month, data = my.data, ylab = unitlab,las =1,
#         xlab = "Month", log = "y", col = "lightgray", axes =F)
# axis(1, at = 1:length(month.lab), labels = month.lab, cex.axis =1)
# axis(2, las = 1)
# #   abline(h=200, lty =2)
# box()
# dev.off()##JPEG.END
#
# ##~~~~~~~~~~~~~~~ Figure 6  Histogram of conductivity   ~~~~~~~~~~~~~~~~
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F06_hist.jpg",sep=""), width = 480, height = 480,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
#
# hist(my.data$lgSO4HCO3, xlab=unitlab, breaks = 61, xaxt = "n", main = "")
#
# max.pow <- ceiling(max(my.data$lgSO4HCO3)); min.pow <- floor(min(my.data$lgSO4HCO3))
# axis.at <- 10 ^ c(min.pow:max.pow)
# axis(1, at = min.pow:max.pow, labels = axis.at)
# axis(1, at = log10(1:10 * rep(axis.at[-1] / 10, each = 10)), tcl = -0.5, labels = FALSE)
#
# axis(2)
# box(bty="l")
# dev.off()##JPEG.END
#
# ###########    Figure_F 7 species sensitivity distribution
# xc95 <- read.table(paste(wd, "/Ion_HC05/", "xc95_nowt.txt",sep=""), header = T)
#
#
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F07_ssd.jpg", sep = ""), width = 480, height = 480,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# plot.ecdf(xc95[,paste("Eco",eco3, "_xc95",sep="")], pch = 19, ylab = "Proportion of genera extirpated", col = 1, ylim = c(0,1),
#           xlim=range(xc95[,paste("Eco",eco3, "_xc95",sep="")], na.rm =T), xlab = unitlab, bg = 8,
#           axes=F, log="x", main="", cex.lab = 1.2)                              # quantile xc95
# hc05 <- round(quantile(xc95[,paste("Eco",eco3, "_xc95",sep="")], prob = 0.05, type = 6, na.rm = T))
# box(bty = "l")
# axis(1)
# axis(2)
# arrows(hc05,0.3, hc05,0.05, length = 0.1,col = "black")
# text(hc05,0.3, labels = paste(hc05,"mg/L"),pos=3, col="black")
# abline(h = 0.05, lty = 3, col = "black")
# box()
# dev.off() ##JPEG.END
#
#
# ############## Figure_F 8 SSD for all year  ######### 1/3 of figure 8
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F08_ssd_part.jpg",sep = ""), width = 500, height = 480,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# plot.ecdf(xc95[,paste("Eco",eco3, "_xc95",sep="")], pch=19, ylab="Proportion of genera extirpated", col=1, ylim = c(0, 0.3),
#           xlim = range(sort(xc95[,paste("Eco",eco3, "_xc95",sep="")])[1:30]), xlab = unitlab,
#           axes = F, log = "x", main = "", cex.lab = 1.2)                              # quantile xc95
# hc05 <- round(quantile(xc95[,paste("Eco",eco3, "_xc95",sep="")], prob = 0.05, type = 6, na.rm = T))
# axis(1)
# axis(2)
# abline(h = 0.05, lty = 3, col = "black")
# box(bty = "l")
# arrows(hc05,0.2, hc05, 0.05, length = 0.1,col = "black")
# text(hc05,0.2, labels = paste(hc05,"mg/L"),pos=3, col="black")
# abline(h = 0.05, lty = 3, col = "black")
# box(bty="l")
# #   mtext(hc05, side = 1, line = 1, at = hc05, cex = 1.5, col = "red")
# dev.off() ##JPEG.END
#
# ###~~~~~~~~~~~~~~ determine sample size #~~~~~~~~~~~~~~~~
# # foo <- function() {
# nn <- seq(5, 60, 5)
# res <- matrix(NA, nrow=length(nn), ncol = 3)
# res[,1] <- nn
# for(nind in 1:length(nn)) {
#   xctmp <- weightcdf(df1 = my.data, ss = ss, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = nn[nind])
#   res[nind, 2] <- quantile(xctmp[,2], prob = 0.05, type = 6, na.rm =T)
#   res[nind, 3] <- nrow(xctmp)
# }
# res <- data.frame(res)
# tiff(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F11_Samplereq_vs_HC.tiff",sep=""), width = 450, height = 400,
#      pointsize = 13,  restoreConsole = TRUE)
#
# threeaxes.plot(my.data = res, xvar=1, yvar1=2, yvar2=3, y2log =F, eq.scale =F,
#                xlab = "Number of Minimum Required Samples", mar = c(5,5,1,5),
#                y1lab = unitlab, y2lab = "Number of Taxa in XCD")
# dev.off()##TIFF.END
#
# #}
# #################  step 1b...HC05 model  using all taxa calcualtion (not only reference taxa)
#
# if(switch0 ==1) all.genera <- taxa69
#
# nonambiguous <- all.genera[regexpr("/", all.genera)==-1]
# species11 <- subset(species, (GENUS %in% nonambiguous)); dim(species11)
# ss11 <- cast(species11, Sample.ID ~ GENUS, sum, value = "RA");dim(ss11)
# all11 <- weightcdf(df1 = my.data, ss = ss11, SampID = "Sample.ID", xvar = "lgSO4HCO3", nt = 25)
# apply(all11[-1], 2, quantile, prob =0.05, type =6, na.rm = T)       # 296 171taxa
#
#
# ######## plot sample size vs. HC05
# hc05.conf <- read.table(paste(wd, ecolab, "hc05_samplesizeissue.txt", sep=""), header =T )
# nt.conf <- read.table(paste(wd, ecolab, "ntaxa_samplesizeissue.txt", sep=""), header =T )
# hc <- cbind(hc05.conf, nt.conf)
# names(hc) <- c("hc025", "hc05", "hc50", "hc95", "hc975","n025", "n05", "n50", "n95", "n975")
# hc.scale <- (hc05.conf- min(hc05.conf))/(max(hc05.conf)-min(hc05.conf))
# nt.scale <- (nt.conf- min(nt.conf))/(max(nt.conf)-min(nt.conf))
# N <- c(100, 200, 300, 500, 800, 925)
# j <- ifelse(max(hc05.conf) - min(hc05.conf)>100, 50, 10)         # increaser steps
# at1min <- floor(min(hc05.conf)/j)*j; at1max <- floor(max(hc05.conf)/j)*j     # min and max scale
# len1 =9                                                                # initial intervals (maximum)
# while((at1max -at1min)/(j*len1) !=round((at1max -at1min)/(j*len1))) len1= len1 -1    # perfect intervals
# at1.lab <- seq(at1min, at1max, length = len1+1)                   # ticks label for axis 2
# at1 <- (at1.lab- min(hc05.conf))/(max(hc05.conf)-min(hc05.conf))  # ticks become original 0 to 1 scale
#
# at2.lab <- round(at1*(max(nt.conf)-min(nt.conf)) + min(nt.conf))        # ticks 4 become axis 4 label
#
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", "Figure_F12_SampleSizeIssue.jpg",sep=""), width = 500, height = 450,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mar = c(5,5, 2, 5))
# plot(hc.scale[,3] ~ rownames(hc.scale), ylab = unitlab,
#      xlab = "Sample Size", col ="red", bg = "gray", pch = 21, ylim = range(hc.scale),
#      las = 1, axes = F)
# box()
# grid(nx=NULL, ny=NA)
# mtext("Number of genera in XCD", side = 4, line = 2.2)
# abline(h=at1, col = "gray", lty =3)
# segments(N, hc.scale[,1], N, hc.scale[,5], col = "red")
# axis(1)
# axis(2, at = at1, labels = at1.lab, las=1)         # round(c(min(hc05.conf)+ axTicks(2) * (max(hc05.conf)-min(hc05.conf)))), las =1)
# axis(4, at  = at1, labels = at2.lab, las=1) # round(c(min(nt.conf)+ axTicks(4) * (max(nt.conf)-min(nt.conf)))), las =1)
# points(N + 10, nt.scale[,3], pch =24, col = "blue", bg = "lightgray")
# segments(N + 10, nt.scale[,1], N+10, nt.scale[,5], col = "blue")
#
# dev.off() ##JPEG.END
#
# ############ updated MG dataset from 201507
# ########################## 1.1 summary statistics
# sumstats <- function(x) {
#   len <- length(x[!is.na(x)])
#   pct <- c(0,10, 25, 50, 75, 100)
#   tiles <- as.vector(quantile(x, prob=pct/100,type = 5, na.rm=T))# type 5 is only for EPA new data
#   final <- c(tiles,len)
#   names(final) <- c( paste(pct,"th", sep="_"),"Length")
#   return(final)
# }
#
# newepa.df <- read.csv(paste(wd, "/Cascades/EPA/mg_l_MG20150729.csv", sep=""), header =T); dim(newepa.df)
# newepa.df$SO4[(newepa.df$SO4==0)] <- NA
# newepa.df$CL[(newepa.df$CL==0)] <- NA
# newepa.df$ratio <- (newepa.df$SO4+newepa.df$HCO3)/newepa.df$CL
#
# newepa.df$Hardness <- (newepa.df$CA/40.1 + newepa.df$MG/24.3)*100
# newepa.df$SO4HCO3 <- newepa.df$HCO3+ newepa.df$SO4
# newepa.69 <- subset(newepa.df, ecoregl3 == 69 & PHSTVL>6 &state!="WV"& SO4>0); dim(newepa.69)
#
# newepa.var <- c("HCO3", "SO4","CL","CA","MG", "NA.","K","PHSTVL","ratio","COND","SO4HCO3")
#
# table(factor(newepa.69$dataset))
# table(factor(newepa.69$state))
# #NRiverStreamsA      OR-Rivers     OR-Streams   R10-Cascades
# #             9              1             33             77
# # R10-Deschutes REMAP_R10_WAOR    WA-Chehalis          WEMAP
# #             2              2              1             27
#
# # KY MD PA TN VA
# # 16  8 49 11 11
#
# listout69 <- t(apply(newepa.69[newepa.var], 2, sumstats))
#
# write.csv(listout69, "_Tables_Figures_Final/Table_F04_newepa_69.csv")
#
#
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# xc95 <- read.table(paste(wd, "/Ion_HC05/", "xc95_nowt.txt", sep=""), header = T); dim(xc95)
#
# tolval.e <- (read.table(paste(wd, ecolab, "bootxc95.txt", sep = ""),header=T)); bootxc95 <- t(tolval.e); dim(bootxc95)
#
# sortxc95 <- data.frame(lapply(1:1000,function(x) sort(as.vector(bootxc95[,x]), na.last =TRUE)))    # sort each iteration by xc95, 1000X
# names(sortxc95) <- 1:1000
# conf <- apply(sortxc95, 1, quantile, prob = c(0.05, 0.95), na.rm =T)          # conf for xc95 not for taxa
# 10^(quantile(apply(bootxc95, 2, quantile, prob = 0.05, type = 6, na.rm =T) ,
#              prob = c(0.025, 0.5, 0.975)) )
#
# hcs1 <- round(10^apply(tolval.e, 1, quantile, prob = seq(0.01, 1, length = 100), type=6, na.rm = T))
# hcs.conf <- data.frame(t(round(apply(hcs1, 1, quantile, prob = c(0.025, 0.975)))))
# hcs.conf$prob <-  seq(0.01, 1, length = 100)
#
# ######## a function to calcualte probability of XC95
# pro.fun <- function(x) {
#   x2 <- sort(x)
#   prob.temp <- ((1:length(x2))/(1+length(x2)))[1:100]
#   return(prob.temp)
# }
#
# prob1 <-  apply(sortxc95, 2, pro.fun)   # probability for each XC95
# conf.prob <- apply(prob1, 1, quantile, prob = c(0.05, 0.95))                    # prob confidence intervals
#
# xc_conf <- t(data.frame(apply(bootxc95, 1, quantile, prob = c(0.05,0.95), na.rm =T) ) )    # taxa confidence interval
# xc95.cur <- na.omit(xc95[c("taxaname", paste("Eco", eco3,"_xc95", sep =""))]);(xc95.cur)
# xc95.cur$taxaname <- factor( xc95.cur$taxaname)
# xc95.cur <- merge(xc95.cur, xc_conf, by.x ="taxaname", by.y = 0); dim(xc95.cur)
# xc95.cur <- xc95.cur[order(xc95.cur[paste("Eco", eco3,"_xc95", sep ="")]),]
# # sens.taxa <- xc95.cur$taxaname[1:10]
# xc95.cur$prb <- 1:nrow(xc95.cur)/nrow(xc95.cur)             # taxa probability
#
# ############## ssd taxa confidence limits
# plot.ecdf(log10(xc95.cur[,paste("Eco", eco3,"_xc95", sep ="")]))
# plot.ecdf((xc95.cur[,"5%"]),add =T,  col ="green")
# plot.ecdf((xc95.cur[,"95%"]),add =T,  col ="green")
# ###### end
# ###~~~~~~~~SSD confidence limits ~~~~~~~~~~~~~~~~~~~
# filename <- "Figure_F13_xc95_ci.jpg"
# jpeg(filename = paste(wd, "/_Tables_Figures_Final/", filename, sep=""), width = 400, height = 400,
#      pointsize = 12,  restoreConsole = TRUE)
#
# plot(log10(xc95.cur[,paste("Eco", eco3,"_xc95", sep ="")]), xc95.cur$prb, pch = 20, ylab="Proportion of genera extirpated", cex = 0.4,
#      col = "green", ylim = range(xc95.cur$prb[1:36]), xlim = c(min(xc95.cur[,"5%"], na.rm=T), 3),
#      xlab = unitlab, axes = F, main = "", xaxt = "n")                              # quantile xc95.cur
# for(index in 1:1000) {
#   sampl2 <- sort(bootxc95[,index])
#   prob2 <- 1:length(sampl2)/length(sampl2)
#   points(sampl2, prob2, pch=20, col = "gray", cex = 0.4)
# }
# lines(log10(hcs.conf[2:50,1]), hcs.conf$prob[2:50], lty = 2)          # SSD confidence interval
# lines(log10(hcs.conf[2:50,2]), hcs.conf$prob[2:50], lty = 2)
# axis(1, at = log10(c(100,200, 500, 1000,2000)), labels = c(100,200, 500, 1000,2000))
# axis(1, at = log10(seq(100,1000,by =100)), tcl = -0.3, labels = F)
# axis(2)
# abline(h = 0.05, lty = 3, col = "black")
# box(bty="l")
# fn <- ecdf(log10(xc95.cur[,paste("Eco", eco3,"_xc95", sep ="")]))
# plot(fn, do.points = TRUE, col = "blue", add =T)
# legend("topleft", legend = "", bty ="n")
# dev.off()
#
# ###~~~~~~~~~~~~~~~~~~ Table 3
#
# Table3 <- table(my.data$Level3, my.data$Month)
# write.table(Table3, paste(wd, ecolab, "Table3_monthcount.txt", sep =""))
# write.table(Table3, paste(wd, "/_Tables_Figures_Final/", "Table_F02_monthcount.txt", sep =""))
#
# ##~~~~~~~~~~~  STEP 5 ~~~~~~~~~~~~~~~~~~~~~
#
# ####### step 5.1  background from random sites
# # radom sample from line 251
# random.eco <- subset(random, Level3 == eco3); dim(random.eco)    #539 site and 592 samples for ecoregions
#
# random.site <- aggregate(random.eco["lgSO4HCO3"], list(STATION_ID = random$STATION_ID), mean); dim(random.site)
# n <- nrow(random.site)
# (crit.g <- 10^(qt(0.25, n-1)*sd(random.site$lgSO4HCO3)*sqrt(1+1/n)+mean(random.site$lgSO4HCO3)  ) ) # regional criteria based on group mean
# a1 <- 10^quantile( random.site$lgSO4HCO3, prob = 0.25)   #80
# a2 <- 10^quantile( random.eco$lgSO4HCO3, prob = 0.25)   #80
# qqq <- c(crit.g, a1, a2)
# names(qqq) <- c("regionalCrit","siteQtile", "sampleQtile")
#
# write.table(qqq, paste(wd, ecolab, "Table_random_lgSO4HCO3_background.txt", sep = ""))
#
# ############## calculation confidence interval of background conductivity
# boot.rand <- function(myset, n2 =1000) {
#   n1 <- max(table(myset$Month))
#   Month <- sort(unique(myset$Month))
#   qtile2 <- rep(NA, n2)      # all sample togetther to computer a quantile
#   for(i in 1:n2) {
#     qtile1 <- rep(NA, length(Month))
#     #   array2 <- rep(NA, length(Month)*n1)
#     for(index in 1:length(Month)) {
#       array1 <- sample(rep(myset$lgSO4HCO3[myset$Month == Month[index]],2),n1, replace = T)
#       #     print(Month[index])
#       #     print(round(10^quantile(array1, prob = 0.25)))
#       if(index ==1) {
#         array2 <- array1
#       } else {
#         array2 <- c(array2, array1)
#       }
#     }
#     qtile2[i] <- quantile(array2, prob = 0.25) # all sample quantile
#   }
#   return(qtile2)      # return all sample 25th quantile
# }##FUNCTION.boot.rand.END
#
# tmp <- rep(NA, 1000)
# for(i in 1:1000) {
#   samp <- sample(random.site$lgSO4HCO3, n, replace =T)
#   tmp[i] <- quantile(samp, prob = 0.25, type = 6)
# }##FOR.i.END
# b0 <- round(10^quantile(tmp, prob = c(0.025, 0.05, 0.5, 0.95, 0.975)))
# ### ecoregion 69
# a1 <- boot.rand(random.eco)
# b1 <- round(10^quantile(a1, prob = c(0.025, 0.05, 0.5, 0.95, 0.975)))  # get confidence for quantile
#
# d <- cbind(rbind(b0, b1))
# Orig25tile <- as.vector(round(c(10^(quantile(random.site$lgSO4HCO3, prob = .25)),
#                                 10^quantile(random.eco$lgSO4HCO3, prob = .25)  )   )   )
# d <- cbind(d, Orig25tile)
# write.table (d, paste(wd, ecolab, "Table_random_cond_25th_CI.txt", sep = ""))
#
# # par(mfrow= c(2,2))
# cutv <- c(1, 20, 50, 100, 150, 300, 5000)
# random$cutf <- cut(random$SO4HCO3, cutv)
# # table(random$season, random$cutf)
# monthchart <-  table(random$Level3, random$Month)
# monthmean <- 10^tapply(random.eco$lgSO4HCO3, random.eco$Month, mean)
# jpeg(filename = paste(wd, ecolab, "Figure_Cond_random_monthly_variation.jpg",sep =""), width = 600, height = 600,
#      units = "px", pointsize = 12, quality = 100, bg = "white")
# par(mfrow = c(1,1), mar = c(4,4,4,1))
# boxplot(SO4HCO3~Month, data = random.eco, xlab = "Month",
#         log = "y", ylab = unitlab, col = "lightgray")
# axis(3, at = axTicks(3), labels = as.numeric(monthchart[1,]))
# mtext("Number of Samples", side =3, line = 2.5)
# graphics.off()
# ######## step 5.2    background for reference sites
# ref  <- subset(ref.sub, REFERENCE == "Level I" & !is.na(lgSO4HCO3))[c("Sample.ID", "STATION_ID","SO4HCO3","Level3",
#                                                                       "Year", "Month","lgSO4HCO3")]; dim(ref)
#
# sprg <- ref$Month <=6 & ref$Month >=3             # month range from 4 to 10
# summ <- ref$Month <=10 & ref$Month >=7
# ref$season <- ref$Month
# ref$season[sprg] <- "Spring"
# ref$season[summ] <- "Summer"
# ref$season[!(summ | sprg)] <- "Winter"
#
# q1 <- aggregate(ref["lgSO4HCO3"], list(season = ref$season, Ecoregion = ref$Level3),
#                 quantile, prob = 0.75)
# q3 <- aggregate(ref["lgSO4HCO3"], list(season = rep("allyear", length(ref$season)),
#                                        Ecoregion = ref$Level3), quantile, prob = 0.75)
# count <- table(ref$season, ref$Level3)
# ref.qq <- rbind(q1, q3)
# ref.qq$SO4HCO3 <- round(10^ref.qq$lgSO4HCO3)
# ref.qq$NSample <- c(as.vector(count), as.vector(sum(count)))
# write.table(ref.qq, paste(wd, ecolab, "Table_ref_cond_75th.txt", sep = ""))
#
# ########   estimate confidence intervals for reference sites
# my.ref <- subset(ref, Level3==eco3); dim(my.ref)             # 43
#
# boot.ref <- function(myset) {
#   n1 <- max(table(myset$Month))
#   Month <- unique(myset$Month)
#   mean2 <- rep(NA, 1000)       # get a mean
#   qtile2 <- rep(NA, 1000)      # all sample togetther to computer a quantile
#   qtile3 <- rep(NA,1000)       # only get a quantile for a month and then compute mean
#   for(i in 1:1000) {
#     mean1 <- rep(NA, length(Month))
#     qtile1 <- rep(NA, length(Month))
#     #   array2 <- rep(NA, length(Month)*n1)
#     for(index in 1:length(Month)) {
#       array1 <- sample(rep(myset$lgSO4HCO3[myset$Month == Month[index]],2),n1, replace = T)
#       mean1[index] <- mean(array1)
#       qtile1[index] <- quantile(array1, prob=0.75)
#       if(index ==1) {
#         array2 <- array1
#       } else {
#         array2 <- c(array2, array1)
#       }
#     }
#     mean2[i] <- mean(mean1)     # month mean o
#     qtile2[i] <- quantile(array2, prob = 0.75) # all sample quantile
#     qtile3[i] <- mean(qtile1)          # quantile month mean
#   }
#   return(c(mean2, qtile2, qtile3))      # return mean, all sample 75th quantile, and 75th quantile mean
# }##FUNCTION.boot.ref.END
#
# ### ecoregion 69
# a1 <- boot.ref(my.ref)
# b1 <- round(10^quantile(a1[1:1000], prob = c(0.025, 0.05, 0.5, 0.95, 0.975)))    # get confidence for mean
# b11 <- round(10^quantile(a1[1001:2000], prob = c(0.025, 0.05, 0.5, 0.95, 0.975)))  # get confidence for quantile
# b12 <-round(10^quantile(a1[2001:3000], prob = c(0.025, 0.05, 0.5, 0.95, 0.975)))  # get confidence for quantile mean
# c1<- round(10^mean(a1[2001:3000]))  #     mean mean
#
#
# d <- rbind(b1,b11, b12)
# rownames(d) <- c("mean","allqtile","qtile")
# write.table(d, paste(wd, ecolab, "Table_ref_cond_75th_CI.txt", sep = ""))
#
# ################### step 5.3    background for all sites   25th percentile
# sprg <- my.data$Month <=6 & my.data$Month >=3             # month range from 4 to 10
# summ <- my.data$Month <=10 & my.data$Month >=7
# my.data$season <- my.data$Month
# my.data$season[sprg] <- "Spring"
# my.data$season[summ] <- "Summer"
# my.data$season[!(summ | sprg)] <- "Winter"
#
# q1 <- aggregate(my.data["SO4HCO3"], list(season = my.data$season), quantile, prob = 0.25)
# q3 <- aggregate(my.data["SO4HCO3"], list(season = rep("allyear", length(my.data$season))), quantile, prob = 0.25)
# count <- table(my.data$season)
# all.qq <- rbind(q1, q3)
# all.qq$NSample <- c(as.vector(count), as.vector(sum(count)))
# write.table(all.qq, paste(wd, ecolab, "Table_alldata_cond_25percentile.txt", sep = ""))
#
# back <- rep(NA, 1000)
# for(i in 1:1000) {
#   back[i] <- quantile(sample(my.data$SO4HCO3, nrow(my.data), replace = T), prob = 0.25)
# }
# quantile(back, prob = c(0.025, 0.5, 0.95))

##~~~~~~~~~~~~~~~~~~~~ Step 10 ~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~ Appendix F ~~~~~~~~~~~~~~~~~~~~

# move here from above, 20170406.
# use "bio.sample" that is already created.
#
# bio.sample <- subset(bio.sample, !duplicated(STATION_ID, fromLast=TRUE))    # 3115
my70 <- subset(bio.sample, Level3==70 &Ratio>1); dim(my70)              # 1042
my69 <- subset(bio.sample, Level3==69 & Ratio>1); dim(my69)              # 925
#
switch0 <- 1
#
if(switch0 == 1) my.data <- my69
if(switch0 == 2) my.data <- my70
#
eco3 <- ifelse(switch0 == 1, 69, 70)
#



#xc95 <- read.table(paste(wd, "/Ion_HC05/", "xc95_nowt.txt",sep=""), header = T)
xc95 <- read.table(paste(wd, "/data-raw/", "xc95_nowt.txt",sep=""), header = T)
xc.sel <- xc95[!is.na(xc95[paste("Eco", eco3, "_xc95", sep ="")]),]; dim(xc.sel)
taxalist <- as.vector(xc.sel[order(xc.sel[paste("Eco", eco3, "_xc95", sep ="")]) , "taxaname"])
df1 <- merge(ss, my.data[c("Sample.ID","lgSO4HCO3","Long_DD","Lat_DD")], by = "Sample.ID"); dim(df1)
#graphics.off()
full.results <- taxon.response.sort(df1 = df1, xvar="lgSO4HCO3", cutoff=25, region = ecolab,
                                    mtype = 3, dense.N = 201, plot.pdf = T, xlabs= unitlab, add.map = F,
                                    maintext = "", GIS.cord = c("Long_DD", "Lat_DD"),
                                    log.x=TRUE, rounder=0, taus=c(0,95,100), nbin =61, sort.vect = taxalist)

# full.results <- merge(full.results, ref.taxa.all, by.x = "taxaname", by.y = "Genus"); dim(full.results)
# quantile(full.results$Observ100, 0.05, type =6)
#
# table.tmp <- full.results[order(full.results$XC95),]
# write.table(table.tmp, paste(wd, ecolab, "Table_xc95_all.txt", sep =""))
#
# TableD <- table.tmp[,c("Order", "Family", "taxaname","Trend", "XC95", "N")]
#
# write.table(TableD, paste(wd, "/_Tables_Figures_Final/", "Table_F07_TableD_xc95.csv", sep =""),sep=",",row.names=F)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Save as RDA for use in package####
# df1 and taxalist.
# already have ss
devtools::use_data(df1)
devtools::use_data(taxalist)

