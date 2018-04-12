# Prepare data for example for weightcdf2.R; weightcdf()
#
# From _Appendix_F_all.R
#
# Erik.Leppo@tetratech.com
# 20170405
##~~~~~~~~~~~~~~~~~~~~

# 0. Prep####
# library(devtools)
# library(reshape)
wd <- getwd()

# 1. Get data and process#####

#################### do with the new dataset
#species.all <- read.delim(file = paste(wd,"/data/BenthicCompositionGenera.txt",sep=""), header = T); dim(species.all)
species.all <- read.delim(file = paste(wd,"/data-raw/BenthicCompositionGenera.txt",sep=""), header = T); dim(species.all)

mydate <- strptime(species.all$SAMPLE_DATE, format="%m/%d/%Y")
species.all$Year <- mydate$year + 1900
species.all$Month <- mydate$mon + 1
species.all$Julian <- as.numeric(format(mydate, "%j"))
species.all$NumMonth <- species.all$Month
species.all$myDate <- as.Date(species.all$SAMPLE_DATE, format="%m/%d/%Y")
species.all$Sta_Year <- as.factor(paste(species.all$STATION_CODE, species.all$Year, sep="_"))
newsum <- aggregate(species.all[,"Benthic.Count"], list(Sample.ID = species.all$Sample.ID), sum, na.rm =T)
species.all <- merge(species.all, newsum)
species.all$RA <- species.all$Benthic.Count/species.all$x
#ITIS.genera <- read.delim(paste(wd, "/data/ITIS_genera.txt", sep=""), header = T); dim(ITIS.genera)

#xc95 <- read.table(paste(wd, "/ecoregion/xc95_nowt.txt",sep=""), header = T)
xc95 <- read.table(paste(wd, "/data-raw/xc95_nowt.txt",sep=""), header = T)
xc95 <- xc95[order(xc95$Eco69_xc95),]
sens.taxa69 <- xc95$taxaname[xc95$Eco69_xc95<310]
xc95 <- xc95[order(xc95$Eco70_xc95),]
sens.taxa70 <- xc95$taxaname[xc95$Eco70_xc95<340]
new.sens <- subset(species.all, GENUS %in% sens.taxa69 & Level3 == 69); dim(new.sens)
new.sens$Sta_Year <- factor(new.sens$Sta_Year)
new.sta69 <- aggregate(new.sens$GENUS, list(Sta_Year= new.sens$Sta_Year, Sample.ID= new.sens$Sample.ID), length); dim(new.sta69)  # 2379
names(new.sta69)[3] <- "SensTaxa"
new.sens <- subset(species.all, GENUS %in% sens.taxa70 & Level3 == 70); dim(new.sens)
new.sens$Sta_Year <- factor(new.sens$Sta_Year)
new.sta70 <- aggregate(new.sens$GENUS, list(Sta_Year= new.sens$Sta_Year, Sample.ID= new.sens$Sample.ID), length); dim(new.sta70)   # 2127
names(new.sta70)[3] <- "SensTaxa"
new.sta <- rbind(new.sta69, new.sta70) ; dim(new.sta)

#fullset <- read.delim(file = paste(wd,"/data/Data_combined.txt", sep =""), header = T); dim(fullset)
fullset <- read.delim(file = paste(wd,"/data-raw/Data_combined.txt", sep =""), header = T); dim(fullset)
fullset$RBP_7Sc <- with(fullset, Fish.Cover + Embeddedness + Channel.Alteration + Sediment.Deposition +
                          Total.Bank.Stability + Total.Bank.Vegetation + Total.Undisturbed.Vegetation)
RBP_avg <- aggregate(fullset["RBP_7Sc"], list(fullset$STATION_ID), mean, na.rm =T); dim(RBP_avg)
names(RBP_avg) <- c("STATION_ID", "RBP_7avgSc")
fullset <- merge(fullset, RBP_avg, by = "STATION_ID"); dim(fullset)

fullset$Conductivity[fullset$Conductivity ==0] <- NA
fullset$Alkalinity[fullset$Alkalinity < 1]  <- NA
fullset$Na_Tot[fullset$Na_Tot < 0.1] <- NA
fullset$Fecal[fullset$Fecal==0] <- 0.5               # set 0 to 0.5
fullset$Ca_Tot[fullset$Ca_Tot < 0.5] <- NA
fullset$Mg_Tot[fullset$Mg_Tot < 0.5] <- NA
fullset$Fe_Dis[fullset$Fe_Dis < 0.01] <- NA
fullset$Fe_Tot[fullset$Fe_Tot < 0.01] <- NA
fullset$Al_Dis[fullset$Al_Dis < 0.01] <- NA
fullset$Al_Tot[fullset$Al_Tot < 0.01] <- NA
fullset$TP[fullset$TP==0] <- NA
fullset$Hardness[fullset$Hardness < 1] <- NA
fullset$Chloride[fullset$Chloride_Total_Q =="<"] <- fullset$Chloride[fullset$Chloride_Total_Q=="<"]/2
fullset$Chloride_Total_Q[fullset$Chloride_Total_Q !="<"&fullset$Chloride_Total_Q !=""] <- NA
fullset$Bicarbonate <- fullset$Alkalinity * 61 *2/100
fullset$SO4HCO3 <- (fullset$Sulfate + fullset$Bicarbonate)
fullset$lgSO4HCO3 <- log10(fullset$SO4HCO3)
fullset$SO4CL <- (fullset$Sulfate + fullset$Chloride)
fullset$HCO3_eq <- fullset$Alkalinity/100*2            # equivlent
fullset$SO4_eq <- fullset$Sulfate/96*2                 # equivlent
fullset$HCO3SO4_eq <- (fullset$HCO3_eq + fullset$SO4_eq)      # equivlent
fullset$Cl_eq <- fullset$Chloride/35.45               # equivlent
fullset$SO4_mol <- fullset$Sulfate/96
fullset$HCO3SO4_mol <- (fullset$HCO3_eq + fullset$SO4_mol)
fullset$lgSO4CL <- log10(fullset$SO4CL)
fullset$Anion <-  fullset$SO4HCO3 + fullset$Chloride
fullset$lgAnion <- log10(fullset$Anion)
fullset$Cation <- fullset$Ca_Tot + fullset$Mg_Tot + fullset$Na_Tot + fullset$K_Tot
fullset$CaMg <-  fullset$Ca_Tot + fullset$Mg_Tot
fullset$cond <- log10(fullset$Conductivity)
fullset$lgFecal <- log10(fullset$Fecal)
fullset$Ratio <-  fullset$SO4HCO3/fullset$Chloride
mydate <- strptime(fullset$SAMPLE_DATE, format="%m/%d/%Y")
fullset$Year <- mydate$year + 1900
fullset$Month <- mydate$mon + 1
fullset$Julian <- as.numeric(format(mydate, "%j"))
fullset$NumMonth <- fullset$Month
fullset$Year2 <- fullset$Year
fullset$Year2[fullset$Month >7 & !is.na(fullset$Month)] <-
  fullset$Year2[fullset$Month >7 & !is.na(fullset$Month)] +1 # move last year summer samples to this year summer
fullset$Sta_Year <- as.factor(paste(fullset$STATION_CODE, fullset$Year, sep="_"))
fullset$Sta_Year2 <- as.factor(paste(fullset$STATION_CODE, fullset$Year2, sep="_"))
fullset$Julian2 <- fullset$Julian
fullset$Julian2[fullset$Julian2>182 & !is.na(fullset$Julian2)] <- fullset$Julian2[fullset$Julian2>182&!is.na(fullset$Julian2)]-365
fullset$myDate <- as.Date(fullset$SAMPLE_DATE, format="%m/%d/%Y")
fullset <- merge(fullset, new.sta[c("Sta_Year","SensTaxa")], all.x =T, by.x = "Sta_Year2", by.y = "Sta_Year"); dim(fullset)
fullset$SensTaxa[is.na(fullset$SensTaxa)] <- 0

################## biological samples with conductivity
envbio <- subset(fullset, !is.na(PCT_OF_THRESHOLD_GLIMPSS_CF) & !is.na(cond) & !is.na(Month) & !duplicated(Sample.ID)&
                   (Level3 == 69 | Level3 == 70)); dim(envbio)

biosta <- subset(envbio, !duplicated(STATION_ID)); dim(biosta)
bio.sta <- aggregate(envbio[c("lgSO4HCO3", "WVSCI", "PCT_OF_THRESHOLD_GLIMPSS_CF")], #  "SensTaxa")],
                     list(Sta_Year = envbio$Sta_Year, Sample.ID = envbio$Sample.ID), mean, na.rm =T); dim(bio.sta)
bio.date <- aggregate(envbio["Month"], list(Sta_Year= envbio$Sta_Year), min, na.rm =T)  # only one biosample in a year
bio.sta <- merge(bio.sta, bio.date, all.x =T)

# ################# chemistry samples in the two ecoregions
chemistry <- subset(fullset, (Level3 == 69 | Level3 == 70)&!is.na(Conductivity)); dim(chemistry)
chem.sta <- subset(chemistry, !duplicated(STATION_ID)); dim(chem.sta)
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
####################### Define biological samples then weighting the conductivity values
bio.sample <- subset(envbio, pH>6 & (Ratio > 1|is.na(Ratio))); dim(bio.sample)      # 3736
sum(bio.sample$Level3==69); sum(bio.sample$Level3==70)          ###     1661 and 2075

# #   write.table(bio.sample, paste(wd,"/IntermediaFile/bio.sample.csv", sep =""), sep =",", row.names =T)
# bio.sample <- bio.sample[order(bio.sample$STATION_ID, bio.sample$Year),]
# # bio.sample <- subset(bio.sample, !duplicated(STATION_ID, fromLast=TRUE))    # 3115
# my70 <- subset(bio.sample, Level3==70 &Ratio>1); dim(my70)              # 1042
# my69 <- subset(bio.sample, Level3==69 & Ratio>1); dim(my69)              # 925
# sum(table(my69$STATION_ID)>1)

#~~~~~~~~~~~~~~~~~~~~~~
# ##~~~~~~~~~~~
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
# ###~~~~~~~~~~~~

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ecoregion begins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ref.site <- subset(bio.sample, REFERENCE == "Level I"); dim(ref.site)      # 133

ref.taxa <- unique((subset(species.all, GENUS !="" & Sample.ID %in% ref.site$Sample.ID))$GENUS) ; length(ref.taxa)
ref.taxa69 <- unique((subset(species.all, GENUS !="" & Sample.ID %in% ref.site$Sample.ID & Level3 == 69 ))$GENUS) ; length(ref.taxa69) # 193
ref.taxa70 <- unique((subset(species.all, GENUS !="" & Sample.ID %in% ref.site$Sample.ID & Level3 == 70 ))$GENUS) ; length(ref.taxa70) # 179
#ref.taxa.all <- ITIS.genera[ITIS.genera$Genus %in% ref.taxa,2:4]
ref.fam <- unique((subset(species.all, FAMILY !="" & Sample.ID %in% ref.site$Sample.ID))$FAMILY) ; length(ref.fam)
all.taxa <- unique((subset(species.all, GENUS !="" & Sample.ID %in% bio.sample$Sample.ID))$GENUS); length(all.taxa) #476
ref.genera <- as.vector(ref.taxa[regexpr("/", ref.taxa)==-1]); length(ref.genera)       # 224

species <- subset(species.all, GENUS %in% ref.genera & Sample.ID %in% envbio$Sample.ID &
                    !is.na(RA)); dim(species)
species <- merge(chem.sta[c("STATION_ID", "STATION_CODE")], species); dim(species)
ss <- reshape::cast(species, Sample.ID ~ GENUS, sum, value = "RA");dim(ss)

##### ecoregion 70 and 69 have similar taxa

# taxa69 <- sort(unique(species[species$Level3==69,"GENUS"]))    # 227
#
# unitlab <- expression(paste("SO"[4]^{2-phantom()}," + HCO"[3]^{-phantom()}, " (mg/L)"))
#
# #### plot two ecoregion seperately
# my69.ref <- subset(ref.site, Level3==69); dim(my69.ref);sum(table(my69.ref$STATION_ID)>0)  # 87 / 64               # 87 samples,   64 stations
# my69.spr <- subset(my69, Month > 2 & Month < 7); dim(my69.spr)           # 627
# my69.sum <- subset(my69, Month > 6 & Month < 11); dim(my69.sum)          # 1016
# length(unique(my69.ref$STATION_ID))


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Save as RDA for use in package#####
# ss and bio.sample
devtools::use_data(ss)
devtools::use_data(bio.sample)

