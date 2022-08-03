# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("ggplot2")
install.packages("cowplot")
install.packages("dplyr")
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/FRAP data"
setwd(mainDir)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/FRAP data/NDECDmScarlet truncation constructs CLEANED"
# Set the number of prebleach frames
PreBleach <- 10

#### Data analysis functions ####
## Calculate the standard error for each row of a data frame 
SEMrows <- function(df){
  dft <- as.data.frame(t(df))
  rm(df)
  dft <- dft[-1,]
  dft <- na.omit(dft)
  output <- matrix(ncol = length(colnames(dft)), nrow = 1)
  for (n in 1:length(colnames(dft))) {
    dftn <- dft[[n]]
    output[,n] <- sqrt(var(dftn)/length(dftn))
    rm(dftn)
  }
  output <- as.data.frame(t(output))
  return(output$V1)
}
#### end ####

#### Models from PRISM as functions ####
WTmodelOne <- function(dft){
  ModelData <- 0.9201*(1-exp(-0.09454*dft))
  return(ModelData)
}
PESTmodelOne <- function(dft){
  ModelData <- 0.9202*(1-exp(-0.08637*dft))
  return(ModelData)
}
TADmodelOne <- function(dft){
  ModelData <- 0.9680*(1-exp(-0.09228*dft))
  return(ModelData)
}
TADPESTmodelOne <- function(dft){
  ModelData <- 0.9387*(1-exp(-0.1266*dft))
  return(ModelData)
}
OPAmodelOne <- function(dft){
  ModelData <- 0.9343*(1-exp(-0.1058*dft))
  return(ModelData)
}

# Two phase model where Y0 is set to 0 = 
# (Plateau*PercentFast*.01)*(1-exp(-KFast*X)) + (Plateau*(100-PercentFast)*.01)*(1-exp(-KSlow*X))
WTmodelTwo <- function(dft){
  ModelData <- (0.9668*40.96*.01)*(1-exp(-0.3027*dft)) + (0.9668*(100-40.96)*.01)*(1-exp(-0.04577*dft))
  return(ModelData)
}
PESTmodelTwo <- function(dft){
  ModelData <- (0.9588*37.13*.01)*(1-exp(-0.2822*dft)) + (0.9588*(100-37.13)*.01)*(1-exp(-0.04646*dft))
  return(ModelData)
}
TADmodelTwo <- function(dft){
  ModelData <- (1.027*42.75*.01)*(1-exp(-0.3174*dft)) + (1.027*(100-42.75)*.01)*(1-exp(-0.03994*dft))
  return(ModelData)
}
TADPESTmodelTwo <- function(dft){
  ModelData <- (0.9806*35.35*.01)*(1-exp(-0.4941*dft)) + (0.9806*(100-35.35)*.01)*(1-exp(-0.07033*dft))
  return(ModelData)
}
OPAmodelTwo <- function(dft){
  ModelData <- (0.9936*45.18*.01)*(1-exp(-0.3332*dft)) + (0.9936*(100-45.18)*.01)*(1-exp(-0.04432*dft))
  return(ModelData)
}
#### end ####

#### Read in the data ####
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
## Read in all the csv files and calculate the double normalisation value (DNV) and normalised DNV
RawData <- list()
for (f in 1:length(FolderNames)){
  # Create a path for each subfolder
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  # List the files in each subfolder
  Files <- list.files(Path, pattern = "*.csv")
  SubFiles <- list()
  for (i in 1:length(Files)){
    # Create the filepath within the subfolder
    FilePath <- paste(Path,Files[i], sep ="/")
    # Read in one experiment
    df <- read.csv(FilePath, skipNul = TRUE, fileEncoding="latin1", skip = 1)
    rm(FilePath)
    # Rename the ROI columns to be descriptive for clarity
    df <- rename(df, BleachROI = ROI.02... , WholeNucROI = ROI.03... , BackgroundROI = ROI.04... , InnerNucROI = ROI.05...)
    # Calculate the mean of the background throughout the experiment
    BgVal <- mean(df$BackgroundROI)
    # Calculate the mean signal in BleachROI in the frames before bleaching
    PreBleachValues <- slice_head(df, n = PreBleach)
    InitialPreBleach <- mean(PreBleachValues$BleachROI)
    # Calculate the mean signal in NucleusROI in the frames before bleaching
    TotalNucPreBleach <- mean(PreBleachValues$WholeNucROI)
    rm(PreBleachValues)
    # Calculate the double normalisation values
    InitialPreBleachSubBg <- InitialPreBleach - BgVal
    TotalNucPreBleachSubBg <- TotalNucPreBleach - BgVal
    rm(InitialPreBleach, TotalNucPreBleach)
    df$DNV <- NULL
    for (d in 1:nrow(df)){
      BleachNorm <- (df[d,"BleachROI"]-BgVal)*TotalNucPreBleachSubBg
      TotalNucNorm <- (df[d,"WholeNucROI"]-BgVal)*InitialPreBleachSubBg
      df[d,"DNV"] <- BleachNorm/TotalNucNorm
      rm(BleachNorm,TotalNucNorm)
    }
    rm(InitialPreBleachSubBg,TotalNucPreBleachSubBg,BgVal)
    # Normalise the DNV to be from 0-1, here assuming the max DNV value is 1
    BleachedDNV = df[PreBleach+1,"DNV"]
    df$DNVNorm <- NULL
    for (d in 1:nrow(df)){
      df[d,"DNVNorm"] <- (df[d,"DNV"]-BleachedDNV)/(1-BleachedDNV)
    }
    df$Nucleus <- paste(parse_number(FolderNames[f]), Files[i], sep = "_")
    df$Genotype <- str_sub(basename(Files[[i]]), end=-10)
    SubFiles[[i]] <- df
    names(SubFiles)[i] <- paste(parse_number(FolderNames[f]), Files[i], sep = "_")
    rm(df,BleachedDNV)
  } 
  RawData[[f]] <- SubFiles
  rm(SubFiles, Path, Files)
  names(RawData)[f] <- FolderNames[f]
}
rm(FolderNames)
AllData <- RawData[[1]]
for (k in 2:length(RawData)){
  Ls <- RawData[[k]]
  AllData <- c(AllData, Ls)
  rm(Ls)
}
rm(RawData)

## Change the time points to all be the same (there is some very slight variation in each experiment)
AllTimeScale <- as.data.frame(AllData[[1]]$Axis..s.)
for (t in 2:length(AllData)){
  NextCol <- as.data.frame(AllData[[t]]$Axis..s.)
  AllTimeScale <- cbind(AllTimeScale, NextCol)
}
MeanTimeScale <- rowMeans(AllTimeScale)
for (t in 1:length(AllData)){
  AllData[[t]]$Axis..s. <- MeanTimeScale
  names(AllData[[t]])[names(AllData[[t]]) == "Axis..s."] <- "Time (s)"
}
rm(NextCol,AllTimeScale,MeanTimeScale)

## Organise the data into different genotypes
AllDataDF <- do.call(rbind, AllData)
Genotypes <- unique(AllDataDF$Genotype)
OrganisedAllData <- list()
for (g in 1:length(Genotypes)){
  df <- filter(AllDataDF,Genotype == Genotypes[g])
  OrganisedAllData[[g]] <- df
  names(OrganisedAllData)[g] <- Genotypes[g]
  rm(df)
}
rm(AllDataDF)

## Keep only the DNVNorm values organised by Genotype
Data <- list()
DataMinusPreBleach <- list()
for (g in 1:length(OrganisedAllData)){
  df <- OrganisedAllData[[g]]
  df <- subset(df, select = c(`Time (s)`,`DNVNorm`,`Nucleus`,`Genotype`))
  df <- spread(df,`Nucleus`,`DNVNorm`)
  Data[[g]] <- df
  names(Data)[g] <- names(OrganisedAllData)[g]
  df <- df[-(1:PreBleach), , drop = FALSE]
  df$`Time (s)` <- df$`Time (s)` - df$`Time (s)`[1]
  DataMinusPreBleach[[g]] <- df
  names(DataMinusPreBleach)[g] <- names(OrganisedAllData)[g]
  rm(df)
}

## Save the DNV files for each genotype
for (d in 1:length(Genotypes)) {
  df <- Data[[d]]
  filename <- paste(Genotypes[d],".csv")
  write.csv(df,filename)
  rm(filename)
}

## Calculate the model data
ModelDataOne <- as.data.frame(DataMinusPreBleach[["WT"]]$`Time (s)`)
names(ModelDataOne)[names(ModelDataOne) == "DataMinusPreBleach[[\"WT\"]]$`Time (s)`"] <- "Time (s)"
ModelDataOne$WT <- WTmodelOne(ModelDataOne$`Time (s)`)
ModelDataOne$PEST <- PESTmodelOne(ModelDataOne$`Time (s)`)
ModelDataOne$TAD <- TADmodelOne(ModelDataOne$`Time (s)`)
ModelDataOne$TADPEST <- TADPESTmodelOne(ModelDataOne$`Time (s)`)
ModelDataOne$OPA <- OPAmodelOne(ModelDataOne$`Time (s)`)

ModelDataTwo <- as.data.frame(DataMinusPreBleach[["WT"]]$`Time (s)`)
names(ModelDataTwo)[names(ModelDataTwo) == "DataMinusPreBleach[[\"WT\"]]$`Time (s)`"] <- "Time (s)"
ModelDataTwo$WT <- WTmodelTwo(ModelDataTwo$`Time (s)`)
ModelDataTwo$PEST <- PESTmodelTwo(ModelDataTwo$`Time (s)`)
ModelDataTwo$TAD <- TADmodelTwo(ModelDataTwo$`Time (s)`)
ModelDataTwo$TADPEST <- TADPESTmodelTwo(ModelDataTwo$`Time (s)`)
ModelDataTwo$OPA <- OPAmodelTwo(ModelDataTwo$`Time (s)`)

## Calculate the mean and SEM for each genotypes timepoints
MeanData <- NULL
for (d in 1:length(Data)){
  df <- Data[[d]]
  Yvalues <- subset(df, select = -c(`Time (s)`,`Genotype`))
  df <- subset(df, select = c(`Time (s)`,`Genotype`))
  df$MeanFluorescence <- rowMeans(Yvalues)
  df$FluorescenceSEM <- SEMrows(Yvalues)
  df$`NΔECD construct` <- paste(df$Genotype," (N = ", ncol(Yvalues), ")", sep = "")
  MeanData <- rbind(MeanData,df)
  rm(df,Yvalues)
}
MeanDataMinusPreBleach <- NULL
for (d in 1:length(DataMinusPreBleach)){
  df <- DataMinusPreBleach[[d]]
  Yvalues <- subset(df, select = -c(`Time (s)`,`Genotype`))
  df <- subset(df, select = c(`Time (s)`,`Genotype`))
  df$MeanFluorescence <- rowMeans(Yvalues)
  df$FluorescenceSEM <- SEMrows(Yvalues)
  df$`NΔECD construct` <- paste(df$Genotype," (N = ", ncol(Yvalues), ")", sep = "")
  df$ModelOneValues <- unlist(ModelDataOne[unique(df$Genotype)])
  df$ModelTwoValues <- unlist(ModelDataTwo[unique(df$Genotype)])
  MeanDataMinusPreBleach <- rbind(MeanDataMinusPreBleach,df)
  rm(df,Yvalues)
}
#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  plot.title = element_text(hjust=0.5, face = 'bold'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  text = element_text(size=15), # make text larger
  axis.line = element_line(colour = "black"), # add axis lines in black
  aspect.ratio = 2/3) # set the relative length of the y/x axis
}

## Where required, labels are the Title, X axis and Y axis labels defined in a character string (in order)
## e.g. Labels <- c(title = "Normalised to Background value", x = "Distance (um)", y = "Normalised fluorescence")
## these are then read into the function.

## Define a function that plots a dataframe of time, mean, sem, grouping by the Genotype

recovery_plot_OnePhase <- function(Data, Labels){
  ggplot(Data, aes(`Time (s)`,MeanFluorescence,group = `NΔECD construct`,col=`NΔECD construct`,shape=`NΔECD construct`, linetype = Genotype,
                   ymin=MeanFluorescence-FluorescenceSEM,ymax=MeanFluorescence+FluorescenceSEM, ylim(0,1.2))) + 
    theme_jmt() + labs(linetype = "Model") + labs(col = "Data") + labs(shape = "Data") + 
    expand_limits(y = 1.2) + scale_y_continuous(breaks = seq(0, 1.2, 0.2)) +
    geom_errorbar(linetype = "solid",alpha = 0.75, size = 0.25) + geom_point(size = 0.8) + 
    geom_line(aes(y=ModelOneValues), size = 0.4, colour = "black") +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    scale_color_brewer(palette = "Set2")
}
recovery_plot_TwoPhase <- function(Data, Labels){
  ggplot(Data, aes(`Time (s)`,MeanFluorescence,group = `NΔECD construct`,col=`NΔECD construct`,shape=`NΔECD construct`, linetype = Genotype,
                   ymin=MeanFluorescence-FluorescenceSEM,ymax=MeanFluorescence+FluorescenceSEM, ylim(0,1.2))) + 
    theme_jmt() + labs(linetype = "Model") + labs(col = "Data") + labs(shape = "Data") +
    expand_limits(y = 1.2) + scale_y_continuous(breaks = seq(0, 1.2, 0.2)) +
    geom_errorbar(linetype = "solid",alpha = 0.75, size = 0.25) + geom_point(size = 0.8) + 
    geom_line(aes(y=ModelTwoValues), size = 0.4, colour = "black") +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    scale_color_brewer(palette = "Set2")
}
#### end ####

#### Saved plots ####

## One phase model
subDir <- "One phase model plots"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "TADPEST")
recovery_plot_OnePhase(df,Labs) + geom_line(aes(y=ModelOneValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔTADPEST.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "PEST")
recovery_plot_OnePhase(df,Labs) + geom_line(aes(y=ModelOneValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔPEST.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "TAD")
recovery_plot_OnePhase(df,Labs) + geom_line(aes(y=ModelOneValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔTAD.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "OPA")
recovery_plot_OnePhase(df,Labs) + geom_line(aes(y=ModelOneValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔOPA.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

setwd(mainDir)

## Two phase model
subDir <- "Two phase model plots"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "TADPEST")
recovery_plot_TwoPhase(df,Labs) + geom_line(aes(y=ModelTwoValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔTADPEST.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "PEST")
recovery_plot_TwoPhase(df,Labs) + geom_line(aes(y=ModelTwoValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔPEST.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "TAD")
recovery_plot_TwoPhase(df,Labs) + geom_line(aes(y=ModelTwoValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔTAD.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

Labs <- c(title = "", x = "Time (s)", y = "Normalised Su(H)::EGFP\nfluorescence")
df <- subset(MeanDataMinusPreBleach, Genotype == "WT" | Genotype == "OPA")
recovery_plot_TwoPhase(df,Labs) + geom_line(aes(y=ModelTwoValues), size = 0.4, colour = "black")
ggsave("Wild type vs ΔOPA.jpg", device = "jpeg", dpi = "retina",
       width = 20, height = 10, units = "cm")

setwd(mainDir)
#### End ####


