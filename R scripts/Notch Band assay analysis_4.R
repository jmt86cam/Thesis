# clear the environment: 
rm(list = ls())
# clear all the plots: 
dev.off()
# clear the console: ctrl+L

#### Packages for this script ####
install.packages("ggplot2")
install.packages("cowplot")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggsignif")
install.packages("ggpubr")
install.packages("gridExtra")
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggsignif)
library(ggpubr)
library(gridExtra)
#### end ####

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

MeanDF <- function(DataFrame){
  MeanList <- list()
  for (l in 1:length(DataFrame)) {
    df <- DataFrame[[l]]
    dfmeanse <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(df[,-1], na.rm = TRUE))
    dfmeanse$SEM <- SEMrows(df)
    Nnumber <- ncol(df)-1
    dfmeanse$Genotype <- paste(names(DataFrame)[l]," (N = ", Nnumber, ")", sep = "")
    MeanList[[l]] <- dfmeanse
    names(MeanList)[l] <- names(DataFrame[l])
    rm(dfmeanse,Nnumber)
  }
  Output <- dplyr::bind_rows(MeanList, .id = 'Name')
  rm(MeanList)
  return(Output)
}

NonListMeanDF <- function(df, Genotype){
  dfmeanse <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(df[,-1], na.rm = TRUE))
  dfmeanse$SEM <- SEMrows(df[,-1])
  Nnumber <- ncol(df)-1
  dfmeanse$Genotype <- paste(Genotype," (N = ", Nnumber, ")", sep = "")
  rm(Nnumber)
  return(dfmeanse)
}

# Create a list of dataframes, named for each truncation, where the data contained
# is the wild type data collected alongside that truncation only.
SepDates <- function(DataList,GenotypeToSep){
  Dates <- list()
  for (g in 1:length(DataList)) {
    df <- DataList[[g]]
    ColNames <- colnames(df)
    ColNames <- ColNames[-1]
    for (c in 1:length(ColNames)) {
      SplitName <- str_split(ColNames[c]," ")
      SplitName <- SplitName[[1]]
      ColNames[c] <- SplitName[1]
    }
    dfDates <- unique(ColNames)
    Dates[[g]] <- dfDates
    names(Dates)[g] <- names(DataList[g])
    rm(dfDates, ColNames, SplitName)
  }
  
  data <- DataList[[GenotypeToSep]]
  dataD <- as.data.frame(data[,1])
  names(dataD)[names(dataD) == "data[, 1]"] <- names(data)[1]
  
  SubDates <- Dates
  SubDates[[GenotypeToSep]] <- NULL
  
  Output <- list()
  for (c in 1:length(SubDates)) {
    GenotypeDates <- SubDates[[c]]
    dataG <- dataD
    for (d in 1:length(GenotypeDates)) {
      df <- data[,grepl(GenotypeDates[d], names(data))]
      dataG <- cbind(dataG, df)
      rm(df)
    }
    Output[[c]] <- dataG
    names(Output)[c] <- names(SubDates[c])
    rm(dataG,GenotypeDates)
  }
  rm(data, dataD, SubDates, Dates)
  return(Output)
}

# Create a list of dataframes, named for each truncation, where each dataframe
# has the data of the truncation and wild type data from the corresponding days
# organised into the means and SEM for all the bands imaged.
WTvsMutByDate <- function(SeparatedWTdata, AllData){
  OutputList <- list()
  for (l in 1:length(SeparatedWTdata)) {
    Genotype <- names(SeparatedWTdata)[l]
    WT <- NonListMeanDF(SeparatedWTdata[[l]],"Wild Type")
    Mutant <- NonListMeanDF(AllData[[Genotype]],Genotype)
    df <- rbind(WT,Mutant)
    OutputList[[l]] <- df
    names(OutputList)[l] <- names(SeparatedWTdata[l])
    rm(df, WT, Mutant)
  }
  return(OutputList)
}

# This function will take the RelInt values and separate them into a list of
# dateframes for each truncation. The dataframes are named for their mutant
# and contain RelInt values for the truncation and wild type data taken alongside.
RelIntSepDates <- function(RelIntDataFrame,GenotypeToSep) {
  Mutants <- dplyr::filter(RelIntDataFrame, Genotype != GenotypeToSep)
  WildTypes <- dplyr::filter(RelIntDataFrame, Genotype == GenotypeToSep)
  OutputList <- list()
  for (t in 1:length(unique(Mutants$Genotype))) {
    MutantData <- dplyr::filter(RelIntDataFrame, Genotype == unique(Mutants$Genotype)[t])
    StrSplit <- str_split(MutantData$Nucleus, " ")
    Dates <- c()
    for (s in 1:length(StrSplit)) {
      SplitName <- StrSplit[[s]]
      Date <- SplitName[1]
      Dates <- append(Dates,Date)
      rm(Date, SplitName)
    }
    UniqueDates <- unique(Dates)
    rm(Dates)
    df <- WildTypes[grep(UniqueDates[1], WildTypes$Nucleus),]
    for (d in 2:length(UniqueDates)) {
      Nextdf <- WildTypes[grep(UniqueDates[d], WildTypes$Nucleus),]
      df <- rbind(df,Nextdf)
    }
    df <- rbind(df,MutantData)
    OutputList[[t]] <- df
    names(OutputList)[t] <- unique(Mutants$Genotype)[t]
    rm(df)
  }
  return(OutputList)
}
#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  text = element_text(family = "sans", size = 15), # set default font to be Arial
  plot.title = element_text(hjust=0.5, face = 'bold'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.line = element_line(colour = "black")) # add axis lines in black
}

## Where required, labels are the Title, X axis and Y axis labels defined in a character string (in order)
## e.g. Labels <- c(title = "Normalised to Background value", x = "Distance (um)", y = "Normalised fluorescence")
## these are then read into the function.

## Define a function that plots a dataframe of distance, mean, sem, grouping by the IncubationTime
## as a fourth column
band_plot <- function(Data, Labels){
  ggplot(Data, aes(Distance,Mean,group = Genotype,col=Genotype,shape=Genotype,ymin=Mean-SEM,ymax=Mean+SEM)) + theme_jmt() +
    geom_ribbon(linetype = "dashed", fill = "grey90",alpha = 0.5, size = 0.25) + geom_line(size = 0.5) + geom_point() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    scale_color_brewer(palette = "Set2")
}

## Define a function to plot the relative intensities of the bands
RelInt_boxviolin <- function(Data, Labels){
  ggplot(Data, aes(Genotype,RelativeBandIntensity,col=Genotype,shape=Genotype)) + theme_jmt() +
    geom_violin() + geom_boxplot(width=0.1, show.legend = FALSE,colour = "Grey60",alpha = 0.5, outlier.shape = NA) + geom_jitter() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

RelInt_crossbar <- function(Data, Labels){
  ggplot(Data,aes(Genotype,RelativeBandIntensity,col=Genotype,shape=Genotype)) + theme_jmt() +
    geom_jitter(position = position_jitterdodge()) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/NDECD"
setwd(mainDir)
# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/NDECD/Background NDECD values.csv", head = TRUE, sep = ",")
BackgroundValues <- data.frame(Nucleus=DF[,1], NDECDBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/NDECD/Cleaned band assay csv files"

#### Read in and normalise the data ####
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
NDECDRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
  } 
  rm(Path)
  rm(Files)
  df <- read.csv(FilePaths[[1]])
  names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[1]]), end=-11)
  for (i in 2:length(FilePaths)){
    Nextdf <- read.csv(FilePaths[[i]])
    df <- merge(df,Nextdf, by="X", all = TRUE)
    names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[i]]), end=-11)
    rm(Nextdf)
  }
  names(df)[names(df) == 'X'] <- "Distance (um)"
  NDECDRawData[[f]] <- df
  rm(df)
  rm(FilePaths)
  names(NDECDRawData)[f] <- FolderNames[f] 
}

rm(FolderNames)

NDECDminsixSub <- list()
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] - NormValue
    rm(NormValue)
  }
  NDECDminsixSub[[s]] <- df
  names(NDECDminsixSub)[s] <- names(NDECDRawData[s])
  rm(df,ColNames)
}

NDECDBgNorm <- list()
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    NormValue <- BackgroundValues$NDECDBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  NDECDBgNorm[[s]] <- df
  names(NDECDBgNorm)[s] <- names(NDECDRawData[s])
  rm(df,ColNames)
}

MeanNDECDRawData <- MeanDF(NDECDRawData)
MeanNDECDBgNorm <- MeanDF(NDECDBgNorm)
MeanNDECDMinSixSub <- MeanDF(NDECDminsixSub)
#### end ####

#### Select only data from the same days for comparing genotypes ####
WTNDECDRDSepDates <- SepDates(NDECDRawData,"Wild Type")
WTNDECDminsixSubSepDates <- SepDates(NDECDminsixSub,"Wild Type")

MeanGenotypesByDateRawData <- WTvsMutByDate(WTNDECDRDSepDates, NDECDRawData)
MeanGenotypesByDateMinSix <- WTvsMutByDate(WTNDECDminsixSubSepDates, NDECDminsixSub)

#### end ####

#### Calculate fold increase of band over background ####
## Here we will use the mean of the 5 middle values and divide by the mean of the 5 values at
## each edge of the band (10 values in total).
RelIntNDECD <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  dfOrdered <- NULL
  for (c in 1:length(colnames(df))) {
    Column <- as.data.frame(df[,c])
    MinSix <- Column[order(Column$`df[, c]`),][1:6]
    dfOrdered <- cbind(dfOrdered,MinSix)
  }
  Background <- colMeans(dfOrdered)
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(NDECDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntNDECD <- rbind(RelIntNDECD,RelInt)
  rm(RelInt, Band, Background, df,dfOrdered,c,Column,MinSix)
}

RelIntNDECDSepDates <- (RelIntSepDates(RelIntNDECD,"Wild Type"))

RelIntNDECDEnds <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(NDECDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntNDECDEnds <- rbind(RelIntNDECDEnds,RelInt)
  rm(RelInt, Band, Background, df)
}

RelIntNDECDEndsSepDates <- (RelIntSepDates(RelIntNDECDEnds,"Wild Type"))

RelIntNDECDBg <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$NDECDBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(NDECDRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background)
    RelIntNDECDBg <- rbind(RelIntNDECDBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntNDECDBgSepDates <- RelIntSepDates(RelIntNDECDBg, "Wild Type")

RelIntNDECDBg0 <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$NDECDBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    BandSubBg <- Band - Background
    RelInt <- as.data.frame(BandSubBg/Band)
    RelInt$Genotype <- names(NDECDRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `BandSubBg/Band`)
    rm(Band, Background, BandSubBg)
    RelIntNDECDBg0 <- rbind(RelIntNDECDBg0,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntNDECDBg0SepDates <- RelIntSepDates(RelIntNDECDBg0, "Wild Type")

RelIntNDECDBgSub <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$NDECDBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band-Background)
    RelInt$Genotype <- names(NDECDRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
    rm(Band, Background)
    RelIntNDECDBgSub <- rbind(RelIntNDECDBgSub,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntNDECDBgSubSepDates <- RelIntSepDates(RelIntNDECDBgSub, "Wild Type")

RelIntNDECDMinSixSub <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  dfOrdered <- NULL
  for (c in 1:length(colnames(df))) {
    Column <- as.data.frame(df[,c])
    MinSix <- Column[order(Column$`df[, c]`),][1:6]
    dfOrdered <- cbind(dfOrdered,MinSix)
  }
  Background <- colMeans(dfOrdered)
  RelInt <- as.data.frame(Band-Background)
  RelInt$Genotype <- names(NDECDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
  RelIntNDECDMinSixSub <- rbind(RelIntNDECDMinSixSub,RelInt)
  rm(RelInt, Band, Background, df,dfOrdered,c,Column,MinSix)
}

RelIntNDECDMinSixSubSepDates <- RelIntSepDates(RelIntNDECDMinSixSub, "Wild Type")

RelIntChangeNDECDBg <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$NDECDBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(NDECDRawData[s])
    RelInt$Nucleus <- rownames(RelInt)
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeNDECDBg <- rbind(RelIntChangeNDECDBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntChangeNDECDBgSepDates <- RelIntSepDates(RelIntChangeNDECDBg, "Wild Type")

RelIntChangeNDECD <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  dfOrdered <- NULL
  for (c in 1:length(colnames(df))) {
    Column <- as.data.frame(df[,c])
    MinSix <- Column[order(Column$`df[, c]`),][1:6]
    dfOrdered <- cbind(dfOrdered,MinSix)
  }
  Background <- colMeans(dfOrdered)
  RelInt <- as.data.frame((Band-Background)/Background)
  RelInt$Genotype <- names(NDECDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeNDECD <- rbind(RelIntChangeNDECD,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
RelIntChangeNDECDSepDates <- (RelIntSepDates(RelIntChangeNDECD, "Wild Type"))
#### end ####

#### Statistics Tests and Plots ####
subDir <- "Stats"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

## The first thing we look for is normality of the data. Note that anything with
## a sample size greater than 30 can be considered normal because of the central
## limit theorem.

## We can do density plots of the band intensities to make sure the data is
## normally distributed

# This function will create a JPEG file of data density on y axis, it will need
# a list of data, the name of the variable density we are looking at, something
# to separate the plots by (e.g. Genotype), and the name of the JPEG file
NormDis_plots <- function(DataList, VariableDens, SepBy, MainTitle){
  Plots <- list()
  for (p in 1:length(DataList)) {
    df <- DataList[[p]]
    Plot <- ggdensity(df, x = VariableDens, facet.by = SepBy)
    Plots[[p]] <- Plot
    rm(Plot)
  }
  ggsave(paste(MainTitle,".jpg"), arrangeGrob(grobs = Plots, ncol = 2),
         device = "jpeg", dpi = "retina", width = 24,
         height = 9*ceiling(length(Plots)/2), units = "cm")
  rm(Plots)
}

NormDis_plots(RelIntNDECDSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelInt plots")
NormDis_plots(RelIntNDECDBgSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntBg plots")
NormDis_plots(RelIntNDECDBg0SepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntBg0 plots")
NormDis_plots(RelIntNDECDBgSubSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntBgSub plots")
NormDis_plots(RelIntNDECDMinSixSubSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntMinSixSub plots")
NormDis_plots(RelIntChangeNDECDSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntChange plots")
NormDis_plots(RelIntChangeNDECDBgSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntChangeBg plots")

## A quantile-quantile plot can also be used to check for normal distribution
## the data should align with the 45 degree line if it is normal.

# This function will create a JPEG file of QQ plots with sample on y axis, it will need
# a list of data, the name of the variable we are looking at, something
# to separate the plots by (e.g. Genotype), and the name of the JPEG file
QQ_plots <- function(DataList, VariableDens, SepBy, MainTitle){
  Plots <- list()
  for (p in 1:length(DataList)) {
    df <- DataList[[p]]
    Plot <- ggqqplot(df, x = VariableDens, facet.by = SepBy)
    Plots[[p]] <- Plot
    rm(Plot)
  }
  ggsave(paste(MainTitle,".jpg"), arrangeGrob(grobs = Plots, ncol = 2),
         device = "jpeg", dpi = "retina", width = 24,
         height = 9*ceiling(length(Plots)/2), units = "cm")
  rm(Plots)
}

QQ_plots(RelIntNDECDSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelInt")
QQ_plots(RelIntNDECDBgSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBg")
QQ_plots(RelIntNDECDBg0SepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBg0")
QQ_plots(RelIntNDECDBgSubSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBgSub")
QQ_plots(RelIntNDECDMinSixSubSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntMinSixSub")
QQ_plots(RelIntChangeNDECDSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntChange")
QQ_plots(RelIntChangeNDECDBgSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntChangeBg")

## Perform the Shapiro-Wilk test of normality on samples. In this test the null
## hypothesis is that the data is normal and therefore if the p value <0.05 the 
## data should be considered to NOT be normal. The test works best on small sample
## sizes and should be considered alongside the graphs.

# This function will create a csv file with some summary statistics for a list
# including the Shapiro value. It also creates a list of these summaries for
# viewing in R.
DataStatSumTestListFun <- function(DataList,Filename){
  OutputList <- list()
  for (p in 1:length(DataList)) {
    df <- na.omit(DataList[[p]])
    dfStatistics <- data.frame()
    for (g in 1:length(unique(df$Genotype))) {
      Testdf <- dplyr::filter(df,Genotype == unique(df$Genotype)[g])
      MeanValue <- mean(Testdf$RelativeBandIntensity)
      MedianValue <- median(Testdf$RelativeBandIntensity)
      VarValue <- var(Testdf$RelativeBandIntensity)
      StandardDeviation <- sd(Testdf$RelativeBandIntensity)
      StandardError <- StandardDeviation/sqrt(length(Testdf$RelativeBandIntensity))
      result <- shapiro.test(Testdf$RelativeBandIntensity)
      SWpvalue <- result[["p.value"]]
      OneSampleT <- t.test(Testdf$RelativeBandIntensity, mu = 1)
      OneSamplep <- OneSampleT[["p.value"]]
      StatsResults <- c(unique(df$Genotype)[g],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,OneSamplep)
      dfStatistics <- rbind(dfStatistics,StatsResults)
      rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
    }
    names(dfStatistics) <- c("Genotype","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value","One sample T-test p value")
    OutputList[[p]] <- dfStatistics
    names(OutputList)[p] <- names(DataList[p])
    rm(df,dfStatistics)
  }
  SaveResult <- as.data.frame(OutputList)
  SaveResult <- t(SaveResult)
  write.csv(SaveResult, file = paste(Filename,".csv"))
  return(OutputList)
}

SumRelIntNDECD <- DataStatSumTestListFun(RelIntNDECDSepDates,"Summary of NDECD RelInt")
SumRelIntNDECDBg <- DataStatSumTestListFun(RelIntNDECDBgSepDates,"Summary of NDECD RelInt Bg")
SumRelIntNDECDBg0 <- DataStatSumTestListFun(RelIntNDECDBg0SepDates,"Summary of NDECD RelInt Bg0")
SumRelIntNDECDBgSub <- DataStatSumTestListFun(RelIntNDECDBgSubSepDates,"Summary of NDECD RelInt BgSub")
SumRelIntNDECDMinSixSub <- DataStatSumTestListFun(RelIntNDECDMinSixSubSepDates,"Summary of NDECD RelInt MinSixSub")
SumRelIntChangeNDECD <- DataStatSumTestListFun(RelIntChangeNDECDSepDates,"Summary of NDECD RelIntChange")
SumRelIntChangeNDECDBg <- DataStatSumTestListFun(RelIntChangeNDECDBgSepDates,"Summary of NDECD RelInt BgChange")

## Finally we want to generate a set of p-values for the data which can be checked for
## significance and compared to the summaries generated above to know which is most
## relevant. For normal data with equal variance it is the t test, unequal variance use
## the Welch's t test. If the data is nonparametric (a.k.a. not normal), then it is
## the Mann-Whitney test.

Stats_tests <- function(DataList, Filename){
  dfStatistics <- data.frame()
  for (p in 1:length(DataList)) {
    df <- na.omit(DataList[[p]])
    Genotype <- names(DataList)[p]
    for (g in 1:length(unique(df$Genotype))) {
      df1 <- dplyr::filter(df,Genotype == unique(df$Genotype)[1])
      df2 <- dplyr::filter(df,Genotype == unique(df$Genotype)[2])
    }
    StudentT <- t.test(df1$RelativeBandIntensity, df2$RelativeBandIntensity, var.equal = TRUE)
    StudentT <- StudentT[["p.value"]]
    WelchT <- t.test(df1$RelativeBandIntensity, df2$RelativeBandIntensity, var.equal = FALSE)
    WelchT <- WelchT[["p.value"]]
    MannWhitney <- wilcox.test(df1$RelativeBandIntensity,df2$RelativeBandIntensity)
    MannWhitney <- MannWhitney[["p.value"]]
    SummaryVec <- c(Genotype, StudentT, WelchT, MannWhitney)
    dfStatistics <- rbind(dfStatistics, SummaryVec)
    rm(SummaryVec,df1,df2,df,Genotype, StudentT, WelchT, MannWhitney)
  }
  names(dfStatistics) <- c("Genotype","StudentT","WelchT","MannWhitney")
  write.csv(dfStatistics, file = paste(Filename,".csv"))
  return(dfStatistics)
}

SigStatsRelIntNDECD <- Stats_tests(RelIntNDECDSepDates,"SigStats of NDECD RelInt")
SigStatsRelIntNDECDBg <- Stats_tests(RelIntNDECDBgSepDates,"SigStats of NDECD RelInt Bg")
SigStatsRelIntNDECDBg0 <- Stats_tests(RelIntNDECDBg0SepDates,"SigStats of NDECD RelInt Bg0")
SigStatsRelIntNDECDBgSub <- Stats_tests(RelIntNDECDBgSubSepDates,"SigStats of NDECD RelInt BgSub")
SigStatsRelIntNDECDMinSixSub <- Stats_tests(RelIntNDECDMinSixSubSepDates,"SigStats of NDECD RelInt MinSixSub")
SigStatsRelIntChangeNDECD <- Stats_tests(RelIntChangeNDECDSepDates,"SigStats of NDECD RelIntChange")
SigStatsRelIntChangeNDECDBg <- Stats_tests(RelIntChangeNDECDBgSepDates,"SigStats of NDECD RelIntChange Bg")

setwd(mainDir)
#### end ####

#### Saved Plots ####

## Create and save band plots of the bands compared to their own wt data 
for (p in 1:length(MeanGenotypesByDateMinSix)) {
  Labs <- c(title = " ", x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
  df <- MeanGenotypesByDateMinSix[[p]]
  band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank())
  ggsave(paste("MinSix Band plot of ",names(MeanGenotypesByDateMinSix)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(MeanGenotypesByDateRawData)) {
  Labs <- c(title = " ", x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
  df <- MeanGenotypesByDateRawData[[p]]
  band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
    scale_y_continuous(breaks = seq(800,1300, by = 100), limits = c(800,1300))
  ggsave(paste("RawData Band plot of ",names(MeanGenotypesByDateRawData)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

## Create and save RelInt plots of the truncations compare to their own wt data
for (p in 1:length(RelIntNDECDSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Fold change of NICD::mScarlet intensity")
  df <- na.omit(RelIntNDECDSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntNDECDSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 1.7) +
    scale_y_continuous(breaks = seq(1, 1.8, by = 0.1), limits = c(1,1.8))
  ggsave(paste("RelInt plot of ", names(RelIntNDECDSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntChangeNDECDSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Relative change of NICD::mScarlet intensity")
  df <- na.omit(RelIntChangeNDECDSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntChangeNDECDSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 0.7) +
    scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0,0.8))
  ggsave(paste("RelIntChange plot of ", names(RelIntChangeNDECDSepDates)[p], ".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntChangeNDECDBgSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Relative change of NICD::mScarlet intensity")
  df <- na.omit(RelIntChangeNDECDBgSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntChangeNDECDBgSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("RelIntChangeBg plot of ", names(RelIntChangeNDECDBgSepDates)[p], ".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntNDECDBgSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Fold change of NICD::mScarlet intensity")
  df <- na.omit(RelIntNDECDBgSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntNDECDBgSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("Fold change plot of ", names(RelIntNDECDBgSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}
for (p in 1:length(RelIntNDECDBg0SepDates)) {
  Labs <- c(title = " ", x = " ", y = "Fold change of NICD::mScarlet intensity")
  df <- na.omit(RelIntNDECDBg0SepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntNDECDBg0SepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("0 plot of ", names(RelIntNDECDBg0SepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}
for (p in 1:length(RelIntNDECDBgSubSepDates)) {
  Labs <- c(title = " ", x = " ", y = "NICD::mScarlet Fluorescence (AU)")
  df <- na.omit(RelIntNDECDBgSubSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntNDECDBgSubSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("Band minus Bg plot of ", names(RelIntNDECDBgSubSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}
for (p in 1:length(RelIntNDECDMinSixSubSepDates)) {
  Labs <- c(title = " ", x = " ", y = "NICD::mScarlet Fluorescence (AU)")
  df <- na.omit(RelIntNDECDMinSixSubSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntNDECDMinSixSubSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 470) +
    scale_y_continuous(breaks = seq(0, 500, by = 100), limits = c(0,500))
  ggsave(paste("Band minus MinSix plot of ", names(RelIntNDECDMinSixSubSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}
## Plot individual bands and save plots
for (p in 1:length(names(NDECDRawData))) {
  Labs <- c(title = paste("All ",names(NDECDRawData)[p]," Raw Data"), x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
  df <- NDECDRawData[[p]]
  df$Distance <- df$`Distance (um)`
  df$`Distance (um)` <- NULL
  PlotNumber <- length(colnames(df))-1
  df <- df %>% tidyr::gather("id", "value", 1:PlotNumber)
  ggplot(df, aes(Distance,value)) + theme_jmt() +
    geom_point()+
    facet_wrap(~id,ncol = 5) +
    labs(title = Labs[[1]], x = Labs[[2]], y = Labs[[3]]) +
    scale_color_brewer(palette = "Set2")
  ggsave(paste("All ",names(NDECDRawData)[p]," Raw Data",".pdf"), device = "pdf", width = 45, height = 5*(ceiling(PlotNumber/5)), units = "cm", dpi = "retina")
  rm(df,PlotNumber)
}

#### end ####

#### Random Plots ####

## Plots of the raw data
Labs <- c(title = "Raw Data", x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
df <- filter(MeanNDECDRawData,Name %in% c("wtNDECD","NDECDTADPEST","NDECDPEST","NDECDTAD","NDECDOPA"))
band_plot(df, Labs)

## Plots of the MinSix data
Labs <- c(title = "Subtract minimum six values", x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
df <- filter(MeanNDECDMinSixSub,Name %in% c("wtNDECD","NDECDTADPEST","NDECDPEST","NDECDTAD","NDECDOPA"))
band_plot(df, Labs)

## Plots of the BgNorm data
Labs <- c(title = "Normalised to background value", x = "Distance (um)", y = "Normalised NDECD Fluorescence")
df <- filter(MeanNDECDBgNorm,Name %in% c("wtNDECD","NDECDTADPEST","NDECDPEST","NDECDTAD","NDECDOPA"))
band_plot(df, Labs)

## RelInt plots
Labs <- c(title = "Relative intensity of band", x = "Genotype", y = "Relative NICD intensity")
df <- filter(RelIntNDECD,Name %in% c("wtNDECD","NDECDTADPEST","NDECDPEST","NDECDTAD","NDECDOPA"))
RelInt_crossbar(df,Labs)

Labs <- c(title = "Fold change of band over background", x = "Genotype", y = "NICD fold change")
df <- filter(RelIntNDECDBg,Name %in% c("wtNDECD","NDECDTADPEST","NDECDPEST","NDECDTAD","NDECDOPA"))
RelInt_crossbar(df,Labs)

#### end ####

