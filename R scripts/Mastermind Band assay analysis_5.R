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
RelIntSepDates <- function(RelIntDataFrame, GenotypeToSep) {
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
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Mam"
setwd(mainDir)
# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Mam/Background Mam levels.csv", head = TRUE, sep = ",")
BackgroundValues <- data.frame(Nucleus=DF[,1], MamBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Mam/Selected Band assay csv files"

#### Read in and normalise the data ####
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
MamRawData <- list()
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
  names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[1]]), end=-9)
  for (i in 2:length(FilePaths)){
    Nextdf <- read.csv(FilePaths[[i]])
    df <- merge(df,Nextdf, by="X", all = TRUE)
    names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[i]]), end=-9)
    rm(Nextdf)
  }
  names(df)[names(df) == 'X'] <- "Distance (um)"
  MamRawData[[f]] <- df
  rm(df)
  rm(FilePaths)
  names(MamRawData)[f] <- FolderNames[f] 
}

rm(FolderNames)

MamminsixSub <- list()
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] - NormValue
    rm(NormValue)
  }
  MamminsixSub[[s]] <- df
  names(MamminsixSub)[s] <- names(MamRawData[s])
  rm(df,ColNames)
}

MamBgNorm <- list()
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    NormValue <- BackgroundValues$MamBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  MamBgNorm[[s]] <- df
  names(MamBgNorm)[s] <- names(MamRawData[s])
  rm(df,ColNames)
}

MeanMamRawData <- MeanDF(MamRawData)
MeanMamBgNorm <- MeanDF(MamBgNorm)
MeanMamMinSixSub <- MeanDF(MamminsixSub)

#### end ####

#### Select only data from the same days for comparing genotypes ####
WTMamRDSepDates <- SepDates(MamRawData,"Wild Type")
WTMamminsixSubSepDates <- SepDates(MamminsixSub,"Wild Type")

MeanGenotypesByDateRawData <- WTvsMutByDate(WTMamRDSepDates, MamRawData)
MeanGenotypesByDateMinSix <- WTvsMutByDate(WTMamminsixSubSepDates, MamminsixSub)

#### end ####

#### Calculate fold increase of band over background ####
## Here we will use the mean of the 5 middle values and divide by the mean of the 5 values at
## each edge of the band (10 values in total).
RelIntMam <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
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
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntMam <- rbind(RelIntMam,RelInt)
  rm(RelInt, Band, Background, df)
}

RelIntMamSepDates <- (RelIntSepDates(RelIntMam,"Wild Type"))

RelIntMamEnds <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntMamEnds <- rbind(RelIntMamEnds,RelInt)
  rm(RelInt, Band, Background, df)
}

RelIntMamEndsSepDates <- (RelIntSepDates(RelIntMamEnds,"Wild Type"))

RelIntMamBg <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$MamBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background, BandSubBg)
    RelIntMamBg <- rbind(RelIntMamBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntMamBgSepDates <- RelIntSepDates(RelIntMamBg, "Wild Type")

RelIntMamMinSixSub <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
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
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
  RelIntMamMinSixSub <- rbind(RelIntMamMinSixSub,RelInt)
  rm(RelInt, Band, Background, df,dfOrdered,c,Column,MinSix)
}
RelIntMamMinSixSubSepDates <- RelIntSepDates(RelIntMamMinSixSub, "Wild Type")

RelIntMamBg0 <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$MamBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    BandSubBg <- Band - Background
    RelInt <- as.data.frame(BandSubBg/Band)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `BandSubBg/Band`)
    rm(Band, Background, BandSubBg)
    RelIntMamBg0 <- rbind(RelIntMamBg0,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntMamBg0SepDates <- RelIntSepDates(RelIntMamBg0, "Wild Type")

RelIntMamBgSub <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$MamBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band-Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
    rm(Band, Background, BandSubBg)
    RelIntMamBgSub <- rbind(RelIntMamBgSub,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntMamBgSubSepDates <- RelIntSepDates(RelIntMamBgSub, "Wild Type")

RelIntChangeMamBg <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$MamBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeMamBg <- rbind(RelIntChangeMamBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntChangeMamBgSepDates <- RelIntSepDates(RelIntChangeMamBg, "Wild Type")

RelIntChangeMam <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
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
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeMam <- rbind(RelIntChangeMam,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
RelIntChangeMamSepDates <- (RelIntSepDates(RelIntChangeMam, "Wild Type"))
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

NormDis_plots(RelIntMamSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelInt plots")
NormDis_plots(RelIntMamBgSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntBg plots")
NormDis_plots(RelIntMamBg0SepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntBg0 plots")
NormDis_plots(RelIntMamBgSubSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntBgSub plots")
NormDis_plots(RelIntMamMinSixSubSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntMinSixSub plots")
NormDis_plots(RelIntChangeMamSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntChange plots")
NormDis_plots(RelIntChangeMamBgSepDates, "RelativeBandIntensity", "Genotype","Normal distribution of RelIntChangeBg plots")

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

QQ_plots(RelIntMamSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelInt")
QQ_plots(RelIntMamBgSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBg")
QQ_plots(RelIntMamBg0SepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBg0")
QQ_plots(RelIntMamBgSubSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBgSub")
QQ_plots(RelIntMamMinSixSubSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntBgSub")
QQ_plots(RelIntChangeMamSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntChange")
QQ_plots(RelIntChangeMamBgSepDates, "RelativeBandIntensity", "Genotype","QQ Plots of RelIntChangeBg")

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

SumRelIntMam <- DataStatSumTestListFun(RelIntMamSepDates,"Summary of Mam RelInt")
SumRelIntMamBg <- DataStatSumTestListFun(RelIntMamBgSepDates,"Summary of Mam RelInt Bg")
SumRelIntMamBg0 <- DataStatSumTestListFun(RelIntMamBg0SepDates,"Summary of Mam RelInt Bg0")
SumRelIntMamBgSub <- DataStatSumTestListFun(RelIntMamBgSubSepDates,"Summary of Mam RelInt BgSub")
SumRelIntMamMinSixSub <- DataStatSumTestListFun(RelIntMamMinSixSubSepDates,"Summary of Mam RelInt MinSixSub")
SumRelIntChangeMam <- DataStatSumTestListFun(RelIntChangeMamSepDates,"Summary of Mam RelIntChange")
SumRelIntChangeMamBg <- DataStatSumTestListFun(RelIntChangeMamBgSepDates,"Summary of Mam RelInt BgChange")

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

SigStatsRelIntMam <- Stats_tests(RelIntMamSepDates,"SigStats of Mam RelInt")
SigStatsRelIntMamBg <- Stats_tests(RelIntMamBgSepDates,"SigStats of Mam RelInt Bg")
SigStatsRelIntMamBg0 <- Stats_tests(RelIntMamBg0SepDates,"SigStats of Mam RelInt Bg0")
SigStatsRelIntMamBgSub <- Stats_tests(RelIntMamBgSubSepDates,"SigStats of Mam RelInt BgSub")
SigStatsRelIntMamMinSixSub <- Stats_tests(RelIntMamMinSixSubSepDates,"SigStats of Mam RelInt MinSixSub")
SigStatsRelIntChangeMam <- Stats_tests(RelIntChangeMamSepDates,"SigStats of Mam RelIntChange")
SigStatsRelIntChangeMamBg <- Stats_tests(RelIntChangeMamBgSepDates,"SigStats of Mam RelIntChange Bg")

setwd(mainDir)
#### end ####

#### Saved Plots ####

## Create and save band plots of the bands compared to their own wt data 
for (p in 1:length(MeanGenotypesByDateMinSix)) {
  Labs <- c(title = " ", x = "Distance (um)", y = "Mam::sfGFP Fluorescence (AU)")
  df <- MeanGenotypesByDateMinSix[[p]]
  band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank())
  ggsave(paste("MinSix Band plot of ",names(MeanGenotypesByDateMinSix)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(MeanGenotypesByDateRawData)) {
  Labs <- c(title = " ", x = "Distance (um)", y = "Mam::sfGFP Fluorescence (AU)")
  df <- MeanGenotypesByDateRawData[[p]]
  band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
    scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0,2000))
  ggsave(paste("RawData Band plot of ",names(MeanGenotypesByDateRawData)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

## Create and save RelInt plots of the truncations compare to their own wt data
for (p in 1:length(RelIntMamSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Fold change of Mam::sfGFP intensity")
  df <- na.omit(RelIntMamSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntMamSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("RelInt plot of ", names(RelIntMamSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntChangeMamSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Relative change of Mam::sfGFP intensity")
  df <- na.omit(RelIntChangeMamSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntChangeMamSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("RelIntChange plot of ", names(RelIntChangeMamSepDates)[p], ".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntChangeMamBgSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Relative change of Mam::sfGFP intensity")
  df <- na.omit(RelIntChangeMamBgSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntChangeMamBgSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 9.5) +
    scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0,10))
  ggsave(paste("RelIntChangeBg plot of ", names(RelIntChangeMamBgSepDates)[p], ".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntMamBgSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Fold change of Mam::sfGFP intensity")
  df <- na.omit(RelIntMamBgSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntMamBgSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 11.5) +
    scale_y_continuous(breaks = seq(0, 12, by = 3), limits = c(0,12))
  ggsave(paste("Fold change plot of ", names(RelIntMamBgSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntMamBg0SepDates)) {
  Labs <- c(title = " ", x = " ", y = "Fold change of Mam::sfGFP intensity")
  df <- na.omit(RelIntMamBg0SepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntMamBg0SepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black")
  ggsave(paste("0 plot of ", names(RelIntMamBg0SepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}

for (p in 1:length(RelIntMamBgSubSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Mastermind::sfGFP fluorescence (AU)")
  df <- na.omit(RelIntMamBgSubSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntMamBgSubSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 3250) +
    scale_y_continuous(breaks = seq(0, 3500, by = 500), limits = c(0,3500))
  ggsave(paste("Band minus Bg plot of ", names(RelIntMamBgSubSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}  
for (p in 1:length(RelIntMamMinSixSubSepDates)) {
  Labs <- c(title = " ", x = " ", y = "Mastermind::sfGFP Fluorescence (AU)")
  df <- na.omit(RelIntMamMinSixSubSepDates[[p]])
  RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
    geom_signif(comparisons = list(c(names(RelIntMamMinSixSubSepDates)[p],"Wild Type")), 
                map_signif_level=TRUE, test = "t.test", color = "black", y_position = 3200) +
    scale_y_continuous(breaks = seq(0, 3500, by = 500), limits = c(0,3500))
  ggsave(paste("Band minus MinSix plot of ", names(RelIntMamMinSixSubSepDates)[p],".jpg"), device = "jpeg", dpi = "retina",
         width = 12, height = 12, units = "cm")
}
## Plot individual bands and save plots
for (p in 1:length(names(MamRawData))) {
  Labs <- c(title = paste("All ",names(MamRawData)[p]," Raw Data"), x = "Distance (um)", y = "Mam::GFP Fluorescence (AU)")
  df <- MamRawData[[p]]
  df$Distance <- df$`Distance (um)`
  df$`Distance (um)` <- NULL
  PlotNumber <- length(colnames(df))-1
  df <- df %>% tidyr::gather("id", "value", 1:PlotNumber)
  ggplot(df, aes(Distance,value)) + theme_jmt() +
    geom_point()+
    facet_wrap(~id,ncol = 5) +
    labs(title = Labs[[1]], x = Labs[[2]], y = Labs[[3]]) +
    scale_color_brewer(palette = "Set2")
  ggsave(paste("All ",names(MamRawData)[p]," Raw Data",".pdf"), device = "pdf", width = 45, height = 5*(ceiling(PlotNumber/5)), units = "cm", dpi = "retina")
  rm(df,PlotNumber)
}

#### end ####

#### Random plots ####

## Plots of the raw data
Labs <- c(title = "Raw Data", x = "Distance (um)", y = "Mam::GFP Fluorescence (AU)")
df <- filter(MeanMamRawData,Name %in% c("Wild Type","ΔTADΔPEST"))
band_plot(df, Labs)

## Plots of the MinSix data
Labs <- c(title = "Subtract minimum six values", x = "Distance (um)", y = "Mam::GFP Fluorescence (AU)")
df <- filter(MeanMamMinSixSub,Name %in% c("wtNDECD","NDECDTADPEST"))
band_plot(df, Labs)

## Plots of the BgNorm data
Labs <- c(title = "Fold enrichment over background", x = "Distance (um)", y = "Normalised Mam Fluorescence")
df <- filter(MeanMamBgNorm,Name %in% c("wtNDECD","NDECDTADPEST","NDECDTAD"))
band_plot(df, Labs)

## Compare the WT control conditions
Labs <- c(title = "Wild Type Raw Data", x = "Distance (um)", y = "Mam::GFP Fluorescence (AU)")
df <- filter(MeanDF(WTMamRDSepDates), Name %in% c("wtNDECD","NDECDTADPEST","NDECDTAD","NDECDPEST"))
band_plot(df, Labs) 

## RelInt plots
Labs <- c(title = "Relative intensity of band", x = "Genotype", y = "Relative Mam intensity")
df <- filter(RelIntMam,Genotype %in% c("wtNDECD","NDECDTADPEST"))
RelInt_crossbar(df,Labs)

#### end ####

