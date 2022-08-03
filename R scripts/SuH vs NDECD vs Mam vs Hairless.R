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

DataStatSumTestListFun <- function(DataList,MuNumber,Filename){
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
      OneSampleT <- t.test(Testdf$RelativeBandIntensity, mu = MuNumber)
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
  ggplot(Data, aes(Distance,Mean,group = Name,col=Molecule,shape=Molecule,ymin=Mean-SEM,ymax=Mean+SEM)) + theme_jmt() +
    geom_ribbon(linetype = "dashed", fill = "grey90",alpha = 0.5, size = 0.25) + geom_line(size = 0.5) + geom_point() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    scale_color_brewer(palette = "Set2")
}

## Define a function to plot the relative intensities of the bands
RelInt_crossbar <- function(Data, Labels){
  ggplot(Data,aes(Molecule,RelativeBandIntensity,col=Molecule,shape=Molecule)) + theme_jmt() +
    geom_jitter(position = position_jitterdodge(jitter.width=0.5)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

RelInt_crossbar2 <- function(Data, Labels){
  ggplot(Data,aes(Molecule,RelativeBandIntensity,col=Molecule,shape=Molecule)) + theme_jmt() +
    geom_jitter(position = position_jitterdodge(jitter.width=1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

Background_crossbar2 <- function(Data, Labels){
  ggplot(Data,aes(Molecule,Background,col=Molecule,shape=Molecule)) + theme_jmt() +
    geom_jitter(position = position_jitterdodge(jitter.width=1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}
#### end ####

mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Comparisons between molecules"
setwd(mainDir)

#### Read in SuH data ####

# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Su(H)/Background SuH levels.csv", head = TRUE, sep = ",")
SuHBackgroundValues <- data.frame(Nucleus=DF[,1], SuHBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Su(H)/Band assay csv files"

# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
SuHRawData <- list()
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
  SuHRawData[[f]] <- df
  rm(df)
  rm(FilePaths)
  names(SuHRawData)[f] <- FolderNames[f] 
}

rm(FolderNames)

SuHminsixSub <- list()
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] - NormValue
    rm(NormValue)
  }
  SuHminsixSub[[s]] <- df
  names(SuHminsixSub)[s] <- names(SuHRawData[s])
  rm(df,ColNames)
}

MeanSuHRawData <- MeanDF(SuHRawData)
MeanSuHMinSixSub <- MeanDF(SuHminsixSub)
#### end ####

#### Read in Notch Data ####

# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/NDECD/Background NDECD values.csv", head = TRUE, sep = ",")
NDECDBackgroundValues <- data.frame(Nucleus=DF[,1], NDECDBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/NDECD/Band assay csv files"

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

MeanNDECDRawData <- MeanDF(NDECDRawData)
MeanNDECDMinSixSub <- MeanDF(NDECDminsixSub)
#### end ####

#### Read in Mam Data ####
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Mam/Background Mam levels.csv", head = TRUE, sep = ",")
MamBackgroundValues <- data.frame(Nucleus=DF[,1], MamBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Mam/Band assay csv files"

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

MeanMamRawData <- MeanDF(MamRawData)
MeanMamMinSixSub <- MeanDF(MamminsixSub)
#### end ####

#### Read in Hairless Data ####
# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Hairless/Background Hairless levels.csv", head = TRUE, sep = ",")
HairlessBackgroundValues <- data.frame(Nucleus=DF[,1], HairlessBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/Truncations/Hairless/Band assay csv files"

# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
HairlessRawData <- list()
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
  names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[1]]), end=-14)
  for (i in 2:length(FilePaths)){
    Nextdf <- read.csv(FilePaths[[i]])
    df <- merge(df,Nextdf, by="X", all = TRUE)
    names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[i]]), end=-14)
    rm(Nextdf)
  }
  names(df)[names(df) == 'X'] <- "Distance (um)"
  HairlessRawData[[f]] <- df
  rm(df)
  rm(FilePaths)
  names(HairlessRawData)[f] <- FolderNames[f] 
}

rm(FolderNames)

HairlessminsixSub <- list()
for (s in 1:length(HairlessRawData)) {
  df <- HairlessRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] - NormValue
    rm(NormValue)
  }
  HairlessminsixSub[[s]] <- df
  names(HairlessminsixSub)[s] <- names(HairlessRawData[s])
  rm(df,ColNames)
}

MeanHairlessRawData <- MeanDF(HairlessRawData)
MeanHairlessMinSixSub <- MeanDF(HairlessminsixSub)

#### end ####

#### Calculate fold increase of band over background ####
RelIntSuHBg <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- SuHBackgroundValues$SuHBG[which(grepl(ColNames[[t]], SuHBackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(SuHRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background, BandSubBg)
    RelIntSuHBg <- rbind(RelIntSuHBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntChangeSuHBg <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- SuHBackgroundValues$SuHBG[which(grepl(ColNames[[t]], SuHBackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(SuHRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeSuHBg <- rbind(RelIntChangeSuHBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntChangeSuHBgSepDates <- RelIntSepDates(RelIntChangeSuHBg, "Wild Type")
DataStatSumTestListFun(RelIntChangeSuHBgSepDates,0,"Summary of SuH RelIntChange Bg")
RelIntSuHBgSepDates <- RelIntSepDates(RelIntSuHBg, "Wild Type")
DataStatSumTestListFun(RelIntSuHBgSepDates,1,"Summary of SuH RelInt Bg")

RelIntNDECDBg <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- NDECDBackgroundValues$NDECDBG[which(grepl(ColNames[[t]], NDECDBackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(NDECDRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background, BandSubBg)
    RelIntNDECDBg <- rbind(RelIntNDECDBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntChangeNDECDBg <- NULL
for (s in 1:length(NDECDRawData)) {
  df <- NDECDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- NDECDBackgroundValues$NDECDBG[which(grepl(ColNames[[t]], NDECDBackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(NDECDRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeNDECDBg <- rbind(RelIntChangeNDECDBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

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
DataStatSumTestListFun(RelIntChangeNDECDSepDates,0,"Summary of NDECD RelIntChange")
RelIntNDECDMinSixSubSepDates <- RelIntSepDates(RelIntNDECDMinSixSub, "Wild Type")
DataStatSumTestListFun(RelIntNDECDMinSixSubSepDates,0,"Summary of NDECD RelIntMinSixSub")
RelIntChangeNDECDBgSepDates <- RelIntSepDates(RelIntChangeNDECDBg, "Wild Type")
DataStatSumTestListFun(RelIntChangeNDECDBgSepDates,0,"Summary of NDECD RelIntChange Bg")
RelIntNDECDBgSepDates <- RelIntSepDates(RelIntNDECDBg, "Wild Type")
DataStatSumTestListFun(RelIntNDECDBgSepDates,1,"Summary of NDECD RelInt Bg")

RelIntMamBg <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- MamBackgroundValues$MamBG[which(grepl(ColNames[[t]], MamBackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background, BandSubBg)
    RelIntMamBg <- rbind(RelIntMamBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntChangeMamBg <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- MamBackgroundValues$MamBG[which(grepl(ColNames[[t]], MamBackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeMamBg <- rbind(RelIntChangeMamBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

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
DataStatSumTestListFun(RelIntChangeMamSepDates,0,"Summary of Mam RelIntChange")
RelIntMamMinSixSubSepDates <- RelIntSepDates(RelIntMamMinSixSub, "Wild Type")
DataStatSumTestListFun(RelIntMamMinSixSubSepDates,0,"Summary of Mam RelIntMinSixSub")
RelIntChangeMamBgSepDates <- RelIntSepDates(RelIntChangeMamBg, "Wild Type")
DataStatSumTestListFun(RelIntChangeMamBgSepDates,0,"Summary of Mam RelIntChange Bg")
RelIntMamBgSepDates <- RelIntSepDates(RelIntMamBg, "Wild Type")
DataStatSumTestListFun(RelIntMamBgSepDates,1,"Summary of Mam RelInt Bg")

RelIntHairlessBg <- NULL
for (s in 1:length(HairlessRawData)) {
  df <- HairlessRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- HairlessBackgroundValues$HairlessBG[which(grepl(ColNames[[t]], HairlessBackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(HairlessRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background, BandSubBg)
    RelIntHairlessBg <- rbind(RelIntHairlessBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntChangeHairlessBg <- NULL
for (s in 1:length(HairlessRawData)) {
  df <- HairlessRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- HairlessBackgroundValues$HairlessBG[which(grepl(ColNames[[t]], HairlessBackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(HairlessRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeHairlessBg <- rbind(RelIntChangeHairlessBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntChangeHairless <- NULL
for (s in 1:length(HairlessRawData)) {
  df <- HairlessRawData[[s]]
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
  RelInt$Genotype <- names(HairlessRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeHairless <- rbind(RelIntChangeHairless,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}

RelIntChangeHairlessSepDates <- (RelIntSepDates(RelIntChangeHairless, "Wild Type"))
DataStatSumTestListFun(RelIntChangeHairlessSepDates,0,"Summary of Hairless RelIntChange")
RelIntChangeHairlessBgSepDates <- RelIntSepDates(RelIntChangeHairlessBg, "Wild Type")
DataStatSumTestListFun(RelIntChangeHairlessBgSepDates,0,"Summary of Hairless RelIntChange Bg")
RelIntHairlessBgSepDates <- RelIntSepDates(RelIntHairlessBg, "Wild Type")
DataStatSumTestListFun(RelIntHairlessBgSepDates,1,"Summary of Hairless RelInt Bg")
#### end ####

#### Combine the relevant data sets for the molecules into useful dataframes ####

## Raw Data
MoleculeWTRawDataList <- list()

SuHRawWT <- filter(MeanSuHRawData, Name == "Wild Type")
Nnumber <- ncol(SuHRawData[["Wild Type"]])-1
SuHRawWT$Name <- "SuH"
SuHRawWT$Molecule <- paste("Su(H)::GFP (N = ", Nnumber, ")", sep = "")
MoleculeWTRawDataList[["SuH"]] <- SuHRawWT

NDECDRawWT <- filter(MeanNDECDRawData, Name == "Wild Type")
Nnumber <- ncol(NDECDRawData[["Wild Type"]])-1
NDECDRawWT$Name <- "NDECD"
NDECDRawWT$Molecule <- paste("NΔECD::mScarlet (N = ", Nnumber, ")", sep = "")
MoleculeWTRawDataList[["NDECD"]] <- NDECDRawWT

MamRawWT <- filter(MeanMamRawData, Name == "Wild Type")
Nnumber <- ncol(MamRawData[["Wild Type"]])-1
MamRawWT$Name <- "Mam"
MamRawWT$Molecule <- paste("Mam::GFP (N = ", Nnumber, ")", sep = "")
MoleculeWTRawDataList[["Mam"]] <- MamRawWT

HairlessRawWT <- filter(MeanHairlessRawData, Name == "Wild Type")
Nnumber <- ncol(HairlessRawData[["Wild Type"]])-1
HairlessRawWT$Name <- "Hairless"
HairlessRawWT$Molecule <- paste("Hairless::GFP (N = ", Nnumber, ")", sep = "")
MoleculeWTRawDataList[["Hairless"]] <- HairlessRawWT

rm(SuHRawWT,NDECDRawWT,MamRawWT,HairlessRawWT)
CombinedSuHNDECDRawWT <- rbind(MoleculeWTRawDataList[["SuH"]],MoleculeWTRawDataList[["NDECD"]])
CombinedSuHMamHairlessRawWT <- rbind(MoleculeWTRawDataList[["SuH"]],MoleculeWTRawDataList[["Mam"]],MoleculeWTRawDataList[["Hairless"]])
CombinedSuHMamHairlessRawWT <- CombinedSuHMamHairlessRawWT %>% mutate(Molecule = fct_relevel(Molecule, c("Hairless::GFP (N = 32)","Su(H)::GFP (N = 82)","Mam::GFP (N = 49)")))

## MinSixSub
MoleculeWTminsixsubList <- list()

SuHminsixsubWT <- filter(MeanSuHMinSixSub, Name == "Wild Type")
SuHminsixsubWT$Name <- "SuH"
Nnumber <- ncol(SuHminsixSub[["Wild Type"]])-1
SuHminsixsubWT$Molecule <- paste("Su(H)::GFP (N = ", Nnumber, ")", sep = "")
MoleculeWTminsixsubList[["SuH"]] <- SuHminsixsubWT

NDECDminsixsubWT <- filter(MeanNDECDMinSixSub, Name == "Wild Type")
NDECDminsixsubWT$Name <- "NDECD"
Nnumber <- ncol(NDECDminsixSub[["Wild Type"]])-1
NDECDminsixsubWT$Molecule <- paste("NΔECD::mScarlet (N = ", Nnumber, ")", sep = "")
MoleculeWTminsixsubList[["NDECD"]] <- NDECDminsixsubWT

MamminsixsubWT <- filter(MeanMamMinSixSub, Name == "Wild Type")
MamminsixsubWT$Name <- "Mam"
Nnumber <- ncol(MamminsixSub[["Wild Type"]])-1
MamminsixsubWT$Molecule <- paste("Mam::GFP (N = ", Nnumber, ")", sep = "")
MoleculeWTminsixsubList[["Mam"]] <- MamminsixsubWT

HairlessminsixsubWT <- filter(MeanHairlessMinSixSub, Name == "Wild Type")
HairlessminsixsubWT$Name <- "Hairless"
Nnumber <- ncol(HairlessminsixSub[["Wild Type"]])-1
HairlessminsixsubWT$Molecule <- paste("Hairless::GFP (N = ", Nnumber, ")", sep = "")
MoleculeWTminsixsubList[["Hairless"]] <- HairlessminsixsubWT

rm(SuHminsixsubWT,NDECDminsixsubWT,MamminsixsubWT,HairlessminsixsubWT)
CombinedSuHNDECDminsixsubWT <- rbind(MoleculeWTminsixsubList[["SuH"]],MoleculeWTminsixsubList[["NDECD"]])
CombinedSuHMamHairlessminsixsubWT <- rbind(MoleculeWTminsixsubList[["SuH"]],MoleculeWTminsixsubList[["Mam"]],MoleculeWTminsixsubList[["Hairless"]])
CombinedSuHMamHairlessminsixsubWT <- CombinedSuHMamHairlessminsixsubWT %>% mutate(Molecule = fct_relevel(Molecule, c("Hairless::GFP (N = 32)","Su(H)::GFP (N = 82)","Mam::GFP (N = 49)")))

##RelIntBg
SuHRelInt <- filter(RelIntSuHBg, Genotype == "Wild Type")
SuHRelInt$Molecule <- "Su(H)::GFP"
NDECDRelInt <- filter(RelIntNDECDBg, Genotype == "Wild Type")
NDECDRelInt$Molecule <- "NDECD::mScarlet"
MamRelInt <- filter(RelIntMamBg, Genotype == "Wild Type")
MamRelInt$Molecule <- "Mam::GFP"
HairlessRelInt <- filter(RelIntHairlessBg, Genotype == "Wild Type")
HairlessRelInt$Molecule <- "Hairless::GFP"

SuHRelIntChange <- filter(RelIntChangeSuHBg, Genotype == "Wild Type")
SuHRelIntChange$Molecule <- "Su(H)::GFP"
NDECDRelIntChange <- filter(RelIntChangeNDECDBg, Genotype == "Wild Type")
NDECDRelIntChange$Molecule <- "NDECD::mScarlet"
MamRelIntChange <- filter(RelIntChangeMamBg, Genotype == "Wild Type")
MamRelIntChange$Molecule <- "Mam::GFP"
HairlessRelIntChange <- filter(RelIntChangeHairlessBg, Genotype == "Wild Type")
HairlessRelIntChange$Molecule <- "Hairless::GFP"

CombinedRelIntList <- list()
CombinedRelIntList[["NDECD"]] <- rbind(SuHRelInt,NDECDRelInt)
CombinedRelIntList[["Mam"]] <- rbind(SuHRelInt,MamRelInt)
CombinedRelIntList[["Hairless"]] <- rbind(SuHRelInt,HairlessRelInt)

CombinedRelIntChangeList <- list()
CombinedRelIntChangeList[["NDECD"]] <- rbind(SuHRelIntChange,NDECDRelIntChange)
CombinedRelIntChangeList[["Mam"]] <- rbind(SuHRelIntChange,MamRelIntChange)
CombinedRelIntChangeList[["Hairless"]] <- rbind(SuHRelIntChange,HairlessRelIntChange)

WTbyMolecule <- list()
WTbyMolecule[["SuH"]] <- SuHRelInt
WTbyMolecule[["NDECD"]] <- NDECDRelInt
WTbyMolecule[["Mam"]] <- MamRelInt
WTbyMolecule[["Hairless"]] <- HairlessRelInt
DataStatSumTestListFun(WTbyMolecule, 1, "Summary of different molecules RelInt Bg")

WTbyMolecule <- list()
WTbyMolecule[["SuH"]] <- SuHRelIntChange
WTbyMolecule[["NDECD"]] <- NDECDRelIntChange
WTbyMolecule[["Mam"]] <- MamRelIntChange
WTbyMolecule[["Hairless"]] <- HairlessRelIntChange
DataStatSumTestListFun(WTbyMolecule, 0, "Summary of different molecules RelIntChange")

## Combine SuH Mam and Hairless in one dataframe with their background values
SuHRelInt <- merge(SuHRelInt,SuHBackgroundValues,by="Nucleus")
names(SuHRelInt)[names(SuHRelInt) == "SuHBG"] <- "Background"
NDECDRelInt <- merge(NDECDRelInt,NDECDBackgroundValues,by="Nucleus")
names(NDECDRelInt)[names(NDECDRelInt) == "NDECDBG"] <- "Background"
MamRelInt <- merge(MamRelInt,MamBackgroundValues,by="Nucleus")
names(MamRelInt)[names(MamRelInt) == "MamBG"] <- "Background"
HairlessRelInt <- merge(HairlessRelInt,HairlessBackgroundValues,by="Nucleus")
names(HairlessRelInt)[names(HairlessRelInt) == "HairlessBG"] <- "Background"

SuHMamHairless <- rbind(SuHRelInt,MamRelInt,HairlessRelInt)
SuHMamHairless <- SuHMamHairless %>% mutate(Molecule = fct_relevel(Molecule, c("Hairless::GFP","Su(H)::GFP","Mam::GFP")))

SuHRelIntChange <- merge(SuHRelIntChange,SuHBackgroundValues,by="Nucleus")
names(SuHRelIntChange)[names(SuHRelIntChange) == "SuHBG"] <- "Background"
NDECDRelIntChange <- merge(NDECDRelIntChange,NDECDBackgroundValues,by="Nucleus")
names(NDECDRelIntChange)[names(NDECDRelIntChange) == "NDECDBG"] <- "Background"
MamRelIntChange <- merge(MamRelIntChange,MamBackgroundValues,by="Nucleus")
names(MamRelIntChange)[names(MamRelIntChange) == "MamBG"] <- "Background"
HairlessRelIntChange <- merge(HairlessRelIntChange,HairlessBackgroundValues,by="Nucleus")
names(HairlessRelIntChange)[names(HairlessRelIntChange) == "HairlessBG"] <- "Background"

SuHMamHairlessRelIntChange <- rbind(SuHRelIntChange,MamRelIntChange,HairlessRelIntChange)
SuHMamHairlessRelIntChange <- SuHMamHairlessRelIntChange %>% mutate(Molecule = fct_relevel(Molecule, c("Hairless::GFP","Su(H)::GFP","Mam::GFP")))


rm(SuHRelInt,MamRelInt,HairlessRelInt,NDECDRelInt)
rm(SuHRelIntChange,MamRelIntChange,HairlessRelIntChange,NDECDRelIntChange)
#### end ####

#### Individual molecule saved plots ####
subDir <- "Individual molecules"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

MoleculeList <- list()

df <- filter(RelIntChangeSuHBg, Genotype == "Wild Type")
df$Molecule <- "Su(H)"
MoleculeList[[1]] <- df
names(MoleculeList)[1] <- "SuH" 
Labs <- c(title = " ", x = " ", y = "Relative change of Su(H)::EGFP intensity")
RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0,5))
ggsave(paste("RelIntChangeBg plot of SuH.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

df <- filter(RelIntChangeHairless, Genotype == "Wild Type")
df$Molecule <- "Hairless"
MoleculeList[[2]] <- df
names(MoleculeList)[2] <- "Hairless" 
Labs <- c(title = " ", x = " ", y = "Relative change of Hairless::EGFP intensity")
RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0,1.2))
ggsave(paste("RelIntChange plot of Hairless.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

df <- filter(RelIntChangeMamBg, Genotype == "Wild Type")
df$Molecule <- "Mastermind"
MoleculeList[[3]] <- df
names(MoleculeList)[3] <- "Mastermind" 
Labs <- c(title = " ", x = " ", y = "Relative change of Mam::sfGFP intensity")
RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0,10))
ggsave(paste("RelIntChangeBg plot of Mastermind.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

df <- filter(RelIntMamMinSixSub, Genotype == "Wild Type")
df$Molecule <- "Mastermind"
MoleculeList[[3]] <- df
names(MoleculeList)[3] <- "Mastermind" 
Labs <- c(title = " ", x = " ", y = "Mam::sfGFP fluorescence (AU)")
RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 3500, by = 500), limits = c(0,3500))
ggsave(paste("RelIntMinSixSub plot of Mastermind.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

df <- filter(RelIntChangeNDECD, Genotype == "Wild Type")
df$Molecule <- "NICD"
MoleculeList[[4]] <- df
names(MoleculeList)[4] <- "NICD" 
Labs <- c(title = " ", x = " ", y = "Relative change of NICD::mScarlet intensity")
RelInt_crossbar(df, Labs) + theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 0.7, by = 0.1), limits = c(0,0.7))
ggsave(paste("RelIntChange plot of NICD.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

DataStatSumTestListFun(MoleculeList,0,"Summary of Individual Molecule graphs")

rm(df,MoleculeList)

Labs <- c(title = " ", x = "Distance (μm)", y = "Su(H)::EGFP Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["SuH"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of raw SuH.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "NICD::mScarlet Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["NDECD"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of raw NICD.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Hairless::EGFP Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["Hairless"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of raw Hairless.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Mam::sfGFP Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["Mam"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of raw Mastermind.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

Labs <- c(title = " ", x = "Distance (μm)", y = "Su(H)::EGFP Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["SuH"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(breaks = seq(0, 2000, by = 400), limits = c(0,2000))
ggsave(paste("Band plot of raw SuH SAME SCALE.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "NICD::mScarlet Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["NDECD"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
scale_y_continuous(breaks = seq(0, 2000, by = 400), limits = c(0,2000))
ggsave(paste("Band plot of raw NICD SAME SCALE.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Hairless::EGFP Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["Hairless"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(breaks = seq(0, 2000, by = 400), limits = c(0,2000))
ggsave(paste("Band plot of raw Hairless SAME SCALE.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Mam::sfGFP Fluorescence (AU)")
band_plot(MoleculeWTRawDataList[["Mam"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(breaks = seq(0, 2000, by = 400), limits = c(0,2000))
ggsave(paste("Band plot of raw Mastermind SAME SCALE.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

setwd(mainDir)
#### end ####

#### Saved Plots ####
## Generate band plot comparing SuH and NDECD
Labs <- c(title = " ", x = "Distance (μm)", y = "Fluorescence (AU)")
band_plot(CombinedSuHNDECDRawWT, Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of raw SuH vs NDECD.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Fluorescence (AU)")
band_plot(CombinedSuHNDECDminsixsubWT, Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of MinSixSub SuH vs NDECD.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Fluorescence (AU)")
band_plot(CombinedSuHMamHairlessRawWT, Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of raw SuH vs Mam vs Hairless.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm") 

Labs <- c(title = " ", x = "Distance (μm)", y = "Fluorescence (AU)")
band_plot(CombinedSuHMamHairlessminsixsubWT, Labs) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste("Band plot of MinSixSub SuH vs Mam vs Hairless.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm") 

## Generate plot for relative intensity to background for SuH vs NDECD
Labs <- c(title = " ", x = " ", y = "Relative intensity")
RelInt_crossbar(CombinedRelIntList[["NDECD"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_signif(comparisons = list(c("Su(H)::GFP","NDECD::mScarlet")), 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("Fold change plot of SuH vs NDECD.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

## Generate plot for relative intensity to background for SuH vs Mam vs Hairless
Labs <- c(title = " ", x = " ", y = "Relative intensity")
RelInt_crossbar2(SuHMamHairless, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_signif(comparisons = list(c("Su(H)::GFP","Hairless::GFP")), y_position = max(filter(SuHMamHairless,Molecule == "Su(H)::GFP")$RelativeBandIntensity)+0.2, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Su(H)::GFP","Mam::GFP")), y_position = max(filter(SuHMamHairless,Molecule == "Mam::GFP")$RelativeBandIntensity)+0.2, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Mam::GFP","Hairless::GFP")), y_position = max(filter(SuHMamHairless,Molecule == "Mam::GFP")$RelativeBandIntensity)+1, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("Fold change plot of SuH vs Mam vs Hairless.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

Labs <- c(title = " ", x = " ", y = "Relative intensity")
RelInt_crossbar(CombinedRelIntChangeList[["NDECD"]], Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_signif(comparisons = list(c("Su(H)::GFP","NDECD::mScarlet")), 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("RelIntBgChange plot of SuH vs NDECD.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

## Generate plot for relative intensity to background for SuH vs Mam vs Hairless
Labs <- c(title = " ", x = " ", y = "Relative intensity")
RelInt_crossbar2(SuHMamHairlessRelIntChange, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_signif(comparisons = list(c("Su(H)::GFP","Hairless::GFP")), y_position = max(filter(SuHMamHairlessRelIntChange,Molecule == "Su(H)::GFP")$RelativeBandIntensity)+0.2, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Su(H)::GFP","Mam::GFP")), y_position = max(filter(SuHMamHairlessRelIntChange,Molecule == "Mam::GFP")$RelativeBandIntensity)+0.2, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Mam::GFP","Hairless::GFP")), y_position = max(filter(SuHMamHairlessRelIntChange,Molecule == "Mam::GFP")$RelativeBandIntensity)+1, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("RelIntBgChange plot of SuH vs Mam vs Hairless.jpg"), device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

#### end ####

#### Random Plots ####

## Generate plot for comparing background levels of each molecule. I need to check the imaging settings!
Labs <- c(title = " ", x = " ", y = "Nuclear fluorescence (AU)")
Background_crossbar2(SuHMamHairless, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_signif(comparisons = list(c("Su(H)::GFP","Hairless::GFP")), y_position = max(filter(SuHMamHairless,Molecule == "Su(H)::GFP")$Background)+0.2, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Su(H)::GFP","Mam::GFP")), y_position = max(filter(SuHMamHairless,Molecule == "Su(H)::GFP")$Background)+100, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Mam::GFP","Hairless::GFP")), y_position = max(filter(SuHMamHairless,Molecule == "Su(H)::GFP")$Background)+200, 
              map_signif_level=TRUE, test = "t.test", color = "black")

#### end ####

#### Stats on saved plots ####
subDir <- "Stats"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

Plots <- list()
for (p in 1:length(CombinedRelIntList)) {
  df <- CombinedRelIntList[[p]]
  Plot <- ggdensity(df, x = "RelativeBandIntensity", facet.by = "Molecule")
  Plots[[p]] <- Plot
  rm(Plot)
}
ggsave("Normal distribution plots of RelIntBg for molecule comparison.jpg", arrangeGrob(grobs = Plots, ncol = 2),
       device = "jpeg", dpi = "retina", width = 24,
       height = 9*ceiling(length(Plots)/2), units = "cm")
rm(Plots)

Plots <- list()
for (p in 1:length(CombinedRelIntChangeList)) {
  df <- CombinedRelIntChangeList[[p]]
  Plot <- ggdensity(df, x = "RelativeBandIntensity", facet.by = "Molecule")
  Plots[[p]] <- Plot
  rm(Plot)
}
ggsave("Normal distribution plots of RelIntChangeBg for molecule comparison.jpg", arrangeGrob(grobs = Plots, ncol = 2),
       device = "jpeg", dpi = "retina", width = 24,
       height = 9*ceiling(length(Plots)/2), units = "cm")
rm(Plots)

Plots <- list()
for (p in 1:length(CombinedRelIntChangeList)) {
  df <- CombinedRelIntChangeList[[p]]
  Plot <- ggqqplot(df, x = "RelativeBandIntensity", facet.by = "Molecule")
  Plots[[p]] <- Plot
  rm(Plot)
}
ggsave("QQ Plots of RelIntChangeBg for molecule comparison.jpg", arrangeGrob(grobs = Plots, ncol = 2),
       device = "jpeg", dpi = "retina", width = 24,
       height = 9*ceiling(length(Plots)/2), units = "cm")
rm(Plots)

Plots <- list()
for (p in 1:length(CombinedRelIntList)) {
  df <- CombinedRelIntList[[p]]
  Plot <- ggqqplot(df, x = "RelativeBandIntensity", facet.by = "Molecule")
  Plots[[p]] <- Plot
  rm(Plot)
}
ggsave("QQ Plots of RelIntBg for molecule comparison.jpg", arrangeGrob(grobs = Plots, ncol = 2),
       device = "jpeg", dpi = "retina", width = 24,
       height = 9*ceiling(length(Plots)/2), units = "cm")
rm(Plots)

SummRelIntBgList <- list()
for (p in 1:length(CombinedRelIntList)) {
  df <- na.omit(CombinedRelIntList[[p]])
  dfStatistics <- data.frame()
  for (g in 1:length(unique(df$Molecule))) {
    Testdf <- dplyr::filter(df,Molecule == unique(df$Molecule)[g])
    MeanValue <- mean(Testdf$RelativeBandIntensity)
    MedianValue <- median(Testdf$RelativeBandIntensity)
    VarValue <- var(Testdf$RelativeBandIntensity)
    StandardDeviation <- sd(Testdf$RelativeBandIntensity)
    StandardError <- StandardDeviation/sqrt(length(Testdf$RelativeBandIntensity))
    result <- shapiro.test(Testdf$RelativeBandIntensity)
    SWpvalue <- result[["p.value"]]
    StatsResults <- c(unique(df$Molecule)[g],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue)
    dfStatistics <- rbind(dfStatistics,StatsResults)
    rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
  }
  names(dfStatistics) <- c("Molecule","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value")
  SummRelIntBgList[[p]] <- dfStatistics
  names(SummRelIntBgList)[p] <- names(CombinedRelIntList[p])
  rm(df,dfStatistics)
}
SaveResult <- as.data.frame(SummRelIntBgList)
SaveResult <- t(SaveResult)
write.csv(SaveResult, file = "Summary of RelIntBg for different molecules.csv")

dfStatistics <- data.frame()
for (p in 1:length(CombinedRelIntList)) {
  df <- na.omit(CombinedRelIntList[[p]])
  Molecule <- names(CombinedRelIntList)[p]
  for (g in 1:length(unique(df$Molecule))) {
    df1 <- dplyr::filter(df,Molecule == unique(df$Molecule)[1])
    df2 <- dplyr::filter(df,Molecule == unique(df$Molecule)[2])
  }
  StudentT <- t.test(df1$RelativeBandIntensity, df2$RelativeBandIntensity, var.equal = TRUE)
  StudentT <- StudentT[["p.value"]]
  WelchT <- t.test(df1$RelativeBandIntensity, df2$RelativeBandIntensity, var.equal = FALSE)
  WelchT <- WelchT[["p.value"]]
  MannWhitney <- wilcox.test(df1$RelativeBandIntensity,df2$RelativeBandIntensity)
  MannWhitney <- MannWhitney[["p.value"]]
  SummaryVec <- c(Molecule, StudentT, WelchT, MannWhitney)
  dfStatistics <- rbind(dfStatistics, SummaryVec)
  rm(SummaryVec,df1,df2,df,Molecule, StudentT, WelchT, MannWhitney)
}
names(dfStatistics) <- c("Molecule","StudentT","WelchT","MannWhitney")
write.csv(dfStatistics, file = "SigStats of RelIntBg for different molecules compared to SuH.csv")
SigStatsRelIntBgSuH <- dfStatistics
rm(dfStatistics)

SummRelIntChangeBgList <- list()
for (p in 1:length(CombinedRelIntChangeList)) {
  df <- na.omit(CombinedRelIntChangeList[[p]])
  dfStatistics <- data.frame()
  for (g in 1:length(unique(df$Molecule))) {
    Testdf <- dplyr::filter(df,Molecule == unique(df$Molecule)[g])
    MeanValue <- mean(Testdf$RelativeBandIntensity)
    MedianValue <- median(Testdf$RelativeBandIntensity)
    VarValue <- var(Testdf$RelativeBandIntensity)
    StandardDeviation <- sd(Testdf$RelativeBandIntensity)
    StandardError <- StandardDeviation/sqrt(length(Testdf$RelativeBandIntensity))
    result <- shapiro.test(Testdf$RelativeBandIntensity)
    SWpvalue <- result[["p.value"]]
    StatsResults <- c(unique(df$Molecule)[g],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue)
    dfStatistics <- rbind(dfStatistics,StatsResults)
    rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
  }
  names(dfStatistics) <- c("Molecule","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value")
  SummRelIntChangeBgList[[p]] <- dfStatistics
  names(SummRelIntChangeBgList)[p] <- names(CombinedRelIntChangeList[p])
  rm(df,dfStatistics)
}
SaveResult <- as.data.frame(SummRelIntChangeBgList)
SaveResult <- t(SaveResult)
write.csv(SaveResult, file = "Summary of RelIntChangeBg for different molecules.csv")

dfStatistics <- data.frame()
for (p in 1:length(CombinedRelIntChangeList)) {
  df <- na.omit(CombinedRelIntChangeList[[p]])
  Molecule <- names(CombinedRelIntChangeList)[p]
  for (g in 1:length(unique(df$Molecule))) {
    df1 <- dplyr::filter(df,Molecule == unique(df$Molecule)[1])
    df2 <- dplyr::filter(df,Molecule == unique(df$Molecule)[2])
  }
  StudentT <- t.test(df1$RelativeBandIntensity, df2$RelativeBandIntensity, var.equal = TRUE)
  StudentT <- StudentT[["p.value"]]
  WelchT <- t.test(df1$RelativeBandIntensity, df2$RelativeBandIntensity, var.equal = FALSE)
  WelchT <- WelchT[["p.value"]]
  MannWhitney <- wilcox.test(df1$RelativeBandIntensity,df2$RelativeBandIntensity)
  MannWhitney <- MannWhitney[["p.value"]]
  SummaryVec <- c(Molecule, StudentT, WelchT, MannWhitney)
  dfStatistics <- rbind(dfStatistics, SummaryVec)
  rm(SummaryVec,df1,df2,df,Molecule, StudentT, WelchT, MannWhitney)
}
names(dfStatistics) <- c("Molecule","StudentT","WelchT","MannWhitney")
write.csv(dfStatistics, file = "SigStats of RelIntChangeBg for different molecules compared to SuH.csv")
SigStatsRelIntChangeBgSuH <- dfStatistics
rm(dfStatistics)

setwd(mainDir)
#### end ####