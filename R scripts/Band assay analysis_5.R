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
library(tidyverse)
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
    dfmeanse$IncubationTime <- paste(names(DataFrame)[l]," (N = ", Nnumber, ")", sep = "")
    MeanList[[l]] <- dfmeanse
    names(MeanList)[l] <- names(DataFrame[l])
    rm(dfmeanse,Nnumber)
  }
  Output <- dplyr::bind_rows(MeanList, .id = 'Name')
  rm(MeanList)
  return(Output)
}
#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator"
setwd(mainDir)
# Load the background values
BackgroundValues <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/SuHGFP and NICDmScarlet nuclear background signals.csv", head = TRUE, sep = ",")
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/Band assay csv files"

#### Read in and normalise the data ####
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
SuHRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    if (grepl("SuH",Files[i])) {
      FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
    }
  } 
  FilePaths <- FilePaths[-which(sapply(FilePaths, is.null))]
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
MeanSuHRawData <- MeanDF(SuHRawData)
SuHXAxisOrderN <- str_sort(unique(MeanSuHRawData$IncubationTime), numeric = TRUE)
MeanSuHRawData <- MeanSuHRawData %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrderN)) 

NICDRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    if (grepl("NICD",Files[i]) & !grepl("not good",Files[i])) {
      FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
    }
  } 
  FilePaths <- FilePaths[-which(sapply(FilePaths, is.null))]
  rm(Path)
  rm(Files)
  df <- read.csv(FilePaths[[1]])
  names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[1]]), end=-10)
  for (i in 2:length(FilePaths)){
    Nextdf <- read.csv(FilePaths[[i]])
    df <- merge(df,Nextdf, by="X", all = TRUE)
    names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[i]]), end=-10)
    rm(Nextdf)
  }
  names(df)[names(df) == 'X'] <- "Distance (um)"
  NICDRawData[[f]] <- df
  rm(df)
  rm(FilePaths)
  names(NICDRawData)[f] <- FolderNames[f] 
}
MeanNICDRawData <- MeanDF(NICDRawData)
NICDXAxisOrderN <- str_sort(unique(MeanNICDRawData$IncubationTime), numeric = TRUE)
MeanNICDRawData <- MeanNICDRawData %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrderN)) 

SuHminsixNorm <- list()
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  SuHminsixNorm[[s]] <- df
  names(SuHminsixNorm)[s] <- names(SuHRawData[s])
  rm(df,ColNames)
}
MeanSuHMinSix <- MeanDF(SuHminsixNorm)
MeanSuHMinSix <- MeanSuHMinSix %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrderN)) 

SuHBgNorm <- list()
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    NormValue <- BackgroundValues$Su.H...GFP.nuclear.background[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  SuHBgNorm[[s]] <- df
  names(SuHBgNorm)[s] <- names(SuHRawData[s])
  rm(df,ColNames)
}
MeanSuHBgNorm <- MeanDF(SuHBgNorm)
MeanSuHBgNorm <- MeanSuHBgNorm %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrderN)) 

NICDminsixNorm <- list()
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  NICDminsixNorm[[s]] <- df
  names(NICDminsixNorm)[s] <- names(NICDRawData[s])
  rm(df,ColNames)
}
MeanNICDMinSix <- MeanDF(NICDminsixNorm)
MeanNICDMinSix <- MeanNICDMinSix %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrderN)) 

NICDBgNorm <- list()
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    NormValue <- BackgroundValues$NICD..mScarlet.nuclear.background[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  NICDBgNorm[[s]] <- df
  names(NICDBgNorm)[s] <- names(NICDRawData[s])
  rm(df,ColNames)
}
MeanNICDBgNorm <- MeanDF(NICDBgNorm)
MeanNICDBgNorm <- MeanNICDBgNorm %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrderN)) 

rm(FolderNames,NICDXAxisOrderN,SuHXAxisOrderN,s,t)
#### end ####

#### Calculate fold increase of band over background ####
## Here we will use the mean of the 5 middle values and divide by the mean of the 6 minimum values
## across the band
SuHXAxisOrder <- str_sort(names(SuHRawData), numeric = TRUE)
NICDXAxisOrder <- str_sort(names(NICDRawData), numeric = TRUE)

RelIntSuH <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
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
  RelInt$IncubationTime <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntSuH <- rbind(RelIntSuH,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
RelIntSuH <- RelIntSuH %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrder)) 

RelIntNICD <- NULL
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
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
  RelInt$IncubationTime <- names(NICDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntNICD <- rbind(RelIntNICD,RelInt)
  rm(RelInt, Band, Background, df,dfOrdered,c,Column,MinSix)
}
RelIntNICD <- RelIntNICD %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrder)) 

## Here we will use the mean of the 5 middle values and divide by the mean of the 5 values at
## each edge of the band (10 values in total).
RelIntSuHEdge <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$IncubationTime <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntSuHEdge <- rbind(RelIntSuHEdge,RelInt)
  rm(RelInt, Band, Background, df)
}
RelIntSuHEdge <- RelIntSuHEdge %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrder)) 

RelIntNICDEdge <- NULL
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$IncubationTime <- names(NICDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntNICDEdge <- rbind(RelIntNICDEdge,RelInt)
  rm(RelInt, Band, Background, df)
}
RelIntNICDEdge <- RelIntNICDEdge %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrder)) 

## Here we will use the mean of the 5 middle values and divide by the background values
RelIntSuHBg <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$Su.H...GFP.nuclear.background[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$IncubationTime <- names(SuHRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background)
    RelIntSuHBg <- rbind(RelIntSuHBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntSuHBg <- RelIntSuHBg %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrder)) 

RelIntNICDBg <- NULL
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$NICD..mScarlet.nuclear.background[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$IncubationTime <- names(NICDRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background)
    RelIntNICDBg <- rbind(RelIntNICDBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntNICDBg <- RelIntNICDBg %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrder)) 

## Here we will calculate the relative change rather than fold change, using the mean of the 
## 5 middle values and the minimum six values across the band
RelIntChangeSuH <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
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
  RelInt$IncubationTime <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeSuH <- rbind(RelIntChangeSuH,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
RelIntChangeSuH <- RelIntSuHBg %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrder)) 

RelIntChangeNICD <- NULL
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
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
  RelInt$IncubationTime <- names(NICDRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeNICD <- rbind(RelIntChangeNICD,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
RelIntChangeNICD <- RelIntNICDBg %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrder)) 

## Here we will calculate the relative change rather than fold change, using the mean of the 
## 5 middle values and the background values
RelIntChangeSuHBg <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$Su.H...GFP.nuclear.background[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$IncubationTime <- names(SuHRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeSuHBg <- rbind(RelIntChangeSuHBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntChangeSuHBg <- RelIntChangeSuHBg %>% mutate(IncubationTime = fct_relevel(IncubationTime, SuHXAxisOrder)) 

RelIntChangeNICDBg <- NULL
for (s in 1:length(NICDRawData)) {
  df <- NICDRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$NICD..mScarlet.nuclear.background[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$IncubationTime <- names(NICDRawData[s])
    RelInt$Nucleus <- rownames(RelInt)
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeNICDBg <- rbind(RelIntChangeNICDBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntChangeNICDBg <- RelIntChangeNICDBg %>% mutate(IncubationTime = fct_relevel(IncubationTime, NICDXAxisOrder)) 

#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  text = element_text(family = "sans",size=15), # set default font to be Arial
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
  ggplot(Data, aes(Distance,Mean,group = IncubationTime,col=IncubationTime,shape=IncubationTime,ymin=Mean-SEM,ymax=Mean+SEM)) + theme_jmt() +
    geom_ribbon(linetype = "dashed", fill = "grey90",alpha = 0.5, size = 0.25) + geom_line(size = 0.5) + geom_point() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    scale_color_brewer(palette = "Set2")
}

## Define a function to plot the relative intensities of the bands
RelInt_plot <- function(Data, Labels){
  ggplot(Data, aes(IncubationTime,RelativeBandIntensity,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_violin() + geom_boxplot(width=0.1, show.legend = FALSE,colour = "Grey60",alpha = 0.5, outlier.shape = NA) + geom_jitter() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

RelInt_crossbar <- function(Data, Labels){
  ggplot(Data,aes(IncubationTime,RelativeBandIntensity,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}
#### end ####

#### Saved Plots ####

## Plots of increasing light incubations
subDir <- "Increasing light incubations"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

Labs <- c(title = " ", x = "Distance (μm)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,Name %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0,900))
ggsave(paste("RawData Band plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = " ", x = "Distance (μm)", y = "NICD::mCherry Fluorescence (AU)")
df <- filter(MeanNICDRawData,Name %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE))
ggsave(paste("RawData NICD Band plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Fold change of Su(H)::EGFP intensity")
df <- filter(RelIntSuHBg,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2), limits = c(0,9)) 
ggsave(paste("Fold SuH over background plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Fold change of NICD::mCherry intensity")
df <- filter(RelIntNICDBg,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Fold NICD over background plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuHBg,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) + guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0,5.35))
ggsave(paste("Relative change SuH over background plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuHBg,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) + guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0,5.35))  +
  geom_signif(comparisons = list(c("2 hours light","5 hours light")), y_position = 3.5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("2 hours light","8 hours light")), y_position = 4,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("2 hours light","12 hours light")), y_position = 4.5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("2 hours light","24 hours light")), y_position = 5,
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("Relative change SuH over background plot of 2, 5, 8, 12 and 24 hours in incubator PLUS STAT.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of NICD::mCherry intensity")
df <- filter(RelIntChangeNICDBg,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change NICD over background plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuH,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change SuH over minsix plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of NICD::mCherry intensity")
df <- filter(RelIntChangeNICD,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change NICD over minsix plot of 2, 5, 8, 12 and 24 hours in incubator.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

setwd(mainDir)

## Plots of two hours followed by dark
subDir <- "Two hours followed by dark"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

Labs <- c(title = " ", x = "Distance (μm)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,Name %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0,900))
ggsave(paste("RawData Band plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = " ", x = "Distance (μm)", y = "NICD::mCherry Fluorescence (AU)")
df <- filter(MeanNICDRawData,Name %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE))
ggsave(paste("RawData NICD Band plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Fold change of Su(H)::EGFP intensity")
df <- filter(RelIntSuHBg,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Fold SuH over background plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Fold change of NICD::mCherry intensity")
df <- filter(RelIntNICDBg,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Fold NICD over background plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuHBg,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0,5.35)) 
ggsave(paste("Relative change SuH over background plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of NICD::mCherry intensity")
df <- filter(RelIntChangeNICDBg,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change NICD over background plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuH,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change SuH over minsix plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of NICD::mCherry intensity")
df <- filter(RelIntChangeNICD,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change NICD over minsix plot of 2 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

setwd(mainDir)

## Plots of five hours followed by dark
subDir <- "Five hours followed by dark"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

Labs <- c(title = " ", x = "Distance (μm)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,Name %in% c("5 hours light","5 hours light 2 hours dark"))
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE))
ggsave(paste("RawData Band plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = " ", x = "Distance (μm)", y = "NICD::mCherry Fluorescence (AU)")
df <- filter(MeanNICDRawData,Name %in% c("5 hours light","5 hours light 2 hours dark"))
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE))
ggsave(paste("RawData NICD Band plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Fold change of Su(H)::EGFP intensity")
df <- filter(RelIntSuHBg,IncubationTime %in% c("5 hours light","5 hours light 2 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Fold SuH over background plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Fold change of NICD::mCherry intensity")
df <- filter(RelIntNICDBg,IncubationTime %in% c("5 hours light","5 hours light 2 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Fold NICD over background plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuHBg,IncubationTime %in% c("5 hours light","5 hours light 2 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change SuH over background plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of NICD::mCherry intensity")
df <- filter(RelIntChangeNICDBg,IncubationTime %in% c("5 hours light","5 hours light 2 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change NICD over background plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuH,IncubationTime %in% c("5 hours light","5 hours light 2 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change SuH over minsix plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

Labs <- c(title = "", x = "", y = "Relative change of NICD::mCherry intensity")
df <- filter(RelIntChangeNICD,IncubationTime %in% c("5 hours light","5 hours light 2 hours dark"))
RelInt_crossbar(df,Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave(paste("Relative change NICD over minsix plot of 5 hours in incubator then dark.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

setwd(mainDir)
## Plots of each incubation time showing the date and gland.
df <- RelIntSuHBg
df$Date <- str_sub(df$Nucleus,1,6)
df$Gland <- str_sub(df$Nucleus,-4,-3)
Labs <- c(title = "", x = "", y = "Relative Su(H)::EGFP intensity")
ggplot(df, aes(Date,RelativeBandIntensity,col=Gland,shape=Gland)) + theme_jmt() +
  geom_point() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  facet_wrap(~IncubationTime,ncol = 3) +
  labs(title = Labs[[1]], x = Labs[[2]], y = Labs[[3]]) +
  scale_color_brewer(palette = "Set2")
ggsave("Separated Incubation Times point.pdf", device = "pdf", width = 30, height = 5*(ceiling(length(unique(df$IncubationTime))/3)), units = "cm", dpi = "retina", limitsize = FALSE)

ggplot(df, aes(Date,RelativeBandIntensity,col=Gland,shape=Gland)) + theme_jmt() +
  geom_jitter() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  facet_wrap(~IncubationTime,ncol = 3) +
  labs(title = Labs[[1]], x = Labs[[2]], y = Labs[[3]]) +
  scale_color_brewer(palette = "Set2")
ggsave("Separated Incubation Times jitter.pdf", device = "pdf", width = 30, height = 5*(ceiling(length(unique(df$IncubationTime))/3)), units = "cm", dpi = "retina", limitsize = FALSE)
rm(df)
#### end ####


## Plots of the raw data

Labs <- c(title = "Raw Data", x = "Distance (um)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
band_plot(df, Labs)

Labs <- c(title = "Raw Data", x = "Distance (um)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,IncubationTime %in% c("2 hours light","5 hours light","8 hours light"))
band_plot(df, Labs)

Labs <- c(title = "Raw Data", x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
df <- filter(MeanNICDRawData,IncubationTime %in% c("2 hours light","5 hours light","8 hours light"))
band_plot(df, Labs)

## Plots of the MinSix data

Labs <- c(title = "Normalised to minimum six values", x = "Distance (um)", y = "Normalised Su(H) Fluorescence)")
df <- filter(MeanSuHMinSix,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to minimum six values", x = "Distance (um)", y = "Normalised Su(H) Fluorescence")
df <- filter(MeanSuHMinSix,IncubationTime %in% c("2 hours light","5 hours light","8 hours light"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to minimum six values", x = "Distance (um)", y = "Normalised NICD Fluorescence")
df <- filter(MeanNICDMinSix,IncubationTime %in% c("2 hours light","5 hours light","8 hours light"))
band_plot(df, Labs)

## Plots of the BgNorm data

Labs <- c(title = "Normalised to background value", x = "Distance (um)", y = "Normalised Su(H) Fluorescence")
df <- filter(MeanSuHBgNorm,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to background value", x = "Distance (um)", y = "Normalised Su(H) Fluorescence")
df <- filter(MeanSuHBgNorm,IncubationTime %in% c("2 hours light","5 hours light","8 hours light"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to background value", x = "Distance (um)", y = "Normalised NICD Fluorescence")
df <- filter(MeanNICDBgNorm,IncubationTime %in% c("2 hours light","5 hours light","8 hours light"))
band_plot(df, Labs)

## RelInt plots

Labs <- c(title = "Relative intensity of band", x = "Incubation Time", y = "Relative Su(H) intensity")
df <- filter(RelIntSuH,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_plot(df,Labs)

Labs <- c(title = "Relative intensity of band", x = "Incubation Time", y = "Relative NICD intensity")
df <- filter(RelIntNICD,IncubationTime %in% c("2 hours light","5 hours light","8 hours light","12 hours light","24 hours light"))
RelInt_plot(df,Labs)

## Plots of dark period after 2 hours of light

Labs <- c(title = "Raw Data", x = "Distance (um)", y = "Su(H)::GFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs)

Labs <- c(title = "Raw Data", x = "Distance (um)", y = "NICD::mScarlet Fluorescence (AU)")
df <- filter(MeanNICDRawData,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to minimum six values", x = "Distance (um)", y = "Normalised Su(H) Fluorescence")
df <- filter(MeanSuHMinSix,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to minimum six values", x = "Distance (um)", y = "Normalised NICD Fluorescence")
df <- filter(MeanNICDMinSix,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to background value", x = "Distance (um)", y = "Normalised Su(H) Fluorescence")
df <- filter(MeanSuHBgNorm,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs)

Labs <- c(title = "Normalised to background value", x = "Distance (um)", y = "Normalised NICD Fluorescence")
df <- filter(MeanNICDBgNorm,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
band_plot(df, Labs)

Labs <- c(title = "Relative intensity of band", x = "Incubation Time", y = "Relative Su(H) intensity")
df <- filter(RelIntSuH,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_plot(df,Labs)

Labs <- c(title = "Relative intensity of band", x = "Incubation Time", y = "Relative NICD intensity")
df <- filter(RelIntNICD,IncubationTime %in% c("2 hours light","2 hours light 2 hours dark","2 hours light 3.5 hours dark","2 hours light 5 hours dark"))
RelInt_plot(df,Labs)