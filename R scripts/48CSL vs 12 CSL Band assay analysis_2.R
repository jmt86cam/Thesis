# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("ggplot2")
install.packages("cowplot")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggsignif")
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggsignif)
#### end ####

#### Data analysis functions ####
## Calculate the standard error for each row of a data frame 
SEMrows <- function(df){
  dft <- as.data.frame(t(df))
  rm(df)
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

# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/48 vs 12 CSL sites"
setwd(mainDir)
# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/48 vs 12 CSL sites/Background SuH values for 48 and 12 CSL sites experiment.csv", head = TRUE, sep = ",")
BackgroundValues <- data.frame(Nucleus=DF[,1], SuHBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Band assay/48 vs 12 CSL sites/csv files"

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
    FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
  } 
  rm(Path)
  rm(Files)
  df <- read.csv(FilePaths[[1]])
  names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[1]]), end=-5)
  for (i in 2:length(FilePaths)){
    Nextdf <- read.csv(FilePaths[[i]])
    df <- merge(df,Nextdf, by="X", all = TRUE)
    names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[i]]), end=-5)
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

SuHBgNorm <- list()
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    NormValue <- BackgroundValues$SuHBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  SuHBgNorm[[s]] <- df
  names(SuHBgNorm)[s] <- names(SuHRawData[s])
  rm(df,ColNames)
}

MeanminsixSubList <- list()
for (l in 1:length(SuHminsixSub)) {
  df <- SuHminsixSub[[l]]
  B1 <- select(df,ends_with("band 1"))
  B2 <- select(df,ends_with("band 2"))
  Band1 <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(B1, na.rm = TRUE))
  Band1$SEM <- SEMrows(B1)
  Band1$Genotype <- paste(names(SuHminsixSub)[l]," (N = ", ncol(B1), ")", sep = "")
  Band2 <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(B2, na.rm = TRUE))
  Band2$SEM <- SEMrows(B2)
  Band2$Genotype <- paste(names(SuHminsixSub)[l]," (N = ", ncol(B2), ")", sep = "")
  MeanminsixSubList[[l]] <- dplyr::bind_rows(list(Band1=Band1, Band2=Band2), .id = 'Band')
  names(MeanminsixSubList)[l] <- names(SuHminsixSub[l])
  rm(B1,B2,Band1,Band2)
}
MeanSuHminsixSub <- dplyr::bind_rows(MeanminsixSubList, .id = 'Name')
rm(MeanminsixSubList)

MeanRawDataList <- list()
for (l in 1:length(SuHRawData)) {
  df <- SuHRawData[[l]]
  B1 <- select(df,ends_with("band 1"))
  B2 <- select(df,ends_with("band 2"))
  Band1 <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(B1, na.rm = TRUE))
  Band1$SEM <- SEMrows(B1)
  Band1$Genotype <- paste(names(SuHRawData)[l]," (N = ", ncol(B1), ")", sep = "")
  Band2 <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(B2, na.rm = TRUE))
  Band2$SEM <- SEMrows(B2)
  Band2$Genotype <- paste(names(SuHRawData)[l]," (N = ", ncol(B2), ")", sep = "")
  MeanRawDataList[[l]] <- dplyr::bind_rows(list(Band1=Band1, Band2=Band2), .id = 'Band')
  names(MeanRawDataList)[l] <- names(SuHRawData[l])
  rm(B1,B2,Band1,Band2)
}
MeanSuHRawData <- dplyr::bind_rows(MeanRawDataList, .id = 'Name')
rm(MeanRawDataList)

MeanBgNormList <- list()
for (l in 1:length(SuHBgNorm)) {
  df <- SuHBgNorm[[l]]
  B1 <- select(df,ends_with("band 1"))
  B2 <- select(df,ends_with("band 2"))
  Band1 <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(B1, na.rm = TRUE))
  Band1$SEM <- SEMrows(B1)
  Band1$Genotype <- paste(names(SuHBgNorm)[l]," (N = ", ncol(B1), ")", sep = "")
  Band2 <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(B2, na.rm = TRUE))
  Band2$SEM <- SEMrows(B2)
  Band2$Genotype <- paste(names(SuHBgNorm)[l]," (N = ", ncol(B2), ")", sep = "")
  MeanBgNormList[[l]] <- dplyr::bind_rows(list(Band1=Band1, Band2=Band2), .id = 'Band')
  names(MeanBgNormList)[l] <- names(SuHBgNorm[l])
  rm(B1,B2,Band1,Band2)
}
MeanSuHBgNorm <- dplyr::bind_rows(MeanBgNormList, .id = 'Name')
rm(MeanBgNormList)

#### end ####

#### Calculate fold increase of band over background ####
## Here we will use the mean of the 5 middle values and divide by the mean of the 5 values at
## each edge of the band (10 values in total).
RelIntSuH <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntSuH <- rbind(RelIntSuH,RelInt)
  rm(RelInt, Band, Background, df)
}

RelIntSuHBg <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  for (t in 1:length(ColNames)) {
    Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
    Background <- BackgroundValues$SuHBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(SuHRawData[s])
    RelInt$Nucleus <- rownames(RelInt)
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  }
  RelIntSuHBg <- rbind(RelIntSuHBg,RelInt)
  rm(RelInt, Band, Background, df)
  rm(ColNames)
}

RelIntSuHEnds <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntSuHEnds <- rbind(RelIntSuHEnds,RelInt)
  rm(RelInt, Band, Background, df)
}

RelIntSuHMinSixSub <- NULL
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
  RelInt <- as.data.frame(Band-Background)
  RelInt$Genotype <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
  RelIntSuHMinSixSub <- rbind(RelIntSuHMinSixSub,RelInt)
  rm(RelInt, Band, Background, df,dfOrdered,c,Column,MinSix)
}

RelIntSuHBgSub <- NULL
for (s in 1:length(SuHRawData)) {
  df <- SuHRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- BackgroundValues$SuHBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band-Background)
    RelInt$Genotype <- names(SuHRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
    rm(Band, Background, BandSubBg)
    RelIntSuHBgSub <- rbind(RelIntSuHBgSub,RelInt)
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
    Background <- BackgroundValues$SuHBG[which(grepl(ColNames[[t]], BackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(SuHRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeSuHBg <- rbind(RelIntChangeSuHBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

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
  RelInt$Genotype <- names(SuHRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeSuH <- rbind(RelIntChangeSuH,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  plot.title = element_text(hjust=0.5, face = 'bold'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.line = element_line(colour = "black"), # add axis lines in black
  text = element_text(size=15)) # make text larger
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

#### Saved Plots ####
Labs <- c(title = " ", x = "Distance (um)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHminsixSub,Band %in% "Band2")
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) 
ggsave("12 vs 48 CSL MinSixSub.jpg", device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

Labs <- c(title = " ", x = "Distance (um)", y = "Su(H)::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,Band %in% "Band2")
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank())+
  scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0,2000))
ggsave("12 vs 48 CSL RawData.jpg", device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

Labs <- c(title = " ", x = " ", y = "Fold change of Su(H)::EGFP intensity")
df <- filter(RelIntSuHBg,Genotype %in% c("12xCSL","48xCSL"))
RelInt_crossbar(df,Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c("12xCSL","48xCSL")), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 9.0) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0,10))
ggsave("12 vs 48 CSL RelIntBg.jpg", device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

Labs <- c(title = " ", x = " ", y = "Relative change of Su(H)::EGFP intensity")
df <- filter(RelIntChangeSuHBg,Genotype %in% c("12xCSL","48xCSL"))
RelInt_crossbar(df,Labs) + theme(legend.position = "none") +
  geom_signif(comparisons = list(c("12xCSL","48xCSL")), 
              map_signif_level=TRUE, test = "t.test", color = "black", y_position = 7.5) +
  scale_y_continuous(breaks = seq(0, 8, by = 2), limits = c(0,8))
ggsave("12 vs 48 CSL RelIntChangeBg.jpg", device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")
#### end ####

#### Random Plots ####

## Plots of the raw data
Labs <- c(title = "Raw Data", x = "Distance (um)", y = "SuH::EGFP Fluorescence (AU)")
df <- filter(MeanSuHRawData,Band %in% "Band2")
band_plot(df, Labs)

## Plots of the BgNorm data

Labs <- c(title = "Fold enrichment above background", x = "Distance (um)", y = "Normalised SuH Fluorescence")
df <- filter(MeanSuHBgNorm,Band %in% "Band2")
band_plot(df, Labs)

Labs <- c(title = "Fold enrichment above background", x = "Distance (um)", y = "Normalised SuH Fluorescence")
df <- SuHBgNorm[["12xCSL"]]
B2 <- select(df,ends_with("band 2"))
B2$Distance <- df$`Distance (um)`
B2 <- B2 %>% tidyr::gather("id", "value", 1:5)
ggplot(B2, aes(Distance,value)) + theme_jmt() +
  geom_point()+
  facet_wrap(~id) +
  labs(title = Labs[[1]], x = Labs[[2]], y = Labs[[3]]) +
  scale_color_brewer(palette = "Set2")
rm(df, B2)

## RelInt plots

Labs <- c(title = "Relative intensity of band", x = "Genotype", y = "Relative SuH intensity")
df <- filter(RelIntSuH,Genotype %in% c("12xCSL","48xCSL"))
RelInt_boxviolin(df,Labs)
