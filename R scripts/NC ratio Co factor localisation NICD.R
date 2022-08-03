# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("ggplot2")
install.packages("dplyr")
install.packages("stringr")
install.packages("forcats")
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
#### end ####

# Co factor to be analysed
CoFactor <- "NICD"
CoFactorLabel <- "NICD::mCherry"
# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/SuH localisation/Gland image quantifications"
setwd(mainDir)
# Create the path where the data is stored
DataPath <- paste(mainDir, CoFactor, sep = "/")
# Set the X axis order
XAxisOrder <- c("No construct controls","wtNDECDmScarlet","EAAAKA","4ALEA4","NECDmCherryNICD",
                "CD8 TEV mCherry NICD", "OldBLITz", "NewBLITz kept in dark", "NewBLITz kept in light",
                "PDGFR NICD","PDGFR 4ALEA4 NICD")

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
  FilePaths[sapply(FilePaths, is.null)] <- NULL
  rm(Path)
  rm(Files)
  Data <- list()
  for (i in 1:length(FilePaths)){
    df <- read.csv(FilePaths[[i]])
    dfOdd <- df[1:nrow(df)%%2!=0,]
    for (c in 2:length(colnames(dfOdd))) {
      names(dfOdd)[names(dfOdd) == colnames(dfOdd)[c]] <- paste("Cell", str_sub(colnames(dfOdd)[c]))
    }
    rm(c)
    dfOdd$Cell <- ceiling(dfOdd$X/2)
    dfOdd$X <- NULL
    dfEven <- df[1:nrow(df)%%2==0,]
    for (c in 2:length(colnames(dfEven))) {
      names(dfEven)[names(dfEven) == colnames(dfEven)[c]] <- paste("Nucleus", str_sub(colnames(dfEven)[c]))
    }
    dfEven$Cell <- ceiling(dfEven$X/2)
    dfEven$X <- NULL
    df <- merge(dfOdd, dfEven, by = "Cell")
    rm(dfOdd,dfEven,c)
    df$'Cytoplasm Area' <- df$`Cell Area` - df$`Nucleus Area`
    df$'Cytoplasm IntDen' <- df$`Cell IntDen` - df$`Nucleus IntDen`
    df$'Cytoplasm Mean' <- df$`Cytoplasm IntDen` / df$`Cytoplasm Area`
    df$NCratio <- df$`Nucleus Mean` / df$`Cytoplasm Mean`
    df$Gland <- str_sub(basename(FilePaths[[i]]), end=-5)
    df$IncubationTime <- FolderNames[f]
    Data[[i]] <- df
    rm(df)
  }
  Data <- do.call(rbind.data.frame, Data)
  SuHRawData[[f]] <- Data
  rm(Data,i)
  rm(FilePaths)
  names(SuHRawData)[f] <- FolderNames[f] 
}

rm(FolderNames,f)

SuHdata <- do.call(rbind.data.frame, SuHRawData)
SuHdata <- dplyr::select(SuHdata,c("Cell","Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio","Gland","IncubationTime"))
SuHdata <- SuHdata %>% mutate(IncubationTime = fct_relevel(IncubationTime, XAxisOrder)) 

ALEA <- SuHdata %>% dplyr::filter(IncubationTime == "4ALEA4")
Control <- SuHdata %>% dplyr::filter(IncubationTime == "No construct controls")
EAAAKA <- SuHdata %>% dplyr::filter(IncubationTime == "EAAAKA")
WT <- SuHdata %>% dplyr::filter(IncubationTime == "wtNDECDmScarlet")
OldBLITz <- SuHdata %>% dplyr::filter(IncubationTime == "OldBLITz")
CD8 <- SuHdata %>% dplyr::filter(IncubationTime == "CD8 TEV mCherry NICD")
NewBLITzDark <- SuHdata %>% dplyr::filter(IncubationTime == "NewBLITz kept in dark")
NewBLITzLight <- SuHdata %>% dplyr::filter(IncubationTime == "NewBLITz kept in light")
PDGFR <- SuHdata %>% dplyr::filter(IncubationTime == "PDGFR NICD")
PDGFRALEA <- SuHdata %>% dplyr::filter(IncubationTime == "PDGFR 4ALEA4 NICD")

#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  text = element_text(family = "sans", size = 15), # set default font to be Arial
  plot.title = element_text(hjust=0.5, face = 'bold'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  # tilt the x axis
  axis.line = element_line(colour = "black")) # add axis lines in black
}

NCratio_crossbar <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,NCratio,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

NucleusMean_crossbar <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,`Nucleus Mean`,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

CytoplasmMean_crossbar <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,`Cytoplasm Mean`,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

#### end ####

#### Saved plots ####
subDir <- paste("Plots for", CoFactor)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

df1 <- rbind(ALEA,EAAAKA,WT,OldBLITz,CD8)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactorLabel))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("4ALEA4","wtNDECDmScarlet")),y_position = 5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("EAAAKA","wtNDECDmScarlet")), y_position = 4.5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("OldBLITz","wtNDECDmScarlet")),y_position = 6, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","wtNDECDmScarlet")),y_position = 5.5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","OldBLITz")),y_position = 1.5, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("WT vs original linkers NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

df1 <- rbind(ALEA,EAAAKA,WT)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactorLabel))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("4ALEA4","wtNDECDmScarlet")),y_position = 1, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("4ALEA4","EAAAKA")), y_position = 0.75,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("EAAAKA","wtNDECDmScarlet")), y_position = 0.5,
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("4ALEA4 vs EAAAKA NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
width = 15, height = 15, units = "cm")

df1 <- rbind(PDGFR,PDGFRALEA)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactorLabel))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("PDGFR 4ALEA4 NICD","PDGFR NICD")), y_position = 5,
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("PDGFR vs PDGFR4ALEA4 NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

df1 <- rbind(NewBLITzDark,NewBLITzLight)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactorLabel))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("NewBLITz kept in dark","NewBLITz kept in light")),y_position = 0.6, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("NewBLITz kept in dark vs light NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

df1 <- rbind(WT,OldBLITz)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactorLabel))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("OldBLITz","wtNDECDmScarlet")),y_position = 0.6, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("OldBLITz vs WT NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

df1 <- rbind(WT,CD8)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactorLabel))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","wtNDECDmScarlet")),y_position = 0.6, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("CD8 TEV vs WT NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

setwd(mainDir)
#### end ####

#### Statistics Tests and Plots ####
subDir <- paste("Stats for", CoFactor)
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

NormDis_plots(SuHRawData, "NCratio", "IncubationTime",paste("Normal distribution of NCratio plots for",CoFactor))

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

QQ_plots(SuHRawData, "NCratio", "IncubationTime",paste("QQ Plots of",CoFactor))

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
    for (g in 1:length(unique(df$IncubationTime))) {
      Testdf <- dplyr::filter(df,IncubationTime == unique(df$IncubationTime)[g])
      MeanValue <- mean(Testdf$NCratio)
      MedianValue <- median(Testdf$NCratio)
      VarValue <- var(Testdf$NCratio)
      StandardDeviation <- sd(Testdf$NCratio)
      StandardError <- StandardDeviation/sqrt(length(Testdf$NCratio))
      result <- shapiro.test(Testdf$NCratio)
      SWpvalue <- result[["p.value"]]
      StatsResults <- c(unique(df$IncubationTime)[g],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue)
      dfStatistics <- rbind(dfStatistics,StatsResults)
      rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
    }
    names(dfStatistics) <- c("IncubationTime","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value")
    OutputList[[p]] <- dfStatistics
    names(OutputList)[p] <- names(DataList[p])
    rm(df,dfStatistics)
  }
  SaveResult <- as.data.frame(OutputList)
  SaveResult <- t(SaveResult)
  write.csv(SaveResult, file = paste(Filename,".csv"))
  return(OutputList)
}

SumRelIntSuH <- DataStatSumTestListFun(SuHRawData,paste("Summary of",CoFactor))

## Finally we want to generate a set of p-values for the data which can be checked for
## significance and compared to the summaries generated above to know which is most
## relevant. For normal data with equal variance it is the t test, unequal variance use
## the Welch's t test. If the data is nonparametric (a.k.a. not normal), then it is
## the Mann-Whitney test.

Stats_tests <- function(DataList, Filename){
  dfStatistics <- data.frame()
  for (p in 1:(length(DataList))) {
    df1 <- na.omit(DataList[[p]])
    for (q in 1:(length(DataList))) { 
      df2 <- na.omit(DataList[[q]])
      Comparison <- paste(unique(df1$IncubationTime), "vs", unique(df2$IncubationTime))
      StudentT <- t.test(df1$NCratio, df2$NCratio, var.equal = TRUE)
      StudentT <- StudentT[["p.value"]]
      WelchT <- t.test(df1$NCratio, df2$NCratio, var.equal = FALSE)
      WelchT <- WelchT[["p.value"]]
      MannWhitney <- wilcox.test(df1$NCratio,df2$NCratio)
      MannWhitney <- MannWhitney[["p.value"]]
      SummaryVec <- c(Comparison, StudentT, WelchT, MannWhitney)
      dfStatistics <- rbind(dfStatistics, SummaryVec)
      rm(SummaryVec,df2,Comparison, StudentT, WelchT, MannWhitney)
    }
    rm(df1,q)
  }
  names(dfStatistics) <- c("Comparison","StudentT","WelchT","MannWhitney")
  write.csv(dfStatistics, file = paste(Filename,".csv"))
  return(dfStatistics)
}

SigStatsRelIntSuH <- Stats_tests(SuHRawData,paste("SigStats of",CoFactor))

setwd(mainDir)
#### end ####

#### Useful plots ####
subDir <- paste("Useful plots for", CoFactor)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

df1 <- rbind(ALEA,EAAAKA,WT,OldBLITz,CD8)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactor))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("4ALEA4","wtNDECDmScarlet")),y_position = 5.5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("4ALEA4","EAAAKA")), y_position = 4,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","EAAAKA")), y_position = 4.5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","4ALEA4")), y_position = 5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("EAAAKA","wtNDECDmScarlet")), y_position = 5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("OldBLITz","wtNDECDmScarlet")),y_position = 6.5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","wtNDECDmScarlet")),y_position = 6, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("CD8 TEV mCherry NICD","OldBLITz")),y_position = 2, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("All linkers NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

df1 <- rbind(WT,OldBLITz,NewBLITzDark,NewBLITzLight)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactor))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("OldBLITz","wtNDECDmScarlet")),y_position = 4.5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("NewBLITz kept in dark","wtNDECDmScarlet")),y_position = 5.5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("NewBLITz kept in dark","OldBLITz")),y_position = 5, 
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("NewBLITz kept in dark","NewBLITz kept in light")),y_position = 4.5, 
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("All optogenetics NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

df1 <- rbind(WT,PDGFR,PDGFRALEA)
Labs <- c(title = "", x = "Incubation Time", y = paste("Nuclear:Cytoplasmic ratio of",CoFactor))
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) +
  geom_signif(comparisons = list(c("PDGFR 4ALEA4 NICD","PDGFR NICD")), y_position = 6,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("PDGFR 4ALEA4 NICD","wtNDECDmScarlet")), y_position = 6.5,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("wtNDECDmScarlet","PDGFR NICD")), y_position = 5.5,
              map_signif_level=TRUE, test = "t.test", color = "black")
ggsave(paste("All PDGFR NCratio of",CoFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

setwd(mainDir)
#### end ####