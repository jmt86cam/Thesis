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

# Set location of data and where plots will be saved to:
setwd("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/Nuclear levels csv files")
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/Nuclear levels csv files"

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
rm(f)

NICDRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    if (grepl("NICD",Files[i])) {
      FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
    }
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
  NICDRawData[[f]] <- Data
  rm(Data, i)
  rm(FilePaths)
  names(NICDRawData)[f] <- FolderNames[f] 
}

XAxisOrder <- str_sort(FolderNames, numeric = TRUE)

rm(FolderNames,f)


NICDdata <- do.call(rbind.data.frame, NICDRawData)
NICDdata <- dplyr::select(NICDdata,c("Cell","Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio","Gland","IncubationTime"))
NICDdata <- NICDdata %>% mutate(IncubationTime = fct_relevel(IncubationTime, XAxisOrder)) 

SuHdata <- do.call(rbind.data.frame, SuHRawData)
SuHdata <- dplyr::select(SuHdata,c("Cell","Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio","Gland","IncubationTime"))
SuHdata <- SuHdata %>% mutate(IncubationTime = fct_relevel(IncubationTime, XAxisOrder)) 

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
DarkControl <- NICDdata %>% dplyr::filter(IncubationTime == "Dark control")
Light2 <- NICDdata %>% dplyr::filter(IncubationTime == "2 hours light")
Light5 <- NICDdata %>% dplyr::filter(IncubationTime == "5 hours light")
Light8 <- NICDdata %>% dplyr::filter(IncubationTime == "8 hours light")
Light12 <- NICDdata %>% dplyr::filter(IncubationTime == "12 hours light")
Light24 <- NICDdata %>% dplyr::filter(IncubationTime == "24 hours light")

df1 <- rbind(DarkControl,Light2,Light5,Light8,Light12,Light24)

Labs <- c(title = "", x = "Incubation Time", y = "Nuclear:Cytoplasmic ratio of BLITz-mCherry::LOV::NICD")
NCratio_crossbar(df1, Labs) + theme(axis.title.x=element_blank(), legend.position = "none") + guides(colour = guide_legend(nrow=2, byrow = TRUE)) 
ggsave("NCratio of NICDmCherry.jpg", device = "jpeg", dpi = "retina",
width = 15, height = 15, units = "cm")

Labs <- c(title = "Mean Nuclear fluorescence", x = "Incubation Time", y = "Su(H)::GFP fluorescence (AU)")
NucleusMean_crossbar(df1, Labs)

Labs <- c(title = "Mean Cytoplasm fluorescence", x = "Incubation Time", y = "Su(H)::GFP fluorescence (AU)")
CytoplasmMean_crossbar(df1, Labs)

