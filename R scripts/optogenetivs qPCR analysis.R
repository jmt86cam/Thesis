# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("tidyr")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("forcats")
install.packages("dplyr")
library(tidyr)
library(ggplot2)
library(cowplot)
library(forcats)
library(dplyr)
#### end ####

# Set where plots will be saved to:
setwd("/Users/jonathantownson/Documents/PhD/qPCR data/Optogenetic experiments")
# Set the path to the folder where the data csv files are stored:
DataPath <- "/Users/jonathantownson/Documents/PhD/qPCR data/Optogenetic experiments/Cleaned data for R"
# Instruct the script what the control primers being used are
ControlPrimers <- "RP49"
# Instruct the script of any Genotype ordering required e.g. control then experimental condition
XAxisOrder <- c("NoCRY Dark", "NoCRY 24hr Light", "Opto Dark", "Opto 1hr Light", "Opto 2hr Light", "Opto 8hr Light", "Opto 24hr Light")

#### Read in the data ####
## Create the folder names for the data
FolderFiles <- list.files(DataPath, pattern = "*.csv")
## Read in all the files in the folder
Files <- list()
# Comment out one of the below, first for loop only keeps files where StDev <1
# Second for loop keeps all files to be merged into raw data
#for (f in 1:length(FolderFiles)){
  df <- read.csv(paste(DataPath,FolderFiles[f], sep ="/"))
  df <- dplyr::filter(df,Cp < 40)
  TechRepNum <- count(df, Name, Primers)
  ## Calculate mean and StDev of technical replicates
  MeanCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),data=df,FUN = mean)
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'Cp'
  StDevCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),FUN = sd)
  MeanCp <- merge(MeanCp,StDevCp,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'StDevCp'
  MeanCp <- merge(MeanCp,TechRepNum,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'n'] <- 'TechRepNum'
  rm(StDevCp)
  rm(TechRepNum)
  rm(df)
  MeanCp$StDevCp[is.na(MeanCp$StDevCp)] <- 0
  ## Select data with StDev<1 and Cp less than 40
  Data <- dplyr::filter(MeanCp,StDevCp < 1)
  Data <- subset(Data, select = -StDevCp)
  Data$qPCRdate <- FolderFiles[[f]]
  Files[[f]] <- Data
  rm(MeanCp)
  rm(Data)

for (f in 1:length(FolderFiles)){
  df <- read.csv(paste(DataPath,FolderFiles[f], sep ="/"))
  df <- dplyr::filter(df,Cp < 40)
  TechRepNum <- count(df, Name, Primers)
  ## Calculate mean and StDev of technical replicates
  MeanCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),data=df,FUN = mean)
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'Cp'
  StDevCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),FUN = sd)
  MeanCp <- merge(MeanCp,StDevCp,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'StDevCp'
  MeanCp <- merge(MeanCp,TechRepNum,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'n'] <- 'TechRepNum'
  rm(StDevCp)
  rm(TechRepNum)
  rm(df)
  MeanCp$StDevCp[is.na(MeanCp$StDevCp)] <- 0
  ## Select data with StDev<1 and Cp less than 40
  Data <- subset(MeanCp, select = -StDevCp)
  Data$qPCRdate <- FolderFiles[[f]]
  Files[[f]] <- Data
  rm(MeanCp)
  rm(Data)
} 
RawData <- bind_rows(Files)
## Get list of primers in dataset
PrimerList <- unique(RawData$Primers)
Genotypes <- unique(as.character(RawData$Genotype))
ExperimentPrimers <- setdiff(PrimerList,ControlPrimers)
## Create date and genotype columns from Name
RawData$NameTwo <- RawData$Name
## I name my samples in format DDMMYY_Genotype but others might not, the if statements should detect this accordingly.
if (grepl("_",RawData$NameTwo)) {
  RawData <- separate(RawData, NameTwo, into = c("NameThree", "NameFour") ,sep = "_")
  if (!is.na(as.numeric(RawData$NameThree))){
    RawData$Date <- RawData$NameThree
    RawData$Genotype <- RawData$NameFour
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
   } else {
    RawData$Date <- RawData$NameFour
    RawData$Genotype <- RawData$NameThree
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
   }
  NameGenotypeDate <- subset(RawData, select = c(Name, Genotype, Date))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Name), ]
} else if (grepl(" ",RawData$NameTwo)) {
  RawData <- separate(RawData, NameTwo, into = c("NameThree", "NameFour") ,sep = " ")
  if (!is.na(as.numeric(RawData$NameThree))){
    RawData$Date <- RawData$NameThree
    RawData$Genotype <- RawData$NameFour
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
  } else {
    RawData$Date <- RawData$NameFour
    RawData$Genotype <- RawData$NameThree
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
  }
  NameGenotypeDate <- subset(RawData, select = c(Name, Genotype, Date))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Name), ]
} else {
  RawData$Genotype <- RawData$NameTwo
  RawData$NameTwo <- NULL
  NameGenotypeDate <- subset(RawData, select = c(Name, Genotype))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Name), ]
}
# this line determines the order factors are plotted on the x axis, any not mentioned are determined automatically
RawData <- RawData %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder)) 
## Make each primer its own column
DifferentPrimers <- list()
for (i in 1:as.numeric(length(PrimerList))){
  df <- subset(RawData, RawData$Primers == PrimerList[i])
  DifferentPrimers[[i]] <- df
  rm(df)
}
for (i in 1:as.numeric(length(DifferentPrimers))){
  colnames(DifferentPrimers[[i]])[colnames(DifferentPrimers[[i]]) == "Cp"] <- unique(DifferentPrimers[[i]]$Primers)
  DifferentPrimers[[i]] <- subset(DifferentPrimers[[i]], select = -c(Primers,Date,Genotype,TechRepNum))
}
SeparatedPrimers <- DifferentPrimers[[1]]
for (i in 2:as.numeric(length(DifferentPrimers))){
  SeparatedPrimers <- merge(SeparatedPrimers,DifferentPrimers[[i]], by= c("Name","qPCRdate"), all = TRUE)
}
rm(DifferentPrimers)
RawDataSeparatedPrimers <- merge(SeparatedPrimers,NameGenotypeDate,by="Name")
rm(SeparatedPrimers)

## Calculate mean of samples across different qPCR runs, getting rid of any outliers
StDevCp <- aggregate(RawData$Cp,by=list(Name=RawData$Name,Primers=RawData$Primers),FUN = sd)
StDevCp$x[is.na(StDevCp$x)] <- 0
df <- merge(StDevCp,RawData,by= c("Name","Primers"), all = TRUE)
Sb1 <- subset(df,x<1)
Sb2 <- subset(df,x>=1)
Sb2 <- subset(Sb2, TechRepNum == 2)
Sb2$x <- NULL
StDevSb2 <- aggregate(Sb2$Cp,by=list(Name=Sb2$Name,Primers=Sb2$Primers),FUN = sd)
StDevSb2$x[is.na(StDevSb2$x)] <- 0
Sb2 <- merge(StDevSb2,Sb2,by= c("Name","Primers"), all = TRUE)
rm(StDevSb2,df,StDevCp)
Sb2 <- dplyr::filter(Sb2, x < 1)
UsefulData <- bind_rows(Sb1,Sb2)
rm(Sb1,Sb2)
UsefulData <- subset(UsefulData, select = -c(x,TechRepNum,qPCRdate))
UsefulData <- aggregate(UsefulData$Cp,by=list(Name=UsefulData$Name,Primers=UsefulData$Primers),data=UsefulData,FUN = mean)
colnames(UsefulData)[colnames(UsefulData) == 'x'] <- 'MeanCp'
UsefulData <- merge(UsefulData,NameGenotypeDate,by="Name")

## Make each primer its own column
DifferentPrimers <- list()
for (i in 1:as.numeric(length(PrimerList))){
  df <- subset(UsefulData, UsefulData$Primers == PrimerList[i])
  DifferentPrimers[[i]] <- df
  rm(df)
}
for (i in 1:as.numeric(length(DifferentPrimers))){
  colnames(DifferentPrimers[[i]])[colnames(DifferentPrimers[[i]]) == "MeanCp"] <- unique(DifferentPrimers[[i]]$Primers)
  DifferentPrimers[[i]] <- subset(DifferentPrimers[[i]], select = -c(Primers,Date,Genotype))
}
SeparatedPrimers <- DifferentPrimers[[1]]
for (i in 2:as.numeric(length(DifferentPrimers))){
  SeparatedPrimers <- merge(SeparatedPrimers,DifferentPrimers[[i]], by= "Name", all = TRUE)
}
rm(DifferentPrimers)
UsefulDataSeparatedPrimers <- merge(SeparatedPrimers,NameGenotypeDate,by="Name")
rm(SeparatedPrimers)

## Convert Cp values to relative mRNA values
Data <- UsefulData
Data$RelativemRNA <- 2^-Data$MeanCp
Data <- subset(Data, select = -MeanCp)
DifferentPrimers <- list()
for (i in 1:as.numeric(length(PrimerList))){
  df <- subset(Data, Data$Primers == PrimerList[i])
  DifferentPrimers[[i]] <- df
  rm(df)
}
for (i in 1:as.numeric(length(DifferentPrimers))){
  colnames(DifferentPrimers[[i]])[colnames(DifferentPrimers[[i]]) == "RelativemRNA"] <- unique(DifferentPrimers[[i]]$Primers)
  DifferentPrimers[[i]] <- subset(DifferentPrimers[[i]], select = -c(Primers,Date,Genotype))
}
SeparatedPrimers <- DifferentPrimers[[1]]
for (i in 2:as.numeric(length(DifferentPrimers))){
  SeparatedPrimers <- merge(SeparatedPrimers,DifferentPrimers[[i]], by= "Name", all = TRUE)
}
rm(DifferentPrimers)
DataSeparatedPrimers <- merge(SeparatedPrimers,NameGenotypeDate,by="Name")
rm(SeparatedPrimers)

## Normalise data with control target genes
# Create two dataframes named after the control primers
NormalisedData <- list()
for (i in 1:as.numeric(length(ControlPrimers))){
  df <- subset(DataSeparatedPrimers, select = c(Name, Date, Genotype))
  NormalisedData[[i]] <- df
  assign(ControlPrimers[i],NormalisedData[[i]])
  rm(df)
}
rm(NormalisedData)
# To each control primer dataframe add the experiment primer mRNA values and normalise by dividing with the control primer mRNA value
for (i in 1:as.numeric(length(ControlPrimers))){
  df <- get(ControlPrimers[i])
  ControlValues <- DataSeparatedPrimers[[ControlPrimers[i]]]
  for (k in 1:as.numeric(length(ExperimentPrimers))){
    df[[ExperimentPrimers[k]]] <- DataSeparatedPrimers[[ExperimentPrimers[k]]]/ControlValues
  }
  df <- df %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder)) # this line determines the order factors are plotted on the x axis, any not mentioned are determined automatically
  assign(ControlPrimers[i],df)
  rm(df)
  rm(ControlValues)
}
rm(NameGenotypeDate)
#### end ####

## Limit data to RNA extracted on certain days (delete depending on experiment)
#RP49 <- dplyr::filter(RP49, Date %in% c("050321","090321","110321","150321","040621","270421"))

#### DeltaDelta Ct Analysis ####
# Analysis of data done the popular method published by Livak and Schmittgen in 2001
#### end ####

theme_jmt <- function(){theme(
  text = element_text(family = "sans", size = 15), # set default font to be Arial
  plot.title = element_text(hjust=0.5, face = 'bold'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.line = element_line(colour = "black"), # add axis lines in black
  aspect.ratio = 2/3, # set the relative length of the y/x axis
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  # tilt the x axis
}

#### Plot out individual repeats for each genotype ####
# Shows different genotypes individual biological repeat values for the target genes
ggplot(RawData, aes(x=Genotype, y=Cp, fill=Genotype)) + theme_jmt() +
  # Plots the boxplot of the data and specified axes, legend removed as it just shows the colours.
  geom_boxplot(show.legend = FALSE) + 
  # Adds a title and axis labels
  labs(title="All Cp values",x="Genotype", y = "Cp value") + 
  geom_jitter()
# Plot different dates for each geneotype
Plots <- list()
for (g in 1:as.numeric(length(Genotypes))){
  df <- subset(RawData, as.character(RawData$Genotype) == Genotypes[g])
  print(ggplot(df, aes(x=Primers, y=Cp, color = Primers)) + theme_jmt() +
    geom_point(stat = "identity")+
    geom_text(aes(label = Date), vjust = 0.5, hjust = 0, angle = 45))
  rm(df)
}
do.call(grid.arrange,Plots)

for (g in 1:as.numeric(length(Genotypes))){
  df <- subset(RawData, as.character(RawData$Genotype) == Genotypes[g])
  df <- df %>% gather("Primers", "Date", -Cp)
  print(ggplot(df, aes(x=Primers, y=Cp, color = Primers)) + theme_jmt() +
          geom_bar(group = "Date", stat = "identity"))
          rm(df)
}

ggplot(MeanCp, aes(fill=Genotype, y=MeanCp, x=Date)) + theme_jmt() +
  geom_bar(position = "dodge", stat = "identity")
ggplot(RP49, aes(fill=Genotype, y=mBeta, x=Date)) + theme_jmt() +
  geom_bar(position = "dodge", stat = "identity")
ggplot(RP49, aes(fill=Genotype, y=m3, x=Date)) + theme_jmt() +
  geom_bar(position = "dodge", stat = "identity")

#### end ####

#### Plot analysed data ####
# Plots the normalised data (relative mRNA amount of experiment genes divided by control genes)
RP49mBeta <- ggplot(RP49, aes(x=Genotype, y=mBeta, fill=Genotype)) + theme_jmt() +
  # Plots the boxplot of the data and specified axes, legend removed as it just shows the colours.
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") + 
  # Adds a title and axis labels
  labs(title="Normalised to RP49", y = "Relative mBeta expression") + 
  geom_jitter(show.legend = FALSE, size = 0.5) 

RP49m3 <- ggplot(RP49, aes(x=Genotype, y=m3, fill=Genotype)) + theme_jmt() +
  # Plots the boxplot of the data and specified axes, legend removed as it just shows the colours.
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  # Adds a title and axis labels
  labs(title="Normalised to RP49", y = "Relatve m3 expression") + 
  geom_jitter(show.legend = FALSE, size = 0.5) 
plot_grid(RP49mBeta, RP49m3, labels = "AUTO")
#### end ####


## Limit data to RNA extracted on certain days (delete depending on experiment)
#RP49 <- dplyr::filter(RP49, Date %in% c("050321","090321","110321","150321"))


#### Plot as bar graph ####

# need table of mean for each genotype with sem column
FinalGenotypes <- unique(RP49$Genotype)

df <- subset(RP49, RP49$Genotype == FinalGenotypes[1])
RP49Mean <- as.data.frame(mean(df[["mBeta"]], na.rm = TRUE))
RP49Mean$mBeta <- RP49Mean$`mean(df[["mBeta"]], na.rm = TRUE)`
RP49Mean$`mean(df[["mBeta"]], na.rm = TRUE)` <- NULL
RP49Mean$m3 <- mean(df[["m3"]], na.rm = TRUE)
RP49Mean$mBetaSEM <- sd(df$mBeta, na.rm=TRUE) / sqrt(length(df$mBeta[!is.na(df$mBeta)]))
RP49Mean$m3SEM <- sd(df$m3, na.rm=TRUE) / sqrt(length(df$m3[!is.na(df$m3)]))
RP49Mean$Genotype <- FinalGenotypes[1]
rm(df)

for (i in 2:as.numeric(length(FinalGenotypes))){
  df <- subset(RP49, RP49$Genotype == FinalGenotypes[i])
  Means <- as.data.frame(mean(df[["mBeta"]], na.rm = TRUE))
  Means$mBeta <- Means$`mean(df[["mBeta"]], na.rm = TRUE)`
  Means$`mean(df[["mBeta"]], na.rm = TRUE)` <- NULL
  Means$m3 <- mean(df[["m3"]], na.rm = TRUE)
  Means$mBetaSEM <- sd(df$mBeta, na.rm=TRUE) / sqrt(length(df$mBeta[!is.na(df$mBeta)]))
  Means$m3SEM <- sd(df$m3, na.rm=TRUE) / sqrt(length(df$m3[!is.na(df$m3)]))
  Means$Genotype <- FinalGenotypes[i]
  RP49Mean <- rbind(RP49Mean, Means)
  rm(df, Means)
}

df <- RP49Mean


RP49mBeta <- ggplot(df, aes(x=Genotype, y=mBeta, fill = Genotype, ymin=mBeta-mBetaSEM, ymax=mBeta+mBetaSEM)) + theme_jmt() +
  # Plots the boxplot of the data and specified axes, legend removed as it just shows the colours.
  geom_bar(stat='identity',show.legend = FALSE) + geom_errorbar(width = 0.4) 
  scale_fill_brewer(palette = "Set2") + 
  # Adds a title and axis labels
  labs(title="Expression of mBeta", y = "mRNA levels normalised to RPL32") 

RP49m3 <- ggplot(df, aes(x=Genotype, y=m3, fill = Genotype, ymin=m3-m3SEM, ymax=m3+m3SEM)) + theme_jmt() +
  # Plots the boxplot of the data and specified axes, legend removed as it just shows the colours.
  geom_bar(stat='identity',show.legend = FALSE) + geom_errorbar(width = 0.4) +
  scale_fill_brewer(palette = "Set2") +
  # Adds a title and axis labels
  labs(title="Expression of m3", y = "mRNA levels normalised to RPL32") 

plot_grid(RP49mBeta, RP49m3, labels = "AUTO")
#### end ####

#### Bar graphs with/out statistics and saved ####
## Remove an outlier from dataset
df <- subset(RP49, Name != c("040821_NoCRY 24hr Light","161021_Opto 2hr Light"))
my_m3_title <- expression(paste("Relative ", italic("m3"), " mRNA expression"))
my_mB_title <- expression(paste("Relative ", italic("mβ"), " mRNA expression"))

ggplot(df, aes(x=Genotype, y=m3, fill = Genotype)) + theme_jmt() + theme(axis.title.x=element_blank()) + 
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_m3_title)
ggsave("m3 expression of optogenetics.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

ggplot(df, aes(x=Genotype, y=mBeta, fill = Genotype)) + theme_jmt() + theme(axis.title.x=element_blank()) + 
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 24hr Light")), y_position = 0.17,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_mB_title)
ggsave("mBeta expression of optogenetics.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

ggplot(df, aes(x=Genotype, y=m3, fill = Genotype)) + theme_jmt() + theme(axis.title.x=element_blank()) + 
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("NoCRY Dark","NoCRY 24hr Light")), y_position = 0.03,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 1hr Light")), y_position = 0.034,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 2hr Light")), y_position = 0.038,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 8hr Light")), y_position = 0.042,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 24hr Light")), y_position = 0.046,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_m3_title)
ggsave("m3 expression of optogenetics PLUS STAT.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

ggplot(df, aes(x=Genotype, y=mBeta, fill = Genotype)) + theme_jmt() + theme(axis.title.x=element_blank()) + 
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("NoCRY Dark","NoCRY 24hr Light")), y_position = 0.06,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 1hr Light")), y_position = 0.04,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 2hr Light")), y_position = 0.06,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 8hr Light")), y_position = 0.10,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Opto Dark","Opto 24hr Light")), y_position = 0.17,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_mB_title)
ggsave("mBeta expression of optogenetics PLUS STAT.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")
#### end ####