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
setwd("/Users/jonathantownson/Documents/PhD/ATAC data/Dimerisation mutant experiments (also includes some PEST)")
# Instruct the script of any Genotype ordering required e.g. control then experimental condition
XAxisOrder <- c("lacZ", "Wild Type")

#### Read in the data ####
RawData <- read.csv("/Users/jonathantownson/Documents/PhD/ATAC data/Dimerisation mutant experiments (also includes some PEST)/ATAC qPCR selected data for analysis.csv")

## Create date and genotype columns from Name
RawData$NameTwo <- RawData$Sample
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
  NameGenotypeDate <- subset(RawData, select = c(Sample, Genotype, Date))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Sample), ]
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
  NameGenotypeDate <- subset(RawData, select = c(Sample, Genotype, Date))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Sample), ]
} else {
  RawData$Genotype <- RawData$NameTwo
  RawData$NameTwo <- NULL
  NameGenotypeDate <- subset(RawData, select = c(Sample, Genotype))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Sample), ]
}

## Change the Genotype to the Δ versions
RawData <- RawData %>% 
  mutate(Genotype = ifelse(as.character(Genotype) == "NoNotch", "lacZ", as.character(Genotype)))
RawData <- RawData %>% 
  mutate(Genotype = ifelse(as.character(Genotype) == "wtNDECD", "Wild Type", as.character(Genotype)))
RawData <- RawData %>% 
  mutate(Genotype = ifelse(as.character(Genotype) == "R2008A", "R2007A", as.character(Genotype)))

# this line determines the order factors are plotted on the x axis, any not mentioned are determined automatically
RawData <- RawData %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder)) 

m3enhancer <- filter(RawData, Primer == "m3 enhancer")
mBetaenhancer <- filter(RawData, Primer == "mBeta enhancer")
#### end ####

theme_jmt <- function(){theme(
  text = element_text(family = "sans", size=15), # set default font to be Arial
  plot.title = element_text(hjust=0.5, face = 'bold'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.line = element_line(colour = "black"), # add axis lines in black
  aspect.ratio = 2/3, # set the relative length of the y/x axis
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  # tilt the x axis
}

#### Saved Plots ####

## Plot as bar graph
my_m3_title <- expression(paste("Enrichment of ", italic("m3"), " enhancer"))
my_mB_title <- expression(paste("Enrichment of ", italic("mβ"), " enhancer"))

ggplot(m3enhancer, aes(x=Genotype, y=Closed.control.normalised, fill = Genotype)) + theme_jmt() + 
  theme(axis.title.x=element_blank()) +
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("lacZ","Wild Type")), y_position = 90,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Wild Type","R2007A")), y_position = 95,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_m3_title)
ggsave("m3enh closed control normalised PLUS STAT.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

ggplot(mBetaenhancer, aes(x=Genotype, y=Closed.control.normalised, fill = Genotype)) + theme_jmt() + 
  theme(axis.title.x=element_blank()) +
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("lacZ","Wild Type")), y_position = 60,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  geom_signif(comparisons = list(c("Wild Type","R2007A")), y_position = 65,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_mB_title)
ggsave("mBetaenh closed control normalised PLUS STAT.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

ggplot(m3enhancer, aes(x=Genotype, y=Closed.control.normalised, fill = Genotype)) + theme_jmt() + 
  theme(axis.title.x=element_blank()) +
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_m3_title)
ggsave("m3enh closed control normalised.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

ggplot(mBetaenhancer, aes(x=Genotype, y=Closed.control.normalised, fill = Genotype)) + theme_jmt() + 
  theme(axis.title.x=element_blank()) +
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = my_mB_title)
ggsave("mBetaenh closed control normalised.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 12, units = "cm")

#### end ####