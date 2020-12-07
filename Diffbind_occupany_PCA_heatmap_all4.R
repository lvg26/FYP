####Install packages####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
a

BiocManager::install("DiffBind")
a
BiocManager::install("org.Hs.eg.db")
a
BiocManager::install("ChIPpeakAnno")
a

###Load packages###
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 

###Data preep###
SampleList <- read.csv("E:/Project/Bioinformatics_R_no_LT_51/Spreadsheets/SampleList_Larapath_no_LT_51.csv", header=T) #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Project/Bioinformatics_R_no_LT_51/Spreadsheets/SampleList_Larapath_no_LT_51.csv")#create dba (an object within which data can be stored)#
TCLTNCNL.dba <- dba(myDBA)

###Find consensus peaks###
myDBA_consensus_TCLTNCNL <- dba.peakset(TCLTNCNL.dba, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_TCLTNCNL <- dba(myDBA_consensus_TCLTNCNL, mask=myDBA_consensus_TCLTNCNL$masks$Consensus, minOverlap=1) #Filter the peakset for consensus peaks. Because we want to find the number of peaks in each group separately, minOverlap = 1#
consensus_peaks_TCLTNCNL <- dba.peakset(myDBA_consensus_TCLTNCNL, bRetrieve=TRUE)  #save consensus peaks as a peak set to be used in the differential analysis#

###Plot venn - occupancy analysis###
dba.plotVenn(myDBA_consensus_TCLTNCNL, myDBA_consensus_TCLTNCNL$masks$Consensus, main = "Occupancy analysis all 4 phenotypes", label1 = "LT", 
             label2 = "NC", label3 = "NL", label4 = "TC")

###Plot heatmaps and PCA###
TCLTNCNL <- dba.count(TCLTNCNL.dba, peaks=consensus_peaks_TCLTNCNL, summits=250)#counts the number of reads in all samples corresponding to consensus peaks#
TCLTNCNL <- dba.normalize(TCLTNCNL) #normalise data by lib method - normalise the number of reads of each sample to the mean number of reads#

hmap <- colorRampPalette(c("white", "gray 64", "black"))(n=13) #make colour scheme#
dba.plotHeatmap(TCLTNCNL, correlations = TRUE, colScheme = hmap, distMethod = "pearson")#plot binding affinity heatmap - shows 5hmC enrichment (number of reads) for each consensus peak found in each sample#
dba.plotHeatmap(TCLTNCNL, correlations = FALSE, colScheme = hmap) #plot correlation heatmap - shows extent of correlation between samples#
dba.plotPCA(TCLTNCNL, label=1, labelSize = 0.6) #plot Principal Component Analysis (PCA) - shows sample relatedness#






