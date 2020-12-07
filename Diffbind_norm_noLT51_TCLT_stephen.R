####Load packages####
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 
library(tidyverse)

###Data prep###
SampleList <- read.csv("E:/Project/Bioinformatics_R_no_LT_51/Spreadsheets/SampleList_Larapath_no_LT_51.csv", header=T)  #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Project/Bioinformatics_R_no_LT_51/Spreadsheets/SampleList_Larapath_no_LT_51.csv")#create dba (an object within which data can be stored)#
myDBA
TCLT <- dba(myDBA, mask=myDBA$masks$Tumour) #select for only liver samples#


###Find consensus peaks###
myDBA_consensus_TCLT <- dba.peakset(TCLT, consensus=c(DBA_CONDITION, DBA_TISSUE), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_TCLT <- dba(myDBA_consensus_TCLT, mask=myDBA_consensus_TCLT$masks$Consensus, minOverlap=1)#Filter the peakset for consensus peaks. Because we want to find the number of peaks in each group separately, minOverlap = 1#
consensus_peaks_TCLT <- dba.peakset(myDBA_consensus_TCLT, bRetrieve=TRUE) #save consensus peaks as a peak set to be used in the differential analysis#


###Differential binding analysis###
TCLT.count <- dba.count(TCLT, consensus_peaks_TCLT, summits = 250) #count number of reads corresponding to consensus peaks for each sample - 5hmC enrichment of consensus peaks. Summits makes the width of all the peaks the same (250 bp up and downstream of peak summit)#
TCLT_normed <- dba.normalize(TCLT.count) #normalise number of peaks based on lib method - normalises by the mean number of reads across all samples being analysed#
TCLT_contrast <- dba.contrast(TCLT_normed, categories=DBA_TISSUE) #Tells diffbind which comparison we are interested in i.e. we want to compare the 2 conditions for Liver: normal and tumour#
TCLT_contrast #displays which group (group1 or group2) each sample is in the differential analysis - need to know for filtering (see below)#
TCLT_analyzed <- dba.analyze(TCLT_contrast, bBlacklist = FALSE, bGreylist = FALSE) #Run default differential binding analysis (statistical sinificance given by p-value) based on negative binomial distribution = DESeq2. Our data has already been filtered against a "greylist#
TCLT.peaks <- dba.report(TCLT_analyzed, th=1) #Shows the differential binding analysis results. Shows the fold change in peak size between group1 and group2. th=1 will return all peak results (not filtered by a p value)#
TCLT.peaks
data("TSS.human.GRCh38") #Load human genome dataset from org.Hs.eg.db package#
TCLT.peaks <- annotatePeakInBatch(TCLT.peaks, AnnotationData=TSS.human.GRCh38) # match the peaks found with genes in the human genome dataset TSS#
TCLT.peaks <- addGeneIDs(TCLT.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol')) #annotate the peaks found with gene names#

write.csv(TCLT.peaks, "E:/Project/Diffbind_R_norm_noLT_51/Spreadsheets/Stephen_method/TCLT\\TCLT.peaks.csv")#output results of differential analysis as csv file#

###Filter for siginificant fold change in 5hmC enrichment###
as.data.frame(TCLT.peaks) %>% #save results of differential analysis as a data frame (data object in R)#
  filter(FDR < 0.1 & Fold > 0) %>% #Filter for significantly greater fold enrichment in group1 (TC)/ significantly less fold enrichment in group2 (LT)#
  select(seqnames, start, end, symbol, feature, FDR) %>% #select variables to present in the csv file#
  write.csv("E:/Project/Diffbind_R_norm_noLT_51/Spreadsheets/Stephen_method\\TumourColon_FC_UP_FDR0.1.csv")#output as csv file#

as.data.frame(TCLT.peaks) %>% #save results of differential analysis as a data frame (data object in R)#
  filter(Fold < 0 & FDR < 0.1) %>% #Filter for significantly  less fold enrichment in group1 (TC)/ significantly greater fold enrichment in group2 (LT)#
  select(seqnames, start, end, symbol, feature, FDR) %>% #select variables to present in the csv file#
  write.csv("E:/Project/Diffbind_R_norm_noLT_51/Spreadsheets/Stephen_method\\LiverTumour_FC_UP.FDR0.1.csv")#output as csv file#


