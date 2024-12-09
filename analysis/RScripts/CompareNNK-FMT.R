#Author: Alicia Sorgen
#Date: 03-22-2021
#Description: p-value versus p-value plots for 16S datasets

#Libraries
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)
library(stringr)

rm(list=ls())

date <- Sys.Date()
date <- format(date, "%Y%b%d")
date = "2021Sep02"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("CompareNNK-FMT"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  # studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  # studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic", "NNK-FMT", "Smoke-FMT")
}

##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/")
moduleDir = dirname(getwd())
outputDir = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

study <- "NNK-FMT"
inputStudy = paste0(inputDir, "/", study, "/")

levels = c("Phylum", "Class", "Order", "Family", "Genus")

for (level in levels) {
  
  levelDir <- paste0(outputDir, "/", level, "/")
  dir.create(levelDir, showWarnings = FALSE)
  
  studies <- c("Pre-tumor_NNKExp", "Post-tumor_NNKExp")
  studyNames <- c("NNK-FMT Pre-tumor: NNK-exposed Recipient", "NNK-FMT Post-tumor: NNK-exposed Recipient")
  
  r <- vector()
  pval <- vector()
  plotList <- list()
  studyPairs <- vector()
  comparison <- vector()
  comparisonPhenotype <- vector()
  pValueComparisons <- data.frame()
  indexNum <- vector()
  index <- 1
  
  for (s in 1:length(studies)) {
    
    if (s!=length(studies)) {
      
      otherStudies <- c((s + 1):length(studies))
      
      for (s1 in otherStudies) {
        
        studyComparison <- paste0(studies[s], " vs ", studies[s1])
        
        ## Name study 1 parameters
        study1Name <- study
        study1Variable <- strsplit(studies[s],"_")[[1]][1]
        study1Term <- strsplit(studies[s],"_")[[1]][2]
        path1 <- dir(pipeRoot, pattern = paste0(study1Name, "LM"), full.names = TRUE)
        path1 <- paste0(path1, "/output/", level, "/")
        file1 <- list.files(path = path1, pattern = paste0(study1Name, "_", study1Variable, "_", study1Term, "_LinearModelResults"), full.names = TRUE)
        
        ## Name study 2 parameters
        study2Name <- study
        study2Variable <- strsplit(studies[s1],"_")[[1]][1]
        study2Term <- strsplit(studies[s1],"_")[[1]][2]
        path2 <- dir(pipeRoot, pattern = paste0(study2Name, "LM"), full.names = TRUE)
        path2 <- paste0(path2, "/output/", level, "/")
        file2 <- list.files(path = path2, pattern = paste0(study2Name, "_", study2Variable, "_", study2Term, "_LinearModelResults"), full.names = TRUE)
        
        df <- compareStudies(level, file1, study1Term, file2, study2Term)
        
        r[index] <- correlationBetweenStudies(df)[[1]]
        pval[index] <- correlationBetweenStudies(df)[[2]]
        
        Spearman_pval <- correlationBetweenStudies(df)[[2]]
        R_squared <- correlationBetweenStudies(df)[[1]]
        pValues <- df
        pValues <- cbind(studyComparison, pValues)
        pValues<-pValues[(pValues[,"pval1"]<log10(0.05) & pValues[,"pval2"]<log10(0.05)) | (pValues[,"pval1"]>-log10(0.05) & pValues[,"pval2"]>-log10(0.05)),]
        pValues <- cbind(pValues, Spearman_pval)
        pValueComparisons <- cbind(pValues, R_squared)
        # pValueComparisons <- rbind(pValueComparisons, pValues)
        
        studyPairs[index]<-paste0(studies[s]," v ",studies[s1])
        
        indexNum[index] <- index
        index <- index + 1
        
      } # for (s1 in otherStudies)
    } # if (s!=length(studies))
  } # for (s in 1:length(studies))
  
  pValueComparisons$pAdj <- p.adjust(pValueComparisons$Spearman_pval, method = "BH")
  write.table(pValueComparisons, paste0(levelDir, level,"_", study, "_Pre_v_Post-tumor_by_", study1Term, "_SignificantTaxa.tsv"), sep="\t", row.names = FALSE)

  pvalAdj <- p.adjust(pval, method = "BH")
  dFrame<-data.frame(indexNum,studyPairs,pval,r, pvalAdj)
  pairOutFile <- paste0(levelDir, level,"_", study, "_Pre_v_Post-tumor_by_", study1Term,"_PairwiseComparison.tsv")
  write.table(dFrame, pairOutFile, sep="\t", row.names = FALSE)


  plotList <- list()
  r <- vector()
  pval <- vector()
  count = 0

  for (s in 1:length(studies)) {

    if (s!=length(studies)) {

      otherStudies <- c((s+1):length(studies))

      for (s1 in otherStudies) {
        count <- count +1

        ## Name study 1 parameters
        study1Name <- study
        study1Variable <- strsplit(studies[s],"_")[[1]][1]
        study1Term <- strsplit(studies[s],"_")[[1]][2]
        path1 <- dir(pipeRoot, pattern = paste0(study1Name, "LM"), full.names = TRUE)
        path1 <- paste0(path1, "/output/", level, "/")
        file1 <- list.files(path = path1, pattern = paste0(study1Name, "_", study1Variable, "_", study1Term, "_LinearModelResults"), full.names = TRUE)
        
        ## Name study 2 parameters
        study2Name <- study
        study2Variable <- strsplit(studies[s1],"_")[[1]][1]
        study2Term <- strsplit(studies[s1],"_")[[1]][2]
        path2 <- dir(pipeRoot, pattern = paste0(study2Name, "LM"), full.names = TRUE)
        path2 <- paste0(path2, "/output/", level, "/")
        file2 <- list.files(path = path2, pattern = paste0(study2Name, "_", study2Variable, "_", study2Term, "_LinearModelResults"), full.names = TRUE)
        
        df <- compareStudies(level, file1, study1Term, file2, study2Term)
          
        r[count] <- correlationBetweenStudies(df)[[1]]
        pval[count] <- correlationBetweenStudies(df)[[2]]

        xlab=paste0(studyNames[s]," vs. NNK-free Recipient")
        ylab=paste0(studyNames[s1]," vs. NNK-free Recipient")

        plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count])
          
        plotList[[count]] <- plot
          
      } # for (s1 in otherStudies)
        
    } # if (s!=length(studies))

  } # for (s in 1:length(studies))

  plotNum <- length(plotList)
    
  scatPlotOut <- paste0(levelDir, level,"_", study, "_Pre_v_Post-tumor_by_", study1Term, "_scatterPlotsSegments.pdf")
  pdf(scatPlotOut, width = 5, height = 5)
  theme_set(theme_classic(base_size = 9))
  for (i in plotNum) {
    print(plotList)
  } # for (i in plotNum)
  dev.off()

  
  
  
  studies <- c("NNK-exposed_TumorImp", "NNK-free_TumorImp")
  studyNames <- c("NNK-exposed FMT recipients: Pre-tumor", "NNK-free FMT recipients: Pre-tumor")
  
  r <- vector()
  pval <- vector()
  plotList <- list()
  studyPairs <- vector()
  comparison <- vector()
  comparisonPhenotype <- vector()
  pValueComparisons <- data.frame()
  indexNum <- vector()
  index <- 1
  
  for (s in 1:length(studies)) {
    
    if (s!=length(studies)) {
      
      otherStudies <- c((s + 1):length(studies))
      
      for (s1 in otherStudies) {
        
        studyComparison <- paste0(studies[s], " vs ", studies[s1])
        
        ## Name study 1 parameters
        study1Name <- study
        study1Variable <- strsplit(studies[s],"_")[[1]][1]
        study1Term <- strsplit(studies[s],"_")[[1]][2]
        path1 <- dir(pipeRoot, pattern = paste0(study1Name, "LM"), full.names = TRUE)
        path1 <- paste0(path1, "/output/", level, "/")
        file1 <- list.files(path = path1, pattern = paste0(study1Name, "_", study1Variable, "_", study1Term, "_LinearModelResults"), full.names = TRUE)
        
        ## Name study 2 parameters
        study2Name <- study
        study2Variable <- strsplit(studies[s1],"_")[[1]][1]
        study2Term <- strsplit(studies[s1],"_")[[1]][2]
        path2 <- dir(pipeRoot, pattern = paste0(study2Name, "LM"), full.names = TRUE)
        path2 <- paste0(path2, "/output/", level, "/")
        file2 <- list.files(path = path2, pattern = paste0(study2Name, "_", study2Variable, "_", study2Term, "_LinearModelResults"), full.names = TRUE)
        
        df <- compareStudies(level, file1, study1Term, file2, study2Term)
        
        r[index] <- correlationBetweenStudies(df)[[1]]
        pval[index] <- correlationBetweenStudies(df)[[2]]
        
        Spearman_pval <- correlationBetweenStudies(df)[[2]]
        R_squared <- correlationBetweenStudies(df)[[1]]
        pValues <- df
        pValues <- cbind(studyComparison, pValues)
        # pValues<-pValues[(pValues[,"pval1"]<log10(0.05) & pValues[,"pval2"]<log10(0.05)) | (pValues[,"pval1"]>-log10(0.05) & pValues[,"pval2"]>-log10(0.05)),]
        # pValues <- cbind(pValues, Spearman_pval)
        # pValueComparisons <- cbind(pValues, R_squared)

        studyPairs[index]<-paste0(studies[s]," v ",studies[s1])
        
        indexNum[index] <- index
        index <- index + 1
        
      } # for (s1 in otherStudies)
    } # if (s!=length(studies))
  } # for (s in 1:length(studies))
  
  # pValueComparisons$pAdj <- p.adjust(pValueComparisons$Spearman_pval, method = "BH")
  # write.table(pValueComparisons, paste0(levelDir, level,"_", study, "_Pre_v_Post_", study1Term, "_SignificantTaxa.tsv"), sep="\t", row.names = FALSE)
  
  pvalAdj <- p.adjust(pval, method = "BH")
  dFrame<-data.frame(indexNum,studyPairs,pval,r, pvalAdj)
  pairOutFile <- paste0(levelDir, level,"_", study, "_NNKExposure_by_", study1Term,"_PairwiseComparison.tsv")
  write.table(dFrame, pairOutFile, sep="\t", row.names = FALSE)
  
  
  plotList <- list()
  r <- vector()
  pval <- vector()
  count = 0
  
  for (s in 1:length(studies)) {
    
    if (s!=length(studies)) {
      
      otherStudies <- c((s+1):length(studies))
      
      for (s1 in otherStudies) {
        count <- count +1
        
        ## Name study 1 parameters
        study1Name <- study
        study1Variable <- strsplit(studies[s],"_")[[1]][1]
        study1Term <- strsplit(studies[s],"_")[[1]][2]
        path1 <- dir(pipeRoot, pattern = paste0(study1Name, "LM"), full.names = TRUE)
        path1 <- paste0(path1, "/output/", level, "/")
        file1 <- list.files(path = path1, pattern = paste0(study1Name, "_", study1Variable, "_", study1Term, "_LinearModelResults"), full.names = TRUE)
        
        ## Name study 2 parameters
        study2Name <- study
        study2Variable <- strsplit(studies[s1],"_")[[1]][1]
        study2Term <- strsplit(studies[s1],"_")[[1]][2]
        path2 <- dir(pipeRoot, pattern = paste0(study2Name, "LM"), full.names = TRUE)
        path2 <- paste0(path2, "/output/", level, "/")
        file2 <- list.files(path = path2, pattern = paste0(study2Name, "_", study2Variable, "_", study2Term, "_LinearModelResults"), full.names = TRUE)
        
        df <- compareStudies(level, file1, study1Term, file2, study2Term)
        
        r[count] <- correlationBetweenStudies(df)[[1]]
        pval[count] <- correlationBetweenStudies(df)[[2]]
        
        xlab=paste0(studyNames[s]," vs. Post-tumor")
        ylab=paste0(studyNames[s1]," vs. Post-tumor")
        
        plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count])
        
        plotList[[count]] <- plot
        
      } # for (s1 in otherStudies)
      
    } # if (s!=length(studies))
    
  } # for (s in 1:length(studies))
  
  plotNum <- length(plotList)
  
  scatPlotOut <- paste0(levelDir, level,"_", study, "NNKExposure_by_", study1Term, "_scatterPlotsSegments.pdf")
  pdf(scatPlotOut, width = 5, height = 5)
  theme_set(theme_classic(base_size = 9))
  for (i in plotNum) {
    print(plotList)
  } # for (i in plotNum)
  dev.off()
  
  
  
  
  
  
  studies.2 <- c("Pre-tumor__MouseType__Recipient_NNKExpYes", "Pre-tumor__MouseType__Recipient_NNKExpNo", "Post-tumor__MouseType__Recipient_NNKExpYes", "Post-tumor__MouseType__Recipient_NNKExpNo")
  studyNames.2 <- c("Pre-tumor: NNK-exposed Recipient", "Pre-tumor: NNK-free Recipient", "Post-tumor: NNK-exposed Recipient", "Post-tumor: NNK-free Recipient")
  
  r <- vector()
  pval <- vector()
  plotList <- list()
  studyPairs <- vector()
  # comparison <- vector()
  # comparisonPhenotype <- vector()
  pValueComparisons <- data.frame()
  indexNum <- vector()
  index <- 1
  
  for (s in 1:length(studies.2)) {
    
    if (s!=length(studies.2)) {
      
      otherStudies <- c((s + 1):length(studies.2))
      
      for (s1 in otherStudies) {
        
        studyComparison <- paste0(studies.2[s], " vs ", studies.2[s1])
        
        ## Name study 1 parameters
        study1Name <- study
        study1Variable <- strsplit(studies.2[s],"__")[[1]][1]
        study1Term <- strsplit(studies.2[s],"__")[[1]][2]
        study1Mouse <- strsplit(studies.2[s],"__")[[1]][3]
        path1 <- dir(pipeRoot, pattern = paste0(study1Name, "LM"), full.names = TRUE)
        path1 <- paste0(path1, "/output/", level, "/")
        file1 <- list.files(path = path1, pattern = paste0(study1Name, "_", study1Variable, "_", study1Term, "_LinearModelResults"), full.names = TRUE)
        
        ## Name study 2 parameters
        study2Name <- study
        study2Variable <- strsplit(studies.2[s1],"__")[[1]][1]
        study2Term <- strsplit(studies.2[s1],"__")[[1]][2]
        study2Mouse <- strsplit(studies.2[s1],"__")[[1]][3]
        path2 <- dir(pipeRoot, pattern = paste0(study2Name, "LM"), full.names = TRUE)
        path2 <- paste0(path2, "/output/", level, "/")
        file2 <- list.files(path = path2, pattern = paste0(study2Name, "_", study2Variable, "_", study2Term, "_LinearModelResults"), full.names = TRUE)
        
        df <- compareStudies(level, file1, study1Mouse, file2, study2Mouse)

        r[index] <- correlationBetweenStudies(df)[[1]]
        pval[index] <- correlationBetweenStudies(df)[[2]]

        Spearman_pval <- correlationBetweenStudies(df)[[2]]
        R_squared <- correlationBetweenStudies(df)[[1]]
        pValues <- df
        pValues <- cbind(studyComparison, pValues)
        pValues<-pValues[(pValues[,"pval1"]<log10(0.05) & pValues[,"pval2"]<log10(0.05)) | (pValues[,"pval1"]>-log10(0.05) & pValues[,"pval2"]>-log10(0.05)),]
        pValues <- cbind(pValues, Spearman_pval)
        pValueComparisons <- cbind(pValues, R_squared)

        studyPairs[index]<-paste0(studies.2[s]," v ",studies.2[s1])

        indexNum[index] <- index
        index <- index + 1
        
      } # for (s1 in otherStudies)
    } # if (s!=length(studies.2))
  } # for (s in 1:length(studies.2))
  
  pValueComparisons$pAdj <- p.adjust(pValueComparisons$Spearman_pval, method = "BH")
  write.table(pValueComparisons, paste0(levelDir, level,"_", study, "_Pre_v_Post-tumor_by_", study1Term, "_SignificantTaxa.tsv"), sep="\t", row.names = FALSE)
  
  pvalAdj <- p.adjust(pval, method = "BH")
  dFrame<-data.frame(indexNum,studyPairs,pval,r, pvalAdj)
  pairOutFile <- paste0(levelDir, level,"_", study, "_Pre_v_Post-tumor_by_", study1Term,"_PairwiseComparison.tsv")
  write.table(dFrame, pairOutFile, sep="\t", row.names = FALSE)
  
  
  plotList <- list()
  r <- vector()
  pval <- vector()
  count = 0
  
  for (s in 1:length(studies.2)) {
    
    if (s!=length(studies.2)) {
      
      otherStudies <- c((s+1):length(studies.2))
      
      for (s1 in otherStudies) {
        count <- count +1
        
        ## Name study 1 parameters
        study1Name <- study
        study1Variable <- strsplit(studies.2[s],"__")[[1]][1]
        study1Term <- strsplit(studies.2[s],"__")[[1]][2]
        study1Mouse <- strsplit(studies.2[s],"__")[[1]][3]
        path1 <- dir(pipeRoot, pattern = study1Name, full.names = TRUE)
        path1 <- paste0(path1, "/output/", level, "/")
        file1 <- list.files(path = path1, pattern = paste0(study1Name, "_", study1Variable, "_", study1Term, "_LinearModelResults"), full.names = TRUE)
        
        ## Name study 2 parameters
        study2Name <- study
        study2Variable <- strsplit(studies.2[s1],"__")[[1]][1]
        study2Term <- strsplit(studies.2[s1],"__")[[1]][2]
        study2Mouse <- strsplit(studies.2[s1],"__")[[1]][3]
        path2 <- dir(pipeRoot, pattern = study2Name, full.names = TRUE)
        path2 <- paste0(path2, "/output/", level, "/")
        file2 <- list.files(path = path2, pattern = paste0(study2Name, "_", study2Variable, "_", study2Term, "_LinearModelResults"), full.names = TRUE)
        
        df <- compareStudies(level, file1, study1Mouse, file2, study2Mouse)
        
        r[count] <- correlationBetweenStudies(df)[[1]]
        pval[count] <- correlationBetweenStudies(df)[[2]]
        
        xlab=paste0(studyNames.2[s]," vs. NNK-exposed Donor")
        ylab=paste0(studyNames.2[s1]," vs. NNK-exposed Donor")
        
        plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count])
        
        plotList[[count]] <- plot
        
      } # for (s1 in otherStudies)
      
    } # if (s!=length(studies.2))
    
  } # for (s in 1:length(studies.2))
  
  plotNum <- length(plotList)
  
  scatPlotOut <- paste0(levelDir, level, "_", study, "_Pre_v_Post-tumor_by_", study1Term, "_scatterPlotsSegments.pdf")
  pdf(scatPlotOut, width = 5, height = 5)
  theme_set(theme_classic(base_size = 9))
  for (i in plotNum) {
    print(plotList)
  } # for (i in plotNum)
  dev.off()
  
    
} # for (level in levels)

