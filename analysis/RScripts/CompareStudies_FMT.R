#Author: Alicia Sorgen
#Date: 09-15-2021
#Description: p-value versus p-value plots for 16S datasets

R <- sessionInfo()
message(R$R.version$version.string)

#Libraries
library(ggplot2); message("ggplot2:", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra:", packageVersion("gridExtra"))
library(ggsignif); message("ggsignif:", packageVersion("ggsignif"))
library(ggrepel); message("ggrepel:", packageVersion("ggrepel"))
library(stringr); message("stringr:", packageVersion("stringr"))

rm(list=ls())

#Setup if running in RStudio
date <- Sys.Date()
date <- format(date, "%Y%b%d")
# date = "2021Oct06"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = paste0("~/BioLockJ_pipelines/FMTmouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("CompareStudies"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
   message("************* Running in BioLockJ *************")
   # studies = commandArgs(trailingOnly = TRUE)
} else {
   setwd(paste0(root, "/script"))
   # studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic", "NNK-FMT", "Smoke-FMT")
}

##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CombineStudies"),"/output/")
moduleDir = dirname(getwd())
outputDir = file.path(moduleDir,"output")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

#Compare study terms from each antibiotic experiment
levels = c("Phylum", "Class", "Order", "Family", "Genus")

for (level in levels) {
   
   # level = "Phylum"
   levelDir <- paste0(outputDir, "/", level, "/")
   dir.create(levelDir, showWarnings = FALSE)
   
   studies <- c("Smoke-FMT_Pre-tumor", "Smoke-FMT_Post-tumor", "NNK-FMT_Pre-tumor", "NNK-FMT_Post-tumor")
   
   studyNames <- c("Pre-tumor: Smoke", "Post-tumor: Smoke", "Pre-tumor: NNK", "Post-tumor: NNK")
   
   r <- vector()
   pval <- vector()
   plotList <- list()
   studyPairs <- vector()
   comparisonStudy <- vector()
   comparisonCancer <- vector()
   comparisonTumorStatus <- vector()
   pValueComparisons <- data.frame()
   indexNum <- vector()
   index <- 1
   
   for (s1 in 1:length(studies)) {
      
      if (s1!=length(studies)) {
         
         otherStudies <- c((s1 + 1):length(studies))
         
         for (s2 in otherStudies) {
            
            studyComparison <- paste0(studies[s1], " vs ", studies[s2])
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Term <- strsplit(studies[s1],"_")[[1]][2]
            path1 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LM_NicotineExposureComparisons_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Term <- strsplit(studies[s2],"_")[[1]][2]
            path2 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LM_NicotineExposureComparisons_ordered"), full.names = TRUE)
            
            comparison <- "Nic"
            TermType <- "TumorStatus"
            
            df <- compareStudiesFMT(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term)
            
            r[index]<-correlationBetweenStudies(df)[[1]]
            pval[index]<-correlationBetweenStudies(df)[[2]]
            index<-index+1
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   pval<-p.adjust(pval,method = "BH")
   count=0
   
   for (s1 in 1:length(studies)){
      
      if (s1!=length(studies)){
         otherStudies<-c((s1+1):length(studies))
         
         for (s2 in otherStudies){
            
            count<-count+1
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Term <- strsplit(studies[s1],"_")[[1]][2]
            path1 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LM_NicotineExposureComparisons_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Term <- strsplit(studies[s2],"_")[[1]][2]
            path2 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LM_NicotineExposureComparisons_ordered"), full.names = TRUE)
            
            comparison <- "Nic"
            TermType <- "TumorStatus"
            
            df <- compareStudiesFMT(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term)
            
            xlab=paste0(studyNames[s1]," vs. Control Recipient")
            ylab=paste0(studyNames[s2]," vs. Control Recipient")
            
            arrowRight <- "increased in control"
            arrowLeft <- "decreased in control"
            
            plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count], arrowRight, arrowLeft)
            
            plotList[[count]]<-plot
            
            studyPairs[count]<-paste0(studies[s1],"_v_",studies[s2])
            
            if(study1Name == study2Name) {
               comparisonStudy[[count]]<-"Same study"
            } else {
               comparisonStudy[[count]]<-"Different study"
            }
            

            if(study1Term == study2Term) {
               comparisonTumorStatus[[count]]<-"Same tumor status"
            } else {
               comparisonTumorStatus[[count]]<-"Different tumor status"
            }
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   df<-data.frame(studyPairs, comparisonStudy, comparisonTumorStatus, pval, r)
   pairOutFile <- file.path(levelDir, paste0(level,"_NicotineExposure_PairwiseComparison.tsv"))
   write.table(df, pairOutFile, sep="\t", row.names = FALSE)
   
   plotNum <- length(plotList)
   scatPlotOut <- file.path(levelDir, paste0(level,"_NicotineExposure_scatterPlots.pdf"))
   pdf(scatPlotOut, width = 10, height = 5)
   theme_set(theme_classic(base_size = 9))
   
   for (i in seq(from=1, to=plotNum, by=2)) {
      grid.arrange(plotList[[i]],plotList[[i+1]],
                   ncol=2,nrow=1)
      # message(i)
   }
   dev.off()
   
   
} # for (level in levels)






for (level in levels) {
   
   # level = "Phylum"
   levelDir <- paste0(outputDir, "/", level, "/")
   dir.create(levelDir, showWarnings = FALSE)
   
   studies <- c("Smoke-FMT_Smoke-exposed", "Smoke-FMT_Smoke-free", "NNK-FMT_NNK-exposed", "NNK-FMT_NNK-free")
   
   studyNames <- c("Smoke-exposed recipient: Pre-tumor", "Smoke-free recipient: Pre-tumor", "NNK-exposed recipient: Pre-tumor", "NNK-free recipient: Pre-tumor")
   
   r <- vector()
   pval <- vector()
   plotList <- list()
   studyPairs <- vector()
   comparisonStudy <- vector()
   comparisonCancer <- vector()
   comparisonNicotineExposure <- vector()
   pValueComparisons <- data.frame()
   indexNum <- vector()
   index <- 1
   
   for (s1 in 1:length(studies)) {
      
      if (s1!=length(studies)) {
         
         otherStudies <- c((s1 + 1):length(studies))
         
         for (s2 in otherStudies) {
            
            studyComparison <- paste0(studies[s1], " vs ", studies[s2])
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Term <- strsplit(studies[s1],"_")[[1]][2]
            path1 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LM_TumorStatusComparisons_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Term <- strsplit(studies[s2],"_")[[1]][2]
            path2 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LM_TumorStatusComparisons_ordered"), full.names = TRUE)
            
            comparison <- "Tumor"
            TermType <- "NicotineExposure"
            
            df <- compareStudiesFMT(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term)
            
            r[index]<-correlationBetweenStudies(df)[[1]]
            pval[index]<-correlationBetweenStudies(df)[[2]]
            index<-index+1
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   pval<-p.adjust(pval,method = "BH")
   count=0
   
   for (s1 in 1:length(studies)){
      
      if (s1!=length(studies)){
         otherStudies<-c((s1+1):length(studies))
         
         for (s2 in otherStudies){
            
            count<-count+1
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Term <- strsplit(studies[s1],"_")[[1]][2]
            path1 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LM_TumorStatusComparisons_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Term <- strsplit(studies[s2],"_")[[1]][2]
            path2 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LM_TumorStatusComparisons_ordered"), full.names = TRUE)
            
            comparison <- "Tumor"
            TermType <- "NicotineExposure"
            
            df <- compareStudiesFMT(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term)
            
            xlab=paste0(studyNames[s1]," vs. Post-tumor")
            ylab=paste0(studyNames[s2]," vs. Post-tumor")
            
            arrowRight <- "decreased post-tumor"
            arrowLeft <- "increased post-tumor"
            
            plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count], arrowRight, arrowLeft)
            
            plotList[[count]]<-plot
            
            studyPairs[count]<-paste0(studies[s1],"_v_",studies[s2])
            
            if(study1Name == study2Name) {
               comparisonStudy[[count]]<-"Same study"
            } else {
               comparisonStudy[[count]]<-"Different study"
            }
            
            
            if(study1Term == study2Term) {
               comparisonNicotineExposure[[count]]<-"Same nicotine exposure"
            } else {
               comparisonNicotineExposure[[count]]<-"Different nicotine exposure"
            }
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   df<-data.frame(studyPairs, comparisonStudy, comparisonNicotineExposure, pval, r)
   pairOutFile <- file.path(levelDir, paste0(level,"_TumorStatus_PairwiseComparison.tsv"))
   write.table(df, pairOutFile, sep="\t", row.names = FALSE)
   
   plotNum <- length(plotList)
   scatPlotOut <- file.path(levelDir, paste0(level,"_TumorStatus_scatterPlots.pdf"))
   pdf(scatPlotOut, width = 10, height = 5)
   theme_set(theme_classic(base_size = 9))
   
   for (i in seq(from=1, to=plotNum, by=2)) {
      grid.arrange(plotList[[i]],plotList[[i+1]],
                   ncol=2,nrow=1)
      # message(i)
   }
   dev.off()
   
   
} # for (level in levels)





for (level in levels) {
   
   # level = "Phylum"
   levelDir <- paste0(outputDir, "/", level, "/")
   dir.create(levelDir, showWarnings = FALSE)
   
   studies <- c("Smoke-FMT_Smoke-exposed", "Smoke-FMT_Smoke-free", "NNK-FMT_NNK-exposed")
   
   studyNames <- c("Smoke-exposed Donor", "Smoke-free Donor", "NNK-exposed Donor")
   
   r <- vector()
   pval <- vector()
   plotList <- list()
   studyPairs <- vector()
   comparisonStudy <- vector()
   comparisonNicotineExposure <- vector()
   pValueComparisons <- data.frame()
   indexNum <- vector()
   index <- 1
   
   for (s1 in 1:length(studies)) {
      
      if (s1!=length(studies)) {
         
         otherStudies <- c((s1 + 1):length(studies))
         
         for (s2 in otherStudies) {
            
            studyComparison <- paste0(studies[s1], " vs ", studies[s2])
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Term <- strsplit(studies[s1],"_")[[1]][2]
            path1 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LM_Donor_RecipientComparisons_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Term <- strsplit(studies[s2],"_")[[1]][2]
            path2 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LM_Donor_RecipientComparisons_ordered"), full.names = TRUE)
            
            comparison <- "FMTstatus"
            TermType <- "NicotineExposure"
            
            df <- compareStudiesFMT(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term)
            
            r[index]<-correlationBetweenStudies(df)[[1]]
            pval[index]<-correlationBetweenStudies(df)[[2]]
            index<-index+1
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   pval<-p.adjust(pval,method = "BH")
   count=0
   
   for (s1 in 1:length(studies)){
      
      if (s1!=length(studies)){
         otherStudies<-c((s1+1):length(studies))
         
         for (s2 in otherStudies){
            
            count<-count+1
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Term <- strsplit(studies[s1],"_")[[1]][2]
            path1 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LM_Donor_RecipientComparisons_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Term <- strsplit(studies[s2],"_")[[1]][2]
            path2 <- dir(pipeRoot, pattern = paste0("FMT-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LM_Donor_RecipientComparisons_ordered"), full.names = TRUE)
            
            comparison <- "FMTstatus"
            TermType <- "NicotineExposure"
            
            df <- compareStudiesFMT(level, file1, TermType, comparison, study1Name, study1Term, file2, study2Name, study2Term)
            
            xlab=paste0(studyNames[s1]," vs. Recipient")
            ylab=paste0(studyNames[s2]," vs. Recipient")
            
            arrowRight <- "increased in recipient"
            arrowLeft <- "decreased in recipient"
            
            plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count], arrowRight, arrowLeft)
            
            plotList[[count]]<-plot
            
            studyPairs[count]<-paste0(studies[s1],"_v_",studies[s2])
            
            if(study1Name == study2Name) {
               comparisonStudy[[count]]<-"Same study"
            } else {
               comparisonStudy[[count]]<-"Different study"
            }
            
            
            if(study1Term == study2Term) {
               comparisonNicotineExposure[[count]]<-"Same nicotine exposure"
            } else {
               comparisonNicotineExposure[[count]]<-"Different nicotine exposure"
            }
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   df<-data.frame(studyPairs, comparisonStudy, comparisonNicotineExposure, pval, r)
   pairOutFile <- file.path(levelDir, paste0(level,"_Donor_Recipient_PairwiseComparison.tsv"))
   write.table(df, pairOutFile, sep="\t", row.names = FALSE)
   
   plotNum <- length(plotList)
   scatPlotOut <- file.path(levelDir, paste0(level,"_Donor_Recipient_scatterPlots.pdf"))
   pdf(scatPlotOut, width = 15, height = 5)
   theme_set(theme_classic(base_size = 9))
   
   for (i in seq(from=1, to=plotNum, by=3)) {
      grid.arrange(plotList[[i]],plotList[[i+1]],plotList[[i+2]],
                   ncol=3,nrow=1)
      # message(i)
   }
   dev.off()
   
   
} # for (level in levels)
