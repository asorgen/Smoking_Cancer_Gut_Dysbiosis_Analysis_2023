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
root = dir(root, pattern=paste0("CompareStudiesPost"), full.names=TRUE)

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
   
   # studies <- c("Abx-Bladder_Post_Smoke", "Abx-Colon_Post_Smoke", "Abx-Pancreatic_Post_Smoke", "Neo-Pancreatic_Post_Smoke", "Abx-Pancreatic_Pre_Smoke")
   # studyNames <- c("Post-Bladder: Smoke only", "Post-Colon: Smoke only", "Post-Pancreatic: Smoke only", "Post-Pancreatic: Smoke only (neo)", "Pre-Pancreatic: Smoke only")

   studies <- c("Abx-Bladder_Post_Smoke", "Abx-Colon_Post_Smoke", "Abx-Pancreatic_Post_Smoke")
   studyNames <- c("Post-Bladder: Smoke only", "Post-Colon: Smoke only", "Post-Pancreatic: Smoke only")
   
   r <- vector()
   pval <- vector()
   plotList <- list()
   studyPairs <- vector()
   comparisonStudy <- vector()
   comparisonCancer <- vector()
   comparisonTreatment <- vector()
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
            study1Tumor <- strsplit(studies[s1],"_")[[1]][2]
            study1Term <- strsplit(studies[s1],"_")[[1]][3]
            path1 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Tumor <- strsplit(studies[s2],"_")[[1]][2]
            study2Term <- strsplit(studies[s2],"_")[[1]][3]
            path2 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            df <- compareStudies(level, file1, study1Name, study1Tumor, study1Term, file2, study2Name, study2Tumor, study2Term)
            
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
            study1Cancer <- strsplit(study1Name,"-")[[1]][2]
            study1Tumor <- strsplit(studies[s1],"_")[[1]][2]
            study1Term <- strsplit(studies[s1],"_")[[1]][3]
            path1 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Cancer <- strsplit(study2Name,"-")[[1]][2]
            study2Tumor <- strsplit(studies[s2],"_")[[1]][2]
            study2Term <- strsplit(studies[s2],"_")[[1]][3]
            path2 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            df <- compareStudies(level, file1, study1Name, study1Tumor, study1Term, file2, study2Name, study2Tumor, study2Term)
            
            xlab=paste0(studyNames[s1]," vs. Control")
            ylab=paste0(studyNames[s2]," vs. Control")
            arrowRight <- "increased in smoke"
            arrowLeft <- "decreased in smoke"
            
            plot<-plotPairwiseSegments(df,xlab,ylab,r[count],pval[count], arrowRight, arrowLeft)
            
            plotList[[count]]<-plot
            
            studyPairs[count]<-paste0(studies[s1],"_v_",studies[s2])
            
            if(study1Name == study2Name) {
               comparisonStudy[[count]]<-"Same study"
            } else {
               comparisonStudy[[count]]<-"Different study"
            }
            
            if(study1Cancer == study2Cancer) {
               comparisonCancer[[count]]<-"Same tumor type"
            } else {
               comparisonCancer[[count]]<-"Different tumor type"
            }
            
            
            if(study1Term == study2Term) {
               comparisonTreatment[[count]]<-"Same treatment"
            } else {
               comparisonTreatment[[count]]<-"Different treatment"
            }
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   df<-data.frame(studyPairs, comparisonStudy, comparisonCancer, comparisonTreatment, pval, r)
   pairOutFile <- file.path(levelDir, paste0(level,"_PairwiseComparison.tsv"))
   write.table(df, pairOutFile, sep="\t", row.names = FALSE)
   
   plotNum <- length(plotList)
   scatPlotOut <- file.path(levelDir, paste0(level,"_scatterPlots.pdf"))
   pdf(scatPlotOut, width = 15, height = 5)
   theme_set(theme_classic(base_size = 9))
   
   for (i in seq(from=1, to=plotNum, by=3)) {
      grid.arrange(plotList[[i]],plotList[[i+1]], plotList[[i+2]],
                   ncol=3,nrow=1)
      # message(i)
   }
   dev.off()
   
   
} # for (level in levels)


































##### This function generate p-value versus p-value plot #####
makePlot<-function(df, xlab, ylab, coeficient, p, arrowRight, arrowLeft){
   
   df1<-df[(df[,"pval1"]<log10(0.05) & df[,"pval2"]<log10(0.05)) | (df[,"pval1"]>-log10(0.05) & df[,"pval2"]>-log10(0.05)),]
   df2<-df[(df[,"pval1"]<log10(0.05) & df[,"pval2"]>-log10(0.05)) | (df[,"pval1"]>-log10(0.05) & df[,"pval2"]<log10(0.05)),]
   
   if (p==0){
      p.1 = "< 2.2e-16"
   } else {
      p.1 = paste("=",format(p,digits = 3))
   }
   
   theme_set(theme_classic())
   
   lim = ceiling(max(abs( c( df$pval1, df$pval2) )))
   inc <- lim / 2
   dec <- -(lim / 2)
   
   plot<-ggplot(data=df,aes(x=pval1,y=pval2))+
      geom_point(size=1, shape = 1)+
      geom_segment(aes(x = log10(0.05), y = -lim, xend = log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "red")+
      geom_segment(aes(x = -lim, y = log10(0.05), xend = log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "red")+
      geom_segment(aes(x = -log10(0.05), y = lim, xend = -log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "red")+
      geom_segment(aes(x = lim, y = -log10(0.05), xend = -log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "red")+
      
      geom_segment(aes(x = log10(0.05), y = lim, xend = log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "black")+
      geom_segment(aes(x = -lim, y = -log10(0.05), xend = log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "black")+
      geom_segment(aes(x = -log10(0.05), y = -lim, xend = -log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "black")+
      geom_segment(aes(x = -log10(0.05), y = log10(0.05), xend = lim, yend = log10(0.05)), linetype = "dashed", color = "black")+
      labs(x=xlab,y=ylab,title = paste0("Spearman coefficient = ",format(coeficient,digits =3),"\nAdjusted p ",p.1))+
      geom_text_repel(data=df1,aes(x=,y=pval2,label=bugName),segment.colour="red",size=2.5,min.segment.length = 0,
                      segment.color="grey",segment.size=0.2)+
      geom_text_repel(data=df2,aes(x=,y=pval2,label=bugName),segment.colour="red",size=2.5,min.segment.length = 0,
                      segment.color="grey",segment.size=0.2)+
      geom_point(data=df1, aes(x=pval1,y=pval2), color='blue', size=2)+
      annotate("segment", x=-log10(0.05), y=-(lim+1), xend=lim, yend=-(lim+1),
               col="black", arrow=arrow(length=unit(0.3, "cm"))) +
      annotate("text", x=inc, y=-(lim+1.5), label = arrowRight) +
      annotate("segment", x=log10(0.05), y=-(lim+1), xend=-lim, yend=-(lim+1),
               col="black", arrow=arrow(length=unit(0.3, "cm")))+
      annotate("text", x=dec, y=-(lim+1.5), label = arrowLeft)
   
   
}

for (level in levels) {
   
   # level = "Phylum"
   levelDir <- paste0(outputDir, "/", level, "/")
   dir.create(levelDir, showWarnings = FALSE)
   
   studies <- c("Abx-Bladder_Post_Smoke", "Abx-Colon_Post_Smoke", "Abx-Pancreatic_Post_Smoke")
   studyNames <- c("Post-Bladder: Smoke only", "Post-Colon: Smoke only", "Post-Pancreatic: Smoke only")
   
   r <- vector()
   pval <- vector()
   plotList <- list()
   indexNum <- vector()
   index <- 1
   
   for (s1 in 1:length(studies)) {
      
      if (s1!=length(studies)) {
         
         otherStudies <- c((s1 + 1):length(studies))
         
         for (s2 in otherStudies) {
            
            studyComparison <- paste0(studies[s1], " vs ", studies[s2])
            
            ## Name study 1 parameters
            study1Name <- strsplit(studies[s1],"_")[[1]][1]
            study1Tumor <- strsplit(studies[s1],"_")[[1]][2]
            study1Term <- strsplit(studies[s1],"_")[[1]][3]
            path1 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Tumor <- strsplit(studies[s2],"_")[[1]][2]
            study2Term <- strsplit(studies[s2],"_")[[1]][3]
            path2 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            df <- compareStudies(level, file1, study1Name, study1Tumor, study1Term, file2, study2Name, study2Tumor, study2Term)
            
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
            study1Cancer <- strsplit(study1Name,"-")[[1]][2]
            study1Tumor <- strsplit(studies[s1],"_")[[1]][2]
            study1Term <- strsplit(studies[s1],"_")[[1]][3]
            path1 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path1 <- paste0(path1, "/output/", level)
            file1 <- list.files(path = path1, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            ## Name study 2 parameters
            study2Name <- strsplit(studies[s2],"_")[[1]][1]
            study2Cancer <- strsplit(study2Name,"-")[[1]][2]
            study2Tumor <- strsplit(studies[s2],"_")[[1]][2]
            study2Term <- strsplit(studies[s2],"_")[[1]][3]
            path2 <- dir(pipeRoot, pattern = paste0("Abx-LinearModeling"), full.names = TRUE)
            path2 <- paste0(path2, "/output/", level)
            file2 <- list.files(path = path2, pattern = paste0(level, "_Combined_LMresults_ordered"), full.names = TRUE)
            
            df <- compareStudies(level, file1, study1Name, study1Tumor, study1Term, file2, study2Name, study2Tumor, study2Term)
            
            xlab=paste0(studyNames[s1]," vs. Control")
            ylab=paste0(studyNames[s2]," vs. Control")
            
            arrowRight <- "increased in smoke"
            arrowLeft <- "decreased in smoke"
            
            plot<-makePlot(df,xlab,ylab,r[count],pval[count], arrowRight, arrowLeft)
            
            plotList[[count]]<-plot
            
         } # for (s2 in otherStudies)
      } # if (s!=length(studies))
   } # for (s1 in 1:length(studies))
   
   plotNum <- length(plotList)
   scatPlotOut <- file.path(levelDir, paste0(level,"_scatterPlots_update.pdf"))
   pdf(scatPlotOut, width = 15, height = 5)
   theme_set(theme_classic(base_size = 9))
   
   for (i in seq(from=1, to=plotNum, by=3)) {
      grid.arrange(plotList[[i]],plotList[[i+1]], plotList[[i+2]],
                   ncol=3,nrow=1)
      # message(i)
   }
   dev.off()
   
   
} # for (level in levels)
