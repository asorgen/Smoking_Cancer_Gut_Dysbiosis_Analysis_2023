#Author: Alicia Sorgen
#Date: 09-14-21
#Description: Performs basic linear modeling and produces simple boxplots

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))
library(tidyr); message("tidyr:", packageVersion("tidyr"))
library(ggplot2); message("ggplot2:", packageVersion("ggplot2"))

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
date = "2021Oct06"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
# root = dir(root, pattern=paste0("Abx-BladderLM"), full.names=TRUE)
root = dir(root, pattern=paste0("Abx-LinearModeling"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic")
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CombineStudies"),"/output/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

levels = c("Phylum", "Class", "Order", "Family", "Genus")

Level <- vector()
Experiment <- vector()
Sig_pVal <- vector()
smry <- 1

for (level in levels) {
  
  outputLevel = paste0(output, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  inputFile <- paste0(level, "_CombineStudies_ALL.tsv")
  fullPath <- paste0(inputDir, inputFile)
  
  myT <- read.table(file = fullPath, sep="\t", header = TRUE)
  myT <- myT[myT$AntibioticExp == "None",]
  
  BugName <- vector()
  Study <- vector()
  Tumor <- vector()
  p_Smoke <- vector()
  s_Smoke <- vector()
  index <- 1
  
  for (study in studies) {
    
    df <- myT[ myT$Study == study, ]
    df <- df[df$TumorImp == "Post",]
    t <- which( colnames(df) == level )
    
    for (bug in unique(df[,t])) {
      
      if ( !(bug == "Other") ) {
        
        df_taxa <- df[df[,t] == bug,]
        NormalizedCounts <- df_taxa$NormalizedCounts
        Treatment <- df_taxa$Treatment
        
        lm <- lm(NormalizedCounts ~ Treatment)
        sm <- summary(lm)
        
        BugName[index] <- bug
        Study[index] <- study
        Tumor[index] <- "Post" 
        p_Smoke[index] <- sm$coefficients[2,4]
        s_Smoke[index] <- sm$coefficients[2,1]
        index <- index + 1
        
      } # if ( !(bug == "Other") )
    } # for (bug in unique(df[,t]))
    
    if (study == "Abx-Pancreatic") {
      df <- myT[ myT$Study == study, ]
      df <- df[df$TumorImp == "Pre",]
      t <- which( colnames(df) == level )
      
      for (bug in unique(df[,t])) {
        
        if ( !(bug == "Other") ) {
          
          df_taxa <- df[df[,t] == bug,]
          NormalizedCounts <- df_taxa$NormalizedCounts
          Treatment <- df_taxa$Treatment
          
          lm <- lm(NormalizedCounts ~ Treatment)
          sm <- summary(lm)
          
          BugName[index] <- bug
          Study[index] <- study
          Tumor[index] <- "Pre" 
          p_Smoke[index] <- sm$coefficients[2,4]
          s_Smoke[index] <- sm$coefficients[2,1]
          index <- index + 1
          
        } # if ( !(bug == "Other") )
      } # for (bug in unique(df[,t]))
      
    }
  } # for (study in studies)
  
  dFrame <- data.frame(BugName, Study, Tumor, s_Smoke, p_Smoke)
  dFrame$pAdj_Smoke <- p.adjust(dFrame$p_Smoke, method = "BH")
  dFrame <- na.omit(dFrame)
  dFrame <- dFrame[order(dFrame$p_Smoke),]
  
  outputFile <- paste0(outputLevel, level, "_Combined_LMresults_ordered.tsv")
  write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
  
  dFrameFiltered <- dFrame[dFrame$pAdj_Smoke < 0.05,]
  
  for (study in studies) {
    dFrame3 <- dFrameFiltered[dFrameFiltered$Study == study,]
    Level[smry] <- level
    Experiment[smry] <- study
    Sig_pVal[smry] <- nrow(dFrame3)
    smry <- smry + 1
    
  }
  
  for (study in studies) {
    
    df <- myT[ myT$Study == study, ]
    dFrame_study <- dFrame[ dFrame$Study == study, ]
    
    model <- paste0("lm(NormalizedCounts ~ Treatment)")
    unadj <- dFrame_study$p_Smoke
    pAdj_Values <- dFrame_study$pAdj_Smoke
    
    plotList <- list()
    plotIndex <- 1
    
    for (j in 1:nrow(dFrame_study)) {
      
      pval <- pAdj_Values[j]
      p.label <- roundP(pval)
      bugName <- dFrame_study$BugName[j]
      tumorTime <- dFrame_study$Tumor[j]
      df2 <- df[df[,t] == bugName,]
      
      NormAbundance <- df2$NormalizedCounts
      Term <- df2$Treatment
      myFrame <- data.frame(NormAbundance, Term)
      title.lab <- paste0(level, " - ", bugName, " (", tumorTime, "-tumor)\nModel: ", model, "\nAdjusted ", p.label)
      
      plot <- ggplot(myFrame, aes(x=Term, y=NormAbundance))+
        geom_boxplot()+
        geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), pch=21, aes(fill=factor(Term)), show.legend = F, size = 2)+
        labs(title=title.lab, x=Treatment, y = "Log Normalized Abundance")
      
      plotList[[plotIndex]] <- plot
      plotIndex <- plotIndex + 1
      
    }
    
    
    plotNum <- length(plotList)
    
    PlotOut <- paste0(outputLevel, level, "_", study, "_Smoke_boxplots.pdf")
    
    pdf(PlotOut, width = 5, height = 5)
    for (i in plotNum) {
      print(plotList)
    } # for (i in plotNum)
    hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0("Treatment p values"))
    dev.off()
    
  }
  
} # for (level in levels)

summary <- data.frame(Level, Experiment, Sig_pVal)

outputFile <- paste0(output, "Smoke_v_NoSmoke_Significance_summary.tsv")
write.table(summary, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
