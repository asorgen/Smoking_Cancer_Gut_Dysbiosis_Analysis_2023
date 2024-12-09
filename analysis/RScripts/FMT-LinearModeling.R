#Author: Alicia Sorgen
#Date: 09-14-21
#Description: Performs basic linear modeling and produces simple boxplots

R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))
library(tidyr); message("tidyr:", packageVersion("tidyr"))
library(ggplot2); message("ggplot2:", packageVersion("ggplot2"))
library(data.table); message("data.table:", packageVersion("data.table"))

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Oct06"

root = paste0("~/BioLockJ_pipelines/FMTmouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("FMT-LinearModeling"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("Smoke-FMT", "NNK-FMT")
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CombineStudies"),"/output/"); message(inputDir)
moduleDir = dirname(getwd()); message(moduleDir)
output = file.path(moduleDir,"output/"); message(output)

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

levels = c("Phylum", "Class", "Order", "Family", "Genus")

Level <- vector()
Experiment <- vector()
NicExposure <- vector()
Sig_pVal <- vector()
smry <- 1

Level2 <- vector()
Experiment2 <- vector()
Tumor_Status <- vector()
Sig_pVal2 <- vector()
smry2 <- 1

Level3 <- vector()
Experiment3 <- vector()
NicExposure3 <- vector()
Sig_pVal3 <- vector()
smry3 <- 1

for (level in levels) {
  
  outputLevel = paste0(output, level, "/"); message(outputLevel)
  dir.create(outputLevel, showWarnings = FALSE)
  
  inputFile <- paste0(level, "_CombineStudies_ALL.tsv")
  fullPath <- paste0(inputDir, inputFile); message(fullPath)
  
  myT <- read.table(file = fullPath, sep="\t", header = TRUE)
  
  ##### Pre- v Post-tumor Comparisons #####
  myT_Recipient <- myT[myT$FMTstatus == "Recipient",]
  
  BugName <- vector()
  Study <- vector()
  p_Tumor <- vector()
  s_Tumor <- vector()
  NicotineExposure <- vector()
  index <- 1
  
  for (study in studies) {
    
    nicType <- strsplit(study,"-")[[1]][1]    
    nicexp<- c(paste0(nicType, "-free"), paste0(nicType, "-exposed"))
    
    df <- myT_Recipient[ myT_Recipient$Study == study, ]
    
    for (nict in nicexp) {
      
      df2 <- df[df$SmokeExposure %like% nict,]
      t <- which( colnames(df2) == level )
      
      for (bug in unique(df2[,t])) {
        
        if ( !(bug == "Other") ) {
          
          df_taxa <- df2[df2[,t] == bug,]
          NormalizedCounts <- df_taxa$NormalizedCounts
          TumorStatus <- df_taxa$FullMouseDescription
          
          lm <- lm(NormalizedCounts ~ TumorStatus)
          sm <- summary(lm)
          
          BugName[index] <- bug
          Study[index] <- study
          NicotineExposure[index] <- nict
          p_Tumor[index] <- sm$coefficients[2,4]
          s_Tumor[index] <- sm$coefficients[2,1]
          index <- index + 1
          
        } # if ( !(bug == "Other") )
        
      } # for (bug in unique(df[,t]))
      
    }
    
  } # for (study in studies)
  
  dFrame <- data.frame(BugName, Study, NicotineExposure, s_Tumor, p_Tumor)
  dFrame$pAdj_Tumor <- p.adjust(dFrame$p_Tumor, method = "BH")
  dFrame <- na.omit(dFrame)
  dFrame <- dFrame[order(dFrame$p_Tumor),]
  
  outputFile <- paste0(outputLevel, level, "_Combined_LM_TumorStatusComparisons_ordered.tsv")
  write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
  
  dFrameFiltered <- dFrame[dFrame$pAdj_Tumor < 0.05,]
  
  for (study in studies) {
    
    dFrame3 <- dFrameFiltered[dFrameFiltered$Study == study,]
    
    for (j in unique(dFrame3$NicotineExposure)) {
      
      dFrame4 <- dFrame3[dFrame3$NicotineExposure == j,]
      Level[smry] <- level
      Experiment[smry] <- study
      NicExposure[smry] <- j
      Sig_pVal[smry] <- nrow(dFrame4)
      smry <- smry + 1
      
    }
  }
  
  for (study in studies) {

    nicType <- strsplit(study,"-")[[1]][1]
    nicexp<- c(paste0(nicType, "-free"), paste0(nicType, "-exposed"))

    df <- myT_Recipient[ myT_Recipient$Study == study, ]
    dFrame_study <- dFrame[ dFrame$Study == study, ]

    for (nict in nicexp) {

      df_exp <- df[ df$SmokeExposure %like% nict, ]
      dFrame_nic <- dFrame_study[dFrame_study$NicotineExposure %like% nict,]
      model <- paste0("lm(NormalizedCounts ~ TumorStatus)")
      unadj <- dFrame_nic$p_Tumor
      pAdj_Values <- dFrame_nic$pAdj_Tumor

      plotList <- list()
      plotIndex <- 1

      for (j in 1:nrow(dFrame_nic)) {

        pval <- pAdj_Values[j]
        p.label <- roundP(pval)
        bugName <- dFrame_nic$BugName[j]
        Nic <- dFrame_nic$NicotineExposure[j]
        df2 <- df_exp[df_exp[,t] == bugName,]

        NormAbundance <- df2$NormalizedCounts
        Term <- df2$FullMouseDescription
        myFrame <- data.frame(NormAbundance, Term)
        title.lab <- paste0(level, " - ", bugName, "\nModel: ", model, "\nAdjusted ", p.label)

        plot <- ggplot(myFrame, aes(x=Term, y=NormAbundance))+
          geom_boxplot()+
          geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), pch=21, aes(fill=factor(Term)), show.legend = F, size = 2)+
          labs(title=title.lab, x="Mouse", y = "Log Normalized Abundance")

        plotList[[plotIndex]] <- plot
        plotIndex <- plotIndex + 1

      }


      plotNum <- length(plotList)

      PlotOut <- paste0(outputLevel, level, "_", study, "_", nict, "_Pre-v-Post-tumor_boxplots.pdf")

      pdf(PlotOut, width = 5, height = 5)
      for (i in plotNum) {
        print(plotList)
      } # for (i in plotNum)
      hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0("Treatment p values"))
      dev.off()

    }

  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##### Smoke-free v Smoke-exposed Comparisons #####
  myT_Recipient <- myT[myT$FMTstatus == "Recipient",]
  
  BugName <- vector()
  Study <- vector()
  p_Nic <- vector()
  s_Nic <- vector()
  TumorStatus <- vector()
  index <- 1
  
  for (study in studies) {
    
    tumorStatus<- unique(myT_Recipient$TumorStatus)
    
    df <- myT_Recipient[ myT_Recipient$Study == study, ]
    
    for (tumor in tumorStatus) {
      
      df2 <- df[df$TumorStatus %like% tumor,]
      t <- which( colnames(df2) == level )
      
      for (bug in unique(df2[,t])) {
        
        if ( !(bug == "Other") ) {
          
          df_taxa <- df2[df2[,t] == bug,]
          NormalizedCounts <- df_taxa$NormalizedCounts
          NicotineExposure <- df_taxa$FullMouseDescription
          
          lm <- lm(NormalizedCounts ~ NicotineExposure)
          sm <- summary(lm)
          
          BugName[index] <- bug
          Study[index] <- study
          TumorStatus[index] <- tumor
          p_Nic[index] <- sm$coefficients[2,4]
          s_Nic[index] <- sm$coefficients[2,1]
          index <- index + 1
          
        } # if ( !(bug == "Other") )
        
      } # for (bug in unique(df[,t]))
      
    }
    
  } # for (study in studies)
  
  dFrame <- data.frame(BugName, Study, TumorStatus, s_Nic, p_Nic)
  dFrame$pAdj_Nic <- p.adjust(dFrame$p_Nic, method = "BH")
  dFrame <- na.omit(dFrame)
  dFrame <- dFrame[order(dFrame$p_Nic),]
  outputFile <- paste0(outputLevel, level, "_Combined_LM_NicotineExposureComparisons_ordered.tsv")
  write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
  
  dFrameFiltered <- dFrame[dFrame$pAdj_Nic < 0.05,]
  
  for (study in studies) {
    dFrame3 <- dFrameFiltered[dFrameFiltered$Study == study,]
    
    for (j in unique(dFrame3$TumorStatus)) {
      
      dFrame4 <- dFrame3[dFrame3$TumorStatus == j,]
      Level2[smry2] <- level
      Experiment2[smry2] <- study
      Tumor_Status[smry2] <- j
      Sig_pVal2[smry2] <- nrow(dFrame4)
      smry2 <- smry2 + 1
      
    }
  }
  
  
  for (study in studies) {

    nicType <- strsplit(study,"-")[[1]][1]

    tumorStatus<- unique(myT_Recipient$TumorStatus)
    tumorStatus <- factor(tumorStatus, levels = c("Pre-tumor", "Post-tumor"))

    df <- myT_Recipient[ myT_Recipient$Study == study, ]
    dFrame_study <- dFrame[ dFrame$Study == study, ]

    for (tumor in tumorStatus) {

      df_exp <- df[ df$FullMouseDescription %like% tumor, ]
      dFrame_nic <- dFrame_study[dFrame_study$TumorStatus %like% tumor,]
      model <- paste0("lm(NormalizedCounts ~ NicotineExposure)")
      unadj <- dFrame_nic$p_Nic
      pAdj_Values <- dFrame_nic$pAdj_Nic

      plotList <- list()
      plotIndex <- 1

      for (j in 1:nrow(dFrame_nic)) {

        pval <- pAdj_Values[j]
        p.label <- roundP(pval)
        bugName <- dFrame_nic$BugName[j]
        # Nic <- dFrame_nic$NicotineExposure[j]
        df2 <- df_exp[df_exp[,t] == bugName,]

        NormAbundance <- df2$NormalizedCounts
        Term <- df2$FullMouseDescription
        myFrame <- data.frame(NormAbundance, Term)
        title.lab <- paste0(level, " - ", bugName, "\nModel: ", model, "\nAdjusted ", p.label)

        plot <- ggplot(myFrame, aes(x=Term, y=NormAbundance))+
          geom_boxplot()+
          geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), pch=21, aes(fill=factor(Term)), show.legend = F, size = 2)+
          labs(title=title.lab, x="Mouse", y = "Log Normalized Abundance")

        plotList[[plotIndex]] <- plot
        plotIndex <- plotIndex + 1

      }


      plotNum <- length(plotList)

      PlotOut <- paste0(outputLevel, level, "_", study, "_", tumor, "_", nicType, "-free_v_", nicType, "-exposed_boxplots.pdf")

      pdf(PlotOut, width = 5, height = 5)
      for (i in plotNum) {
        print(plotList)
      } # for (i in plotNum)
      hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0("Treatment p values"))
      dev.off()

    }

  }

  
  
  
  
  
  ##### Donor v Recipient Comparisons #####
  myT_2 <- myT[!(myT$TumorStatus == "Post-tumor"),]
  
  BugName <- vector()
  Study <- vector()
  p_FMTstatus <- vector()
  s_FMTstatus <- vector()
  NicotineExposure <- vector()
  index <- 1
  
  for (study in studies) {
    
    
    df <- myT_2[ myT_2$Study == study, ]
    
    nicType <- strsplit(study,"-")[[1]][1]    
    
    if (nicType == "Smoke") {
      Terms1<- c(paste0(nicType, "-free"), paste0(nicType, "-exposed"))
    } else {
      Terms1<- c(paste0(nicType, "-exposed"))
    }
    
    for (term1 in Terms1) {
      
      df2 <- df[df$SmokeExposure %like% term1,]
      t <- which( colnames(df2) == level )
      
      for (bug in unique(df2[,t])) {
        
        if ( !(bug == "Other") ) {
          
          df_taxa <- df2[df2[,t] == bug,]
          NormalizedCounts <- df_taxa$NormalizedCounts
          Variable <- df_taxa$FMTstatus
          
          lm <- lm(NormalizedCounts ~ Variable)
          sm <- summary(lm)
          
          BugName[index] <- bug
          Study[index] <- study
          NicotineExposure[index] <- term1
          p_FMTstatus[index] <- sm$coefficients[2,4]
          s_FMTstatus[index] <- sm$coefficients[2,1]
          index <- index + 1
          
        } # if ( !(bug == "Other") )
        
      } # for (bug in unique(df[,t]))
      
    }
    
  } # for (study in studies)
  
  dFrame <- data.frame(BugName, Study, NicotineExposure, p_FMTstatus, s_FMTstatus)
  dFrame$pAdj_FMTstatus <- p.adjust(dFrame$p_FMTstatus, method = "BH")
  dFrame <- na.omit(dFrame)
  dFrame <- dFrame[order(dFrame$p_FMTstatus),]
  outputFile <- paste0(outputLevel, level, "_Combined_LM_Donor_RecipientComparisons_ordered.tsv")
  write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
  
  dFrameFiltered <- dFrame[dFrame$pAdj_FMTstatus < 0.05,]
  
  for (study in studies) {
    dFrame3 <- dFrameFiltered[dFrameFiltered$Study == study,]
    
    for (j in unique(dFrame3$NicotineExposure)) {
      
      dFrame4 <- dFrame3[dFrame3$NicotineExposure == j,]
      Level3[smry3] <- level
      Experiment3[smry3] <- study
      NicExposure3[smry3] <- j
      Sig_pVal3[smry3] <- nrow(dFrame4)
      smry3 <- smry3 + 1
      
    }
  }
  
  
  for (study in studies) {
    
    nicType <- strsplit(study,"-")[[1]][1]
    
    if (nicType == "Smoke") {
      Terms1<- c(paste0(nicType, "-free"), paste0(nicType, "-exposed"))
    } else {
      Terms1<- c(paste0(nicType, "-exposed"))
    }
    
    df <- myT_2[ myT_2$Study == study, ]
    df <- df[!(df$TumorStatus == "Post-tumor"),]
    
    dFrame_study <- dFrame[ dFrame$Study == study, ]
    
    for (term1 in Terms1) {
      
      df_exp <- df[ df$SmokeExposure %like% term1, ]
      dFrame_2 <- dFrame_study[dFrame_study$NicotineExposure %like% term1,]
      model <- paste0("lm(NormalizedCounts ~ FMTstatus)")
      unadj <- dFrame_2$p_FMTstatus
      pAdj_Values <- dFrame_2$pAdj_FMTstatus
      
      plotList <- list()
      plotIndex <- 1
      
      for (j in 1:nrow(dFrame_2)) {
        
        pval <- pAdj_Values[j]
        p.label <- roundP(pval)
        bugName <- dFrame_2$BugName[j]
        # Nic <- dFrame_2$NicotineExposure[j]
        df2 <- df_exp[df_exp[,t] == bugName,]
        
        NormAbundance <- df2$NormalizedCounts
        Term <- df2$FMTstatus
        myFrame <- data.frame(NormAbundance, Term)
        title.lab <- paste0(level, " - ", bugName, "\nModel: ", model, "\nAdjusted ", p.label)
        
        plot <- ggplot(myFrame, aes(x=Term, y=NormAbundance))+
          geom_boxplot()+
          geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), pch=21, aes(fill=factor(Term)), show.legend = F, size = 2)+
          labs(title=title.lab, x="Mouse", y = "Log Normalized Abundance")
        
        plotList[[plotIndex]] <- plot
        plotIndex <- plotIndex + 1
        
      }
      
      
      plotNum <- length(plotList)
      
      PlotOut <- paste0(outputLevel, level, "_", study, "_", term1, "_Donor_v_Recipient_boxplots.pdf")
      
      pdf(PlotOut, width = 5, height = 5)
      for (i in plotNum) {
        print(plotList)
      } # for (i in plotNum)
      hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0("Treatment p values"))
      dev.off()
      
    }
    
  }
  
  
} # for (level in levels)

summary <- data.frame(Level, Experiment, NicExposure, Sig_pVal)
outputFile <- paste0(output, "Pre_v_Post-tumor_Significance_summary.tsv")
write.table(summary, outputFile, sep="\t", quote = FALSE, row.names = FALSE)


summary <- data.frame(Level2, Experiment2, Tumor_Status, Sig_pVal2)
outputFile <- paste0(output, "Nic_v_NoNic_Significance_summary.tsv")
write.table(summary, outputFile, sep="\t", quote = FALSE, row.names = FALSE)

summary <- data.frame(Level3, Experiment3, NicExposure3, Sig_pVal3)
outputFile <- paste0(output, "Donor_v_Recipient_Significance_summary.tsv")
write.table(summary, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
