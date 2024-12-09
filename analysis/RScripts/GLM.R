#Author: Alicia Sorgen
#Date: 09-21-21
#Description: Performs basic linear modeling and produces simple boxplots

## Libraries
library(stringr)
library(tidyr)
library(ggplot2)

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Sep15"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("Neo-PancreaticLM"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  study = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  study <- "Neo-Pancreatic"
  
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

inputStudy = paste0(inputDir, study, "/")

levels = c("Phylum", "Class", "Order", "Family", "Genus")

Level <- vector()
SigTaxa <- vector()
TotalTaxa <- vector()
count <- 1

for (level in levels) {
  
  
  outputLevel = paste0(output, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  inputFile <- paste0(inputStudy, level, "_normCounts.tsv")
  norm <- read.table(file = inputFile, sep="\t", header = TRUE)
  norm$Treatment <- gsub(pattern = " only", replacement = "", norm$Treatment)
  
  if (study == "Neo-Pancreatic") {
    norm$Treatment <- gsub(pattern = "Neomycin", replacement = "Neo", norm$Treatment)
    norm$Treatment <- factor(norm$Treatment, levels = c("Control", "Smoke", "Neo", "Smoke+Neo"))
    abx <- "Neomycin"
  } else {
    norm$Treatment <- gsub(pattern = "Antibiotic", replacement = "Abx", norm$Treatment)
    norm$Treatment <- factor(norm$Treatment, levels = c("Control", "Smoke", "Abx", "Smoke+Abx"))
    abx <- "Cocktail"
    }
  
  norm.post <- norm[norm$TumorImp == "Post",]
  
  startTaxaIndex <- which(colnames(norm.post)=="Study")+1
  
  bugName <- vector()
  p_SmokeExp <- vector()
  p_AntibioticExp <- vector()
  p_Interaction <- vector()
  SmokeExp_direction <- vector()
  AntibioticExp_direction <- vector()
  index <- 1
  
  ########## Perform linear modeling ##########
  for (i in startTaxaIndex:ncol(norm.post)) {
    
    if (sum(norm.post[,i]!=0) > nrow(norm.post)*0.1) {
      
      SmokeExp <- norm.post$SmokeExp
      AntibioticExp <- norm.post$AntibioticExp
      
      myLm <- lm(norm.post[,i] ~ SmokeExp * AntibioticExp)
      myAnova <- anova(myLm)
      
      p_SmokeExp[index] <- myAnova$"Pr(>F)"[1]
      p_AntibioticExp[index] <- myAnova$"Pr(>F)"[2]
      p_Interaction[index] <- myAnova$"Pr(>F)"[3]
      
      bugName[index] <- names(norm.post)[i]
      
      # Determine direction of higher abundance in SmokeExp
      if(mean(norm.post[norm.post$SmokeExp=="Yes",i])>mean(norm.post[norm.post$SmokeExp=="No",i])){
        SmokeExp_direction[index]="Higher with smoke"
      }
      else if(mean(norm.post[norm.post$SmokeExp=="Yes",i])<mean(norm.post[norm.post$SmokeExp=="No",i])){
        SmokeExp_direction[index]="Lower with smoke"
      }
      else{
        SmokeExp_direction[index]="Same"
      }
      
      # Determine direction of higher abundance in AntibioticExp
      if(mean(norm.post[norm.post$AntibioticExp==abx,i])>mean(norm.post[norm.post$AntibioticExp=="None",i])){
        AntibioticExp_direction[index]="Higher with antibiotics"
      }
      else if(mean(norm.post[norm.post$AntibioticExp==abx,i])<mean(norm.post[norm.post$AntibioticExp=="None",i])){
        AntibioticExp_direction[index]="Lower with antibiotics"
      }
      else{
        AntibioticExp_direction[index]="Same"
      }
      
      
      index <- index + 1
      
    } # if (sum(norm.post[,i]!=0) > nrow(norm.post)*0.1)
    
  } # for (x in startTaxaIndex:ncol(norm.post))
  
  dFrame <- data.frame(bugName, p_SmokeExp,SmokeExp_direction,p_AntibioticExp,p_Interaction)
  
  # Performing BH adjustments on p values
  dFrame$pAdj_SmokeExp<- p.adjust( dFrame$p_SmokeExp, method = "BH")
  dFrame$pAdj_AntibioticExp<- p.adjust( dFrame$p_AntibioticExp, method = "BH")
  dFrame$pAdj_Interaction<- p.adjust( dFrame$p_Interaction, method = "BH")
  
  dFrame <- parseQIIME2Taxa(dFrame)
  
  dFilt1 <- dFrame[!(dFrame$bugName == "Other"),]
  TotalTaxa[count] <- nrow(dFilt1)
  dFilt <- dFilt1[dFilt1$pAdj_Interaction < 0.05,]
  SigTaxa[count] <- nrow(dFilt)
  Level[count] <- level
  
  write.table(dFrame, paste0(outputLevel, level, "_", study, "_Smoke_Antibiotic_Interaction_LM.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
  
  ########## Smoke boxplots ##########
  count <- count + 1
  rownames(dFrame) <- dFrame[,count]
  
  term1 <- "SmokeExp"
  term2 <- "AntibioticExp"
  unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term1))]
  
  dFrame.ordered <- dFrame[order(dFrame[,which(colnames(dFrame) == paste0("p_", term1))]),]
  
  plotList <- list()
  plotIndex <- 1
  
  for (j in 1:nrow(dFrame.ordered)) {
    pAdj_Values <- dFrame.ordered[,which(colnames(dFrame.ordered) == paste0("pAdj_", term1))]
    if (pAdj_Values[j]<0.05) {
      
      pval <- pAdj_Values[j]
      p.label <- roundP(pval)
      bugName <- rownames(dFrame.ordered)[j]
      bugFull <- dFrame.ordered$bugName[j]
      NormAbundance <- norm.post[,which(colnames(norm.post) == bugFull)]
      Term1 <- norm.post[,which(colnames(norm.post) == term1)]
      Term2 <- norm.post[,which(colnames(norm.post) == term2)]
      myFrame <- data.frame(NormAbundance, Term1, Term2)
      title.lab <- paste0(level, " - ", bugName, " (Post-tumor)\nModel: lm(bug ~ SmokeExp * AntibioticExp)\nAdjusted ", p.label)
      
      plot <- ggplot(myFrame, aes(x=Term1, y=NormAbundance))+
        geom_boxplot()+
        geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
                   pch=21, aes(fill=factor(Term2)), show.legend = T, size = 2)+
        labs(title=title.lab, x=term1, y = "Log Normalized Abundance")+
        guides(fill=guide_legend(title=term2))
      
      plotList[[plotIndex]] <- plot
      plotIndex <- plotIndex + 1
    }
  }
  
  
  plotNum <- length(plotList)
  
  PlotOut <- paste0(outputLevel, level, "_", study, "_", term1, "_Post-tumor_GLM_boxplots.pdf")
  
  pdf(PlotOut, width = 5, height = 5)
  for (i in plotNum) {
    print(plotList)
  } # for (i in plotNum)
  hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term1, " p values"))
  dev.off()
  
  ########## Antibiotic boxplots ##########
  term2 <- "SmokeExp"
  term1 <- "AntibioticExp"
  unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term1))]
  
  dFrame.ordered <- dFrame[order(dFrame[,which(colnames(dFrame) == paste0("p_", term1))]),]
  
  plotList <- list()
  plotIndex <- 1
  
  for (j in 1:nrow(dFrame.ordered)) {
    pAdj_Values <- dFrame.ordered[,which(colnames(dFrame.ordered) == paste0("pAdj_", term1))]
    if (pAdj_Values[j]<0.05) {
      
      pval <- pAdj_Values[j]
      p.label <- roundP(pval)
      bugName <- rownames(dFrame.ordered)[j]
      bugFull <- dFrame.ordered$bugName[j]
      NormAbundance <- norm.post[,which(colnames(norm.post) == bugFull)]
      Term1 <- norm.post[,which(colnames(norm.post) == term1)]
      Term2 <- norm.post[,which(colnames(norm.post) == term2)]
      myFrame <- data.frame(NormAbundance, Term1, Term2)
      title.lab <- paste0(level, " - ", bugName, " (Post-tumor)\nModel: lm(bug ~ SmokeExp * AntibioticExp)\nAdjusted ", p.label)
      
      plot <- ggplot(myFrame, aes(x=Term1, y=NormAbundance))+
        geom_boxplot()+
        geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
                   pch=21, aes(fill=factor(Term2)), show.legend = T, size = 2)+
        labs(title=title.lab, x=term1, y = "Log Normalized Abundance")+
        guides(fill=guide_legend(title=term2))
      
      plotList[[plotIndex]] <- plot
      plotIndex <- plotIndex + 1
    }
  }
  
  
  plotNum <- length(plotList)
  
  PlotOut <- paste0(outputLevel, level, "_", study, "_", term1, "_Post-tumor_GLM_boxplots.pdf")
  
  pdf(PlotOut, width = 5, height = 5)
  for (i in plotNum) {
    print(plotList)
  } # for (i in plotNum)
  hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term1, " p values"))
  dev.off()
  
  
  ########## Generate interaction boxplots ##########
  term1 <- "Interaction"
  term2 <- "Treatment"
  
  unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term1))]
  dFrame.ordered <- dFrame[order(dFrame[,which(colnames(dFrame) == paste0("p_", term1))]),]
  
  plotList <- list()
  plotIndex <- 1
  
  for (j in 1:nrow(dFrame.ordered)) {
    pAdj_Values <- dFrame.ordered[,which(colnames(dFrame.ordered) == paste0("pAdj_", term1))]
    if (pAdj_Values[j]<0.05) {
      
      pval <- pAdj_Values[j]
      p.label <- roundP(pval)
      bugName <- rownames(dFrame.ordered)[j]
      bugFull <- dFrame.ordered$bugName[j]
      NormAbundance <- norm.post[,which(colnames(norm.post) == bugFull)]
      # Term1 <- norm.post[,which(colnames(norm.post) == term1)]
      Term2 <- norm.post[,which(colnames(norm.post) == term2)]
      myFrame <- data.frame(NormAbundance, Term2)
      title.lab <- paste0(level, " - ", bugName, " (Post-tumor)\nModel: lm(bug ~ SmokeExp * AntibioticExp)\nAdjusted ", p.label)
      
      plot <- ggplot(myFrame, aes(x=Term2, y=NormAbundance))+
        geom_boxplot()+
        geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), pch=21, aes(fill=factor(Term2)), show.legend = F, size = 2)+
        labs(title=title.lab, x=term2, y = "Log Normalized Abundance")
      
      plotList[[plotIndex]] <- plot
      plotIndex <- plotIndex + 1
    }
  }
  
  
  plotNum <- length(plotList)
  
  PlotOut <- paste0(outputLevel, level, "_", study, "_", term1, "_Post-tumor_GLM_boxplots.pdf")
  
  pdf(PlotOut, width = 5, height = 5)
  for (i in plotNum) {
    print(plotList)
  } # for (i in plotNum)
  hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term1, " p values"))
  dev.off()
  
  
  
} # for (level in levels)

summary <- data.frame(Level, SigTaxa, TotalTaxa)
tsvName <- "Interaction_lm_Summary.tsv"
write.table(summary, paste0(output, tsvName),sep="\t",quote = FALSE, row.names = FALSE)

