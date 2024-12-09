#Author: Alicia Sorgen
#Date: 08-19-21
#Description: Performs basic linear modeling and produces simple boxplots

## Libraries
library(stringr)
library(tidyr)
library(ggplot2)

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
date = "2021Aug16"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("Abx-BladderLM"), full.names=TRUE)

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
output = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

study <- "Abx-Bladder"
inputStudy = paste0(inputDir, "/", study, "/")

levels = c("Phylum", "Class", "Order", "Family", "Genus")
count <- 1

for (level in levels) {
  
  outputLevel = paste0(output, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  count <- 1 + count
  inputFile <- paste0(inputStudy, level, "_normCounts.tsv")
  norm <- read.table(file = inputFile, sep="\t", header = TRUE)
  
  startTaxaIndex <- which(colnames(norm)=="Study")+1
  
  bugName <- vector()
  p_SmokeExp <- vector()
  p_AntibioticExp <- vector()
  p_Interaction <- vector()
  SmokeExp_direction <- vector()
  AntibioticExp_direction <- vector()
  index <- 1
  
  for (i in startTaxaIndex:ncol(norm)) {
    
    if (sum(norm[,i]!=0) > nrow(norm)*0.1) {
      
      SmokeExp <- norm$SmokeExp
      AntibioticExp <- norm$AntibioticExp
      
      myLm <- lm(norm[,i] ~ SmokeExp * AntibioticExp)
      myAnova <- anova(myLm)
      
      p_SmokeExp[index] <- myAnova$"Pr(>F)"[1]
      p_AntibioticExp[index] <- myAnova$"Pr(>F)"[2]
      p_Interaction[index] <- myAnova$"Pr(>F)"[3]

      bugName[index] <- names(norm)[i]
      
      if(mean(norm[norm$SmokeExp=="Yes",i])>mean(norm[norm$SmokeExp=="No",i])){
        SmokeExp_direction[index]="Higher in SmokeExp"
      }
      else if(mean(norm[norm$SmokeExp=="Yes",i])<mean(norm[norm$SmokeExp=="No",i])){
        SmokeExp_direction[index]="Lower in SmokeExp"
      }
      else{
        SmokeExp_direction[index]="Same"
      }
      
      if(mean(norm[norm$AntibioticExp=="Cocktail",i])>mean(norm[norm$AntibioticExp=="None",i])){
        AntibioticExp_direction[index]="Higher in AntibioticExp"
      }
      else if(mean(norm[norm$AntibioticExp=="Cocktail",i])<mean(norm[norm$AntibioticExp=="None",i])){
        AntibioticExp_direction[index]="Lower in AntibioticExp"
      }
      else{
        AntibioticExp_direction[index]="Same"
      }
      
      index <- index + 1

    } # if (sum(norm[,i]!=0) > nrow(norm)*0.1)
    
  } # for (x in startTaxaIndex:ncol(norm))
  
  dFrame <- data.frame(bugName, p_SmokeExp,SmokeExp_direction,p_AntibioticExp,AntibioticExp_direction,p_Interaction)
  # Performing BH adjustments on p values
  dFrame$pAdj_SmokeExp<- p.adjust( dFrame$p_SmokeExp, method = "BH")
  dFrame$pAdj_AntibioticExp<- p.adjust( dFrame$p_AntibioticExp, method = "BH")
  dFrame$pAdj_Interaction<- p.adjust( dFrame$p_Interaction, method = "BH")
  
  dFrame <- parseQIIME2Taxa(dFrame)
  
  write.table(dFrame, paste0(outputLevel, level, "_", study, "_LinearModelResults.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
  
  rownames(dFrame) <- dFrame[,count]
  terms <- c("SmokeExp", "AntibioticExp")
  for (term in terms) {
    
    unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term))]
    
    scatPlotOut <- paste0(outputLevel, level, "_", study, "_", term, "_boxplots.pdf")
    pdf(scatPlotOut, width = 5, height = 5)
    for (j in 1:nrow(dFrame)) {
      pAdj_Values <- dFrame[,which(colnames(dFrame) == paste0("pAdj_", term))]
      if (pAdj_Values[j]<0.05) {
        pval <- round(pAdj_Values[j], digits = 3)
        
        pval <- roundP(pval)
        p.label <- paste0("(adjusted ", pval, ")")
        
        bugName <- rownames(dFrame)[j]
        bugFull <- dFrame$bugName[j]
        NormAbundance <- norm[,which(colnames(norm) == bugFull)]
        Term <- norm[,which(colnames(norm) == term)]
        myFrame <- data.frame(NormAbundance, Term)
        par(cex.axis=0.75, cex.main = 0.75)
        boxplot(NormAbundance ~ Term, main = paste(level, "-", bugName, "\n", p.label, sep = " "), cex.main = 0.75,
                xlab = term, ylab = "Normalized Abundance")
        stripchart(NormAbundance ~ Term,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
      }
    }
    hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term, " p values"))
    dev.off()
    
    } # for (term in terms)
  
  
  
  
  terms <- c("Interaction")
  int1 <- "SmokeExp"
  int2 <- "AntibioticExp"
  
  for (term in terms) {
    term2 <- "Treatment"
    
    unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term))]
    
    scatPlotOut <- paste0(outputLevel, level, "_", study, "_", int1, "_", int2, "_", term, "_boxplots.pdf")
    pdf(scatPlotOut, width = 5, height = 5)
    for (j in 1:nrow(dFrame)) {
      pAdj_Values <- dFrame[,which(colnames(dFrame) == paste0("pAdj_", term))]
      if (pAdj_Values[j]<0.05) {
        pval <- round(pAdj_Values[j], digits = 3)
        
        pval <- roundP(pval)
        p.label <- paste0("(adjusted ", pval, ")")
        
        bugName <- rownames(dFrame)[j]
        bugFull <- dFrame$bugName[j]
        NormAbundance <- norm[,which(colnames(norm) == bugFull)]
        Term <- norm[,which(colnames(norm) == term2)]
        myFrame <- data.frame(NormAbundance, Term)
        title.lab <- paste0(level, " - ", bugName, "\n", int1, " x ", int2, " ", term, " ", p.label)
        par(cex.axis=0.75, cex.main = 0.75)
        boxplot(NormAbundance ~ Term, main = title.lab, cex.main = 0.75,
                xlab = term2, ylab = "Normalized Abundance")
        stripchart(NormAbundance ~ Term,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
      }
    }
    hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term, " p values"))
    dev.off()
    
  } # for (term in terms)
  
  

} # for (level in levels)

# Error


