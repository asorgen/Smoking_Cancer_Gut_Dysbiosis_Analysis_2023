#Author: Alicia Sorgen
#Date: 09-23-21
#Description: Performs basic linear modeling and produces simple boxplots

## Libraries
library(stringr); packageVersion("stringr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Sep15"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("Abx-Pancreatic_Antibiotic_lm"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  study = commandArgs(trailingOnly = TRUE)
} else {
  study <- "Abx-Bladder"
  setwd(paste0(root, "/script"))
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

inputStudy = paste0(inputDir, study, "/")
message(inputStudy)

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
  norm.post <- norm[norm$TumorImp == "Post",]
  norm.NoSmoke <- norm.post[norm.post$SmokeExp == "No",]
  
  if (study == "Neo-Pancreatic") {
    norm.NoSmoke$Treatment <- factor(norm.NoSmoke$Treatment, levels = c("Control", "Neomycin only"))
    abx <- "Neomycin"
  } else {
    norm.NoSmoke$Treatment <- factor(norm.NoSmoke$Treatment, levels = c("Control", "Antibiotic only"))
    abx <- "Cocktail"
  }
  
  
  startTaxaIndex <- which(colnames(norm.NoSmoke)=="Study")+1
  
  bugName <- vector()
  p_AntibioticExp <- vector()
  AntibioticExp_direction <- vector()
  index <- 1
  
  ########## Perform linear modeling ##########
  for (i in startTaxaIndex:ncol(norm.NoSmoke)) {
    
    if (sum(norm.NoSmoke[,i]!=0) > nrow(norm.NoSmoke)*0.1) {
      
      bug <- norm.NoSmoke[,i]
      AntibioticExp <- norm.NoSmoke$AntibioticExp
      
      abx.lm <- lm(bug ~ AntibioticExp)
      sum.lm <- summary(abx.lm)
      
      p_AntibioticExp[index] <- sum.lm$coefficients[2,4]      
      
      bugName[index] <- names(norm.NoSmoke)[i]
      
      # Determine direction of higher abundance in AntibioticExp
      if(mean(norm.NoSmoke[norm.NoSmoke$AntibioticExp==abx,i])>mean(norm.NoSmoke[norm.NoSmoke$AntibioticExp=="None",i])){
        AntibioticExp_direction[index]="Higher with abx"
      }
      else if(mean(norm.NoSmoke[norm.NoSmoke$AntibioticExp==abx,i])<mean(norm.NoSmoke[norm.NoSmoke$AntibioticExp=="None",i])){
        AntibioticExp_direction[index]="Lower with abx"
      }
      else{
        AntibioticExp_direction[index]="Same"
      }
      
      index <- index + 1
      
    } # if (sum(norm.NoSmoke[,i]!=0) > nrow(norm.NoSmoke)*0.1)
    
  } # for (x in startTaxaIndex:ncol(norm.NoSmoke))
  
  dFrame <- data.frame(bugName, p_AntibioticExp, AntibioticExp_direction)
  
  # Performing BH adjustments on p values
  dFrame$pAdj_AntibioticExp<- p.adjust( dFrame$p_AntibioticExp, method = "BH")
  
  dFrame <- parseQIIME2Taxa(dFrame)
  
  dFilt1 <- dFrame[!(dFrame$bugName == "Other"),]
  TotalTaxa[count] <- nrow(dFilt1)
  dFilt <- dFilt1[dFilt1$pAdj_AntibioticExp < 0.05,]
  SigTaxa[count] <- nrow(dFilt)
  Level[count] <- level
  
  write.table(dFrame, paste0(outputLevel, level, "_", study, "_Antibiotic_lm.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
  
  ########## Generate boxplots ##########
  count <- count + 1
  rownames(dFrame) <- dFrame[,count]
  
  terms <- c("AntibioticExp")
  for (term in terms) {
    
    model <- paste0("lm(bug ~ ", term, ")")
    unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term))]
    dFrame.ordered <- dFrame[order(dFrame[,which(colnames(dFrame) == paste0("p_", term))]),]
    
    plotList <- list()
    plotIndex <- 1
    
    for (j in 1:nrow(dFrame.ordered)) {
      pAdj_Values <- dFrame.ordered[,which(colnames(dFrame.ordered) == paste0("pAdj_", term))]
      if (pAdj_Values[j]<0.05) {
        
        pval <- pAdj_Values[j]
        p.label <- roundP(pval)
        bugName <- rownames(dFrame.ordered)[j]
        bugFull <- dFrame.ordered$bugName[j]
        NormAbundance <- norm.NoSmoke[,which(colnames(norm.NoSmoke) == bugFull)]
        Term <- norm.NoSmoke[,which(colnames(norm.NoSmoke) == "Treatment")]
        myFrame <- data.frame(NormAbundance, Term)
        title.lab <- paste0(level, " - ", bugName, " (Post-tumor)\nModel: ", model, "\nAdjusted ", p.label)
        
        plot <- ggplot(myFrame, aes(x=Term, y=NormAbundance))+
          geom_boxplot()+
          geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), pch=21, aes(fill=factor(Term)), show.legend = F, size = 2)+
          labs(title=title.lab, x=term, y = "Log Normalized Abundance")
        
        plotList[[plotIndex]] <- plot
        plotIndex <- plotIndex + 1
      }
    }
    
    
    plotNum <- length(plotList)
    
    PlotOut <- paste0(outputLevel, level, "_", study, "_", term, "_Post-tumor_lm_boxplots.pdf")
    
    pdf(PlotOut, width = 5, height = 5)
    for (i in plotNum) {
      print(plotList)
    } # for (i in plotNum)
    hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term, " p values"))
    dev.off()
    
  } # for (term in terms)
  
} # for (level in levels)


summary <- data.frame(Level, SigTaxa, TotalTaxa)
tsvName <- "Antibiotic_lm_Summary.tsv"
write.table(summary, paste0(output, tsvName),sep="\t",quote = FALSE, row.names = FALSE)



