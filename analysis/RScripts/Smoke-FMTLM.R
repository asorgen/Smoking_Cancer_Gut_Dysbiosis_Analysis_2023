#Date: 08-19-21
#Description: Performs basic linear modeling and produces simple boxplots

## Libraries
library(stringr)
library(tidyr)
library(ggplot2)

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Aug16"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("Smoke-FMTLM"), full.names=TRUE)

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

study <- "Smoke-FMT"
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
  
  norm.tumor <- norm[(norm$FMTstatus == "Recipient"), ]
  for (tumor.status in c("Pre", "Post")) {
    
    norm.tumor2 <- norm.tumor[(norm.tumor$TumorImp == tumor.status), ]
    
    bugName <- vector()
    p_SmokeExp <- vector()
    s_SmokeExp <- vector()
    SmokeExp_direction <- vector()
    index <- 1
    term <- "SmokeExp"
    for (i in startTaxaIndex:ncol(norm.tumor2)) {
      
      if (sum(norm.tumor2[,i]!=0) > nrow(norm.tumor2)*0.1) {
        
        SmokeExp <- norm.tumor2$SmokeExp
        
        myLm <- lm(norm.tumor2[,i] ~ SmokeExp)
        myAnova <- anova(myLm)
        
        s_SmokeExp[index] <- myLm$coefficients[2]
        p_SmokeExp[index] <- myAnova$"Pr(>F)"[1]
        
        bugName[index] <- names(norm.tumor2)[i]
        
        if(mean(norm.tumor2[norm.tumor2$SmokeExp=="Yes",i])>mean(norm.tumor2[norm.tumor2$SmokeExp=="No",i])){
          SmokeExp_direction[index]="Higher in SmokeExp"
        }
        else if(mean(norm.tumor2[norm.tumor2$SmokeExp=="Yes",i])<mean(norm.tumor2[norm.tumor2$SmokeExp=="No",i])){
          SmokeExp_direction[index]="Lower in SmokeExp"
        }
        else{
          SmokeExp_direction[index]="Same"
        }
        
        index <- index + 1
        
      } # if (sum(norm.tumor2[,i]!=0) > nrow(norm.tumor2)*0.1)
      
    } # for (x in startTaxaIndex:ncol(norm.tumor2))
    
    dFrame <- data.frame(bugName, s_SmokeExp, p_SmokeExp,SmokeExp_direction)
    # Performing BH adjustments on p values
    dFrame$pAdj_SmokeExp<- p.adjust( dFrame$p_SmokeExp, method = "BH")
    dFrame <- parseQIIME2Taxa(dFrame)
    write.table(dFrame, paste0(outputLevel, level, "_", study, "_", tumor.status, "-tumor_", term, "_LinearModelResults.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    
    rownames(dFrame) <- dFrame[,count]
    
    terms <- c("SmokeExp")
    for (term in terms) {
      
      unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term))]
      
      boxPlotOut <- paste0(outputLevel, level, "_", study, "_", tumor.status, "-tumor_", term, "_boxplots.pdf")
      pdf(boxPlotOut, width = 5, height = 5)
      for (j in 1:nrow(dFrame)) {
        pAdj_Values <- dFrame[,which(colnames(dFrame) == paste0("pAdj_", term))]
        if (pAdj_Values[j]<0.05) {
          pval <- round(pAdj_Values[j], digits = 3)
          
          pval <- roundP(pval)
          p.label <- paste0("(adjusted ", pval, ")")
          
          bugName <- rownames(dFrame)[j]
          bugFull <- dFrame$bugName[j]
          NormAbundance <- norm.tumor2[,which(colnames(norm.tumor2) == bugFull)]
          Term <- norm.tumor2[,which(colnames(norm.tumor2) == term)]
          myFrame <- data.frame(NormAbundance, Term)
          title.lab <- paste0(level, " - ", bugName, 
                             "\n", tumor.status, "-tumor ", term, " ", p.label)
          
          par(cex.axis=0.5, cex.main = 0.75)
          boxplot(NormAbundance ~ Term, main = title.lab, cex.main = 0.75,
                  xlab = term, ylab = "Normalized Abundance")
          stripchart(NormAbundance ~ Term,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
        }
      }
      hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term, " p values"))
      dev.off()
      
    } # for (term in terms)
  }

  term <- "TumorImp"
  Smoke.status <- c("Yes", "No")
  norm.rec <- norm[!(norm$FMTstatus == "Donor"),]
  for (i in Smoke.status) {
    
    if (i == "Yes") {
      Smoke.label <- "Smoke-exposed"
    } else {
      Smoke.label <- "Smoke-free"
    }
    norm.rec2 <- norm.rec[norm.rec$SmokeExp %in% i,]
    
    bugName <- vector()
    p_TumorImp <- vector()
    s_TumorImp <- vector()
    TumorImp_direction <- vector()
    index <- 1
    
    for (i in startTaxaIndex:ncol(norm.rec2)) {
      
      if (sum(norm.rec2[,i]!=0) > nrow(norm.rec2)*0.1) {
        
        TumorImp <- factor(norm.rec2$TumorImp, levels = c("Pre", "Post"))
        
        myLm <- lm(norm.rec2[,i] ~ TumorImp)
        myAnova <- anova(myLm)
        
        s_TumorImp[index] <- myLm$coefficients[2]
        p_TumorImp[index] <- myAnova$"Pr(>F)"[1]
        
        bugName[index] <- names(norm.rec2)[i]
        
        if(mean(norm.rec2[norm.rec2$TumorImp=="Pre",i])>mean(norm.rec2[norm.rec2$TumorImp=="Post",i])){
          TumorImp_direction[index]="Higher in Pre-tumor"
        }
        else if(mean(norm.rec2[norm.rec2$TumorImp=="Pre",i])<mean(norm.rec2[norm.rec2$TumorImp=="Post",i])){
          TumorImp_direction[index]="Lower in Pre-tumor"
        }
        else{
          TumorImp_direction[index]="Same"
        }
        
        index <- index + 1
        
      } # if (sum(norm.rec2[,i]!=0) > nrow(norm.rec2)*0.1)
      
    } # for (x in startTaxaIndex:ncol(norm.rec2))
    
    dFrame <- data.frame(bugName, s_TumorImp, p_TumorImp,TumorImp_direction)
    # Performing BH adjustments on p values
    dFrame$pAdj_TumorImp<- p.adjust( dFrame$p_TumorImp, method = "BH")
    dFrame <- parseQIIME2Taxa(dFrame)
    write.table(dFrame, paste0(outputLevel, level, "_", study, "_", Smoke.label, "_", term, "_LinearModelResults.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    rownames(dFrame) <- dFrame[,count]
    
    terms <- c("TumorImp")
    for (term in terms) {
      
      unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term))]
      
      boxPlotOut <- paste0(outputLevel, level, "_", study, "_", Smoke.label, "_", term, "_boxplots.pdf")
      pdf(boxPlotOut, width = 5, height = 5)
      for (j in 1:nrow(dFrame)) {
        pAdj_Values <- dFrame[,which(colnames(dFrame) == paste0("pAdj_", term))]
        if (pAdj_Values[j]<0.05) {

          pval <- roundP(pAdj_Values[j])
          p.label <- paste0("(adjusted ", pval, ")")
          
          bugName <- rownames(dFrame)[j]
          bugFull <- dFrame$bugName[j]
          NormAbundance <- norm.rec2[,which(colnames(norm.rec2) == bugFull)]
          Term <- norm.rec2[,which(colnames(norm.rec2) == term)]
          Term <- factor(Term, levels = c("Pre", "Post"))
          myFrame <- data.frame(NormAbundance, Term)
          title.lab <- paste(level, "-", bugName, 
                             "\n", Smoke.label, term, p.label, sep = " ")
          
          par(cex.axis=0.5, cex.main = 0.75)
          boxplot(NormAbundance ~ Term, main = title.lab, cex.main = 0.75,
                  xlab = term, ylab = "Normalized Abundance")
          stripchart(NormAbundance ~ Term,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
        } # if (pAdj_Values[j]<0.05)
      } # for (j in 1:nrow(dFrame))
      hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term, " p values"))
      dev.off()
      
    } # for (term in terms)
    
    
  } # for (i in Smoke.status)
  
  for (time in c("Pre", "Post")) {
    
    norm.FMTstatus <- norm[!(norm$TumorImp %in% time),]
    MouseType <- paste(norm.FMTstatus$FMTstatus, norm.FMTstatus$SmokeExp, sep = "\n")
    MouseType <- gsub(pattern = "Yes", replacement = "Smoke-exposed", MouseType)
    MouseType <- gsub(pattern = "No", replacement = "Smoke-free", MouseType)
    norm.FMTstatus <- cbind(MouseType, norm.FMTstatus)
    norm.FMTstatus$MouseType <- factor(norm.FMTstatus$MouseType, levels = c("Donor\nSmoke-exposed", "Recipient\nSmoke-exposed", "Recipient\nSmoke-free"))
    startTaxaIndex <- which(colnames(norm.FMTstatus)=="Study")+1
    
    term <- "MouseType"
    bugName <- vector()
    p_MouseType <- vector()
    p_Recipient_SmokeExpYes <- vector()
    p_Recipient_SmokeExpNo <- vector()
    s_Recipient_SmokeExpYes <- vector()
    s_Recipient_SmokeExpNo <- vector()
    index <- 1
    
    for (i in startTaxaIndex:ncol(norm.FMTstatus)) {
      
      if (sum(norm.FMTstatus[,i]!=0) > nrow(norm.FMTstatus)*0.1) {
        
        MouseType <- norm.FMTstatus$MouseType
        
        myLm <- lm(norm.FMTstatus[,i] ~ MouseType)
        myAnova <- anova(myLm)
        mySum <- summary(myLm)
        
        p_MouseType[index] <- myAnova$"Pr(>F)"[1]
        s_Recipient_SmokeExpYes[index] <- myLm$coefficients[2]
        s_Recipient_SmokeExpNo[index] <- myLm$coefficients[3]
        p_Recipient_SmokeExpYes[index] <- mySum$coefficients[2,4]
        p_Recipient_SmokeExpNo[index] <- mySum$coefficients[3,4]
        
        bugName[index] <- names(norm.FMTstatus)[i]
        
        index <- index + 1
        
      } # if (sum(norm.FMTstatus[,i]!=0) > nrow(norm.FMTstatus)*0.1)
      
    } # for (x in startTaxaIndex:ncol(norm.FMTstatus))
    
    dFrame <- data.frame(bugName, s_Recipient_SmokeExpYes, s_Recipient_SmokeExpNo, p_MouseType, p_Recipient_SmokeExpYes, p_Recipient_SmokeExpNo)
    # Performing BH adjustments on p values
    dFrame$pAdj_MouseType<- p.adjust( dFrame$p_MouseType, method = "BH")
    dFrame$pAdj_Recipient_SmokeExpYes<- p.adjust( dFrame$p_Recipient_SmokeExpYes, method = "BH")
    dFrame$pAdj_Recipient_SmokeExpNo<- p.adjust( dFrame$p_Recipient_SmokeExpNo, method = "BH")
    dFrame <- parseQIIME2Taxa(dFrame)
    write.table(dFrame, paste0(outputLevel, level, "_", study, "_", time, "-tumor", "_", term, "_LinearModelResults.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    rownames(dFrame) <- dFrame[,count]
    
    terms <- c("MouseType")
    for (term in terms) {
      
      unadj <- dFrame[,which(colnames(dFrame) == paste0("p_", term))]
      
      boxPlotOut <- paste0(outputLevel, level, "_", study, "_", term, "_", time, "-tumor_boxplots.pdf")
      pdf(boxPlotOut, width = 5, height = 5)
      for (j in 1:nrow(dFrame)) {
        pAdj_Values <- dFrame[,which(colnames(dFrame) == paste0("pAdj_", term))]
        if (pAdj_Values[j]<0.05) {
          
          pval <- roundP(pAdj_Values[j])
          p.label <- paste0("(adjusted ", pval, ")")
          
          bugName <- rownames(dFrame)[j]
          bugFull <- dFrame$bugName[j]
          NormAbundance <- norm.FMTstatus[,which(colnames(norm.FMTstatus) == bugFull)]
          Term <- norm.FMTstatus[,which(colnames(norm.FMTstatus) == term)]
          Term <- factor(Term, levels = c("Donor\nSmoke-exposed", "Recipient\nSmoke-exposed", "Recipient\nSmoke-free"))
          myFrame <- data.frame(NormAbundance, Term)
          title.lab <- paste(level, "-", bugName, 
                             "\n", time, "-tumor" , term, p.label, sep = " ")
          
          par(cex.axis=0.75, cex.main = 0.75)
          boxplot(NormAbundance ~ Term, main = title.lab, cex.main = 0.75,
                  xlab = term, ylab = "Normalized Abundance")
          stripchart(NormAbundance ~ Term,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
        } # if (pAdj_Values[j]<0.05)
      } # for (j in 1:nrow(dFrame))
      hist(unadj, main = paste0("Histogram of ", level, " level unadjusted p values"), cex.main = 0.75, xlab = paste0(term, " p values"))
      dev.off()
      
    } # for (term in terms)
    
  }
  
  
} # for (level in levels)

# Error


