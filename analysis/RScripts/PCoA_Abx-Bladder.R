#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 09-02-21
#Description: 


## .) Set libraries
library(vegan)
library(stringr)

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Aug16"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("PCoA_Abx-Bladder"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic", "NNK-FMT", "Smoke-FMT")
}

study <- "Abx-Bladder"

## .) Set directories
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/", study, "/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")


## .) Create PCoAs for each taxonomic level
levels=c("Phylum","Class","Order","Family","Genus")

for (level in levels) {
  # outputDir <- paste0(output, level)
  
  myT=read.table(paste(inputDir,level,"_normCounts.tsv",sep=""),sep='\t',header = TRUE)
  myT$Treatment <- factor(myT$Treatment, levels = c("Control", "Smoke only", "Antibiotic only", "Smoke+Antibiotic"))
  taxaStart <- which(colnames(myT)=="Study")+1
  
  myT[, taxaStart:ncol(myT)][is.na(myT[, taxaStart:ncol(myT)])] <- 0
  
  myMDS <- capscale(myT[,taxaStart:(ncol(myT))]~1,distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  color <- c("steelblue", "tomato", "pink", "purple")
  col2=adjustcolor(color[factor(myT$Treatment)], alpha.f = 1)
  
  adon.results<-adonis(myT[, taxaStart:ncol(myT)] ~ myT$Treatment, method="bray",perm=999)
  # print(adon.results)
  capture.output(adon.results, file = paste0(output, study, "_", level, "_adonis_Treatment.txt"))
  Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
  
  pdf(paste(output, study, "_", level, "_Treatment_PCoA.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:4) {
    ordiellipse(pcoa12, myT$Treatment, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$Treatment))[n], label = T, 
                font = 2, cex = 1)
  }
  # legend("bottomright", c("Control", "Smoke Only", "Antibiotic Only", "Smoke + Antibiotic"),
  #        col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  
}
