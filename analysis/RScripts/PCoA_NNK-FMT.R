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
root = dir(root, pattern=paste0("PCoA_NNK-FMT"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic", "NNK-FMT", "Smoke-FMT")
}

study <- "NNK-FMT"

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
  
  tumorStatus <- c("Pre", "Post")
  
  # for (t in tumorStatus) {
  
  # myT2 <- myT[myT$TumorImp == t,]
  myT2 <- myT
  Mouse <- paste0(myT2$NNKexp, myT2$FMTstatus, myT$TumorImp)
  Mouse <- gsub(pattern ="NoRecipientPost",x = Mouse, replacement = "Post-tumor NNK-free recipient")
  Mouse <- gsub(pattern ="YesRecipientPost",x = Mouse, replacement = "Post-tumor NNK-exposed recipient")
  Mouse <- gsub(pattern ="YesRecipientPre",x = Mouse, replacement = "Pre-tumor NNK-exposed recipient")
  Mouse <- gsub(pattern ="NoRecipientPre",x = Mouse, replacement = "Pre-tumor NNK-free recipient")
  Mouse <- gsub(pattern ="YesDonorNA",x = Mouse, replacement = "NNK-exposed donor")
  myT2 <- cbind(Mouse, myT2)
  myT2$Mouse <- factor(myT2$Mouse, levels = c("NNK-exposed donor", "Pre-tumor NNK-free recipient", "Pre-tumor NNK-exposed recipient", "Post-tumor NNK-free recipient", "Post-tumor NNK-exposed recipient"))
  taxaStart <- which(colnames(myT2)=="Study")+1
  
  myT2[, taxaStart:ncol(myT2)][is.na(myT2[, taxaStart:ncol(myT2)])] <- 0
  
  myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  color <- c("steelblue", "tomato", "pink", "purple", "orange")
  col2=adjustcolor(color[factor(myT2$Mouse)], alpha.f = 1)
  
  adon.results<-adonis(myT2[, taxaStart:ncol(myT2)] ~ myT2$Mouse, method="bray",perm=999)
  # print(adon.results)
  capture.output(adon.results, file = paste0(output, study, "_", level, "_adonis_Mouse.txt"))
  Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
  
  pdf(paste(output, study, "_", level, "_Mouse_PCoA.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:5) {
    ordiellipse(pcoa12, myT2$Mouse, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT2$Mouse))[n], label = T, 
                font = 2, cex = 1)
  }
  # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
  #        col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  # }
  
}
