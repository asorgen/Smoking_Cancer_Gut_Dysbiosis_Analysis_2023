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
root = dir(root, pattern=paste0("PCoA_Abx-Pancreatic"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic", "NNK-FMT", "Smoke-FMT")
}

study <- "Neo-Pancreatic"

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
  
  # tumorStatus <- c("Pre", "Post")
  
  # for (t in tumorStatus) {
    
    # myT2 <- myT[myT$TumorImp == t,]
    myT2 <- myT
    myT2$Treatment <- factor(myT2$Treatment, levels = c("Control", "Smoke only", "Neomycin only", "Smoke+Neomycin"))
    taxaStart <- which(colnames(myT2)=="Study")+1
    
    myT2[, taxaStart:ncol(myT2)][is.na(myT2[, taxaStart:ncol(myT2)])] <- 0
    
    myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("steelblue", "tomato", "pink", "purple")
    col2=adjustcolor(color[factor(myT2$Treatment)], alpha.f = 1)
    
    adon.results<-adonis(myT2[, taxaStart:ncol(myT2)] ~ myT2$Treatment, method="bray",perm=999)
    # print(adon.results)
    capture.output(adon.results, file = paste0(output, study, "_", level, "_adonis_Treatment.txt"))
    Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
    
    pdf(paste(output, study, "_", level, "_Treatment_PCoA.pdf",sep = ""))
    par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
    
    pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:4) {
      ordiellipse(pcoa12, myT2$Treatment, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = color[n], show.groups = levels(factor(myT2$Treatment))[n], label = T, 
                  font = 2, cex = 1)
    }
    # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
    #        col = color, cex = 1.5, pch = 16, bty = "n")
    dev.off()
    
  # }
  
}
