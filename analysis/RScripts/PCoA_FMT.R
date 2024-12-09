#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 09-02-21
#Description: 


## .) Set libraries
library(vegan); message("vegan:", packageVersion("vegan"))
library(stringr); message("stringr:", packageVersion("stringr"))
library(data.table); message("data.table:", packageVersion("data.table"))

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Sep02"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = paste0("~/BioLockJ_pipelines/FMTmouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("PCoA_FMT"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("NNK-FMT")
}

# study <- "Smoke-FMT"

for (study in studies) {
  
  message("\n***** Starting ", study, " *****\n")
  
  ## .) Set directories
  pipeRoot = dirname(dirname(getwd()))
  inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/", study, "/"); message(inputDir)
  moduleDir = dirname(getwd()); message(moduleDir)
  output = file.path(moduleDir,"output/")
  output <- paste0(output, study, "/"); message(output)
  dir.create(output, showWarnings = FALSE)
  
  ## .) Create PCoAs for each taxonomic level
  levels=c("Phylum","Class","Order","Family","Genus")
  
  for (level in levels) {
    
    message("\n***** Starting ", level, " *****\n")
    
    outputDir <- paste0(output, level, "/"); message(outputDir)
    dir.create(outputDir, showWarnings = FALSE)
    
    myT=read.table(paste(inputDir,level,"_normCounts.tsv",sep=""),sep='\t',header = TRUE)
    
    taxaStart <- which(colnames(myT)=="Study")+1
    
    myT[, taxaStart:ncol(myT)][is.na(myT[, taxaStart:ncol(myT)])] <- 0
    
    message("PCoA")
    
    # PCoA showing all cohorts
    myMDS <- capscale(myT[,taxaStart:(ncol(myT))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("steelblue", "tomato", "pink", "purple", "orange", "green", "gray")
    col2=adjustcolor(color[factor(myT$FullMouseDescription)], alpha.f = 1)
    num <- length(unique(myT$FullMouseDescription))
    
    message("adonis")
    adon.results<-adonis(myT[, taxaStart:ncol(myT)] ~ myT$FullMouseDescription, method="bray",perm=999)
    capture.output(adon.results, file = paste0(outputDir, study, "_", level, "_adonis_Mouse.txt"))
    
    Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
    
    message("Figure")
    pdf(paste(outputDir, study, "_", level, "_Mouse_PCoA.pdf",sep = ""))
    par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
    
    pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:num) {
      ordiellipse(pcoa12, myT$FullMouseDescription, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = color[n], show.groups = levels(factor(myT$FullMouseDescription))[n], label = T, 
                  font = 2, cex = 1)
    }
    # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
    #        col = color, cex = 1.5, pch = 16, bty = "n")
    dev.off()
    
    
    
    # PCoA showing donor v. pre-tumor recipients
    name <- "Donor_v_Pre-tumor"
    myT2 <- myT[!(myT$TumorStatus == "Post-tumor"),]
    myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("steelblue", "tomato", "pink", "purple", "orange", "green", "gray")
    col2=adjustcolor(color[factor(myT2$FullMouseDescription)], alpha.f = 1)
    num <- length(unique(myT2$FullMouseDescription))
    
    adon.results<-adonis(myT2[, taxaStart:ncol(myT2)] ~ myT2$FullMouseDescription, method="bray",perm=999)
    capture.output(adon.results, file = paste0(outputDir, study, "_", level, "_adonis_", name, ".txt"))
    
    Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
    
    pdf(paste(outputDir, study, "_", level, "_", name, "_PCoA.pdf",sep = ""))
    par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
    
    pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:num) {
      ordiellipse(pcoa12, myT2$FullMouseDescription, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = color[n], show.groups = levels(factor(myT2$FullMouseDescription))[n], label = T, 
                  font = 2, cex = 1)
    }
    # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
    #        col = color, cex = 1.5, pch = 16, bty = "n")
    dev.off()
    
    
    
    
    # PCoA showing pre-tumor recipients v. post-tumor
    name <- "Pre-tumor_Post-tumor"
    myT2 <- myT[!(myT$FMTstatus == "Donor"),]
    myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("steelblue", "tomato", "pink", "purple", "orange", "green", "gray")
    col2=adjustcolor(color[factor(myT2$FullMouseDescription)], alpha.f = 1)
    num <- length(unique(myT2$FullMouseDescription))
    
    adon.results<-adonis(myT2[, taxaStart:ncol(myT2)] ~ myT2$FullMouseDescription, method="bray",perm=999)
    capture.output(adon.results, file = paste0(outputDir, study, "_", level, "_adonis_", name, ".txt"))
    
    Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
    
    pdf(paste(outputDir, study, "_", level, "_", name, "_PCoA.pdf",sep = ""))
    par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
    
    pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:num) {
      ordiellipse(pcoa12, myT2$FullMouseDescription, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = color[n], show.groups = levels(factor(myT2$FullMouseDescription))[n], label = T, 
                  font = 2, cex = 1)
    }
    # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
    #        col = color, cex = 1.5, pch = 16, bty = "n")
    dev.off()
    
    
    
    
    
    # PCoA showing Smoke-free mice
    nicType <- strsplit(study,"-")[[1]][1]    
    name <- paste0(nicType, "-free_Mice")
    myT2 <- myT[(myT$SmokeExposure %like% paste0(nicType, "-free")),]
    myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("steelblue", "tomato", "pink", "purple", "orange", "green", "gray")
    col2=adjustcolor(color[factor(myT2$FullMouseDescription)], alpha.f = 1)
    num <- length(unique(myT2$FullMouseDescription))
    
    adon.results<-adonis(myT2[, taxaStart:ncol(myT2)] ~ myT2$FullMouseDescription, method="bray",perm=999)
    capture.output(adon.results, file = paste0(outputDir, study, "_", level, "_adonis_", name, ".txt"))
    
    Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
    
    pdf(paste(outputDir, study, "_", level, "_", name, "_PCoA.pdf",sep = ""))
    par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
    
    pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:num) {
      ordiellipse(pcoa12, myT2$FullMouseDescription, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = color[n], show.groups = levels(factor(myT2$FullMouseDescription))[n], label = T, 
                  font = 2, cex = 1)
    }
    # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
    #        col = color, cex = 1.5, pch = 16, bty = "n")
    dev.off()
    

    
    
    
    # PCoA showing Smoke-exposed mice
    nicType <- strsplit(study,"-")[[1]][1]    
    name <- paste0(nicType, "-exposed_Mice")
    myT2 <- myT[(myT$SmokeExposure %like% paste0(nicType, "-exposed")),]
    myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("steelblue", "tomato", "pink", "purple", "orange", "green", "gray")
    col2=adjustcolor(color[factor(myT2$FullMouseDescription)], alpha.f = 1)
    num <- length(unique(myT2$FullMouseDescription))
    
    adon.results<-adonis(myT2[, taxaStart:ncol(myT2)] ~ myT2$FullMouseDescription, method="bray",perm=999)
    capture.output(adon.results, file = paste0(outputDir, study, "_", level, "_adonis_", name, ".txt"))
    
    Title <- paste0(level, " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
    
    pdf(paste(outputDir, study, "_", level, "_", name, "_PCoA.pdf",sep = ""))
    par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
    
    pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:num) {
      ordiellipse(pcoa12, myT2$FullMouseDescription, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = color[n], show.groups = levels(factor(myT2$FullMouseDescription))[n], label = T, 
                  font = 2, cex = 1)
    }
    # legend("bottomright", c("Control", "Smoke Only", "Neomycin Only", "Smoke + Neomycin"),
    #        col = color, cex = 1.5, pch = 16, bty = "n")
    dev.off()
    
    # }
    
  }
}

