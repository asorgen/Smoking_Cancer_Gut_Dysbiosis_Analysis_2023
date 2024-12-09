#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 02-07-23
#Description: Generate PCoA plots for pre- and post-tumor Case v. Control FMT mice


##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(vegan); message("vegan: Version", packageVersion("vegan"))
library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "mouseSmoke"
params <- vector()
params <- c(params, "~/git/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023")
params <- c(params, "NNK-FMT")

moduleRoot <- paste0("PCoA_CasevCtrl")

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "data")
  gitScripts <- file.path(gitRoot, "analysis", "RScripts")
  rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  rm(pipeline)
  
  if (any(dir(root) == "input") == FALSE) {
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  rm(gitInput)
  module <- moduleRoot
  
  if (args[2] %in% c("NNK-FMT", "Smoke-FMT")) {
    mod1 <- sapply(strsplit(module, "_"), "[", 1)
    mod2 <- sapply(strsplit(module, "_"), "[", 2)
    module <- paste0(mod1, "_", args[2], "_", mod2)
    rm(mod1, mod2)
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
  }
  rm(scriptDir)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  rm(outputDir, files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    rm(script)
    
  }
  rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))
  rm(moduleDir)
  
}
rm(params, moduleRoot, ANALYSIS)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

##### Set up functions file #####
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str, funcScript)

##### Script variables #####
levels="Phylum"
# levels=c(levels, "Class")
# levels=c(levels, "Order")
# levels=c(levels, "Family")
levels=c(levels, "Genus")
levels=c(levels, "Species")

study <- args[2]
CSE <- sapply(strsplit(study, "-"), "[", 1)

##### Set up input #####
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/", study, "/")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### PCoA Antibiotic-free Smoke v. Control (Pre- and Post-tumor) #####
for (level in levels) {
  outputLevel <- paste0(outputDir, level)
  dir.create(outputLevel, showWarnings = FALSE)
  
  myT=read.table(paste(inputDir,level,"_normCounts.tsv",sep=""),sep='\t',header = TRUE)
  myT <- myT[myT$FMTstatus == "Recipient",]
  
  tumorStatus <- unique(myT$TumorStatus)
  
  for (t in tumorStatus) {
    
    myT2 <- myT[myT$TumorStatus == t,]
    myT2$SmokeExposure <- gsub("-exposed feces", "", myT2$SmokeExposure)
    myT2$SmokeExposure <- gsub(paste0(CSE, "-free feces"), "Control", myT2$SmokeExposure)
    myT2$SmokeExposure <- factor(myT2$SmokeExposure, levels = c(CSE, "Control"))
    taxaStart <- which(colnames(myT2)=="Study")+1
    
    myT2[, taxaStart:ncol(myT2)][is.na(myT2[, taxaStart:ncol(myT2)])] <- 0
    
    myMDS <- capscale(myT2[,taxaStart:(ncol(myT2))]~1,distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    color <- c("forestgreen", "red")
    # color <- c("steelblue", "tomato", "pink", "purple", "orange", "green", "grey", "maroon")
    col2=adjustcolor(color[factor(myT2$SmokeExposure)], alpha.f = 1)
    
    adon.results<-adonis2(myT2[, taxaStart:ncol(myT2)] ~ myT2$SmokeExposure, method="bray",perm=999)
    # print(adon.results)
    capture.output(adon.results, file = paste0(outputLevel, "/", study, "_", level, "_", t, "_SmokeExposure_adonis.txt"))
    # Title <- paste0(level, " Bray-Curtis dissimilarity\nPERMANOVA p = ", adon.results$`Pr(>F)`[1])
    # Title <- paste0(CSE, " FMT: KPC ", t, "\nPERMANOVA p = ", adon.results$`Pr(>F)`[1])
    Title <- paste0(CSE, " exposure\nPERMANOVA p = ", adon.results$`Pr(>F)`[1])
    
    if (level == "Genus" & study == "Smoke-FMT" & t == "Post-tumor") {
      pdf(paste(outputLevel, "/Figure3B_", study, "_", level, "_", t, "_SmokeExposure_PCoA.pdf",sep = ""))
    } else if (level == "Genus" & study == "NNK-FMT" & t == "Pre-tumor") {
      pdf(paste(outputLevel, "/FigureS3C_", study, "_", level, "_", t, "_SmokeExposure_PCoA.pdf",sep = ""))
    } else {
      pdf(paste(outputLevel, "/", study, "_", level, "_", t, "_SmokeExposure_PCoA.pdf",sep = ""))
    }
    par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(1,1))
    
    all.combindations <- combn(1:5, 2)
    for (combination in 1:ncol(all.combindations)) {
      PCoA_a <- all.combindations[1,combination]
      PCoA_b <- all.combindations[2,combination]
      
      pcoa12 <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2,
                         xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = Title)
      points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
      
      for (n in 1:length(unique(na.omit(myT2$SmokeExposure)))) {
        ordiellipse(pcoa12, myT2$SmokeExposure, display = "sites", 
                    kind = "se", conf = 0.95, 
                    lwd = 2, # line width
                    lty = 2, # line type
                    # draw = "polygon",
                    draw = "line",
                    col = color[n], 
                    show.groups = levels(factor(myT2$SmokeExposure))[n], 
                    label = T,
                    alpha = 50, # transparency of ellipse fill
                    font = 2, 
                    cex = 1)
      }
      # legend("topright", unique(myT2$SmokeExposure),
      #        col = color[1:length(unique(na.omit(myT2$SmokeExposure)))], cex = 1.5, pch = 16, bty = "n")
      
    }
    
    
    dev.off()
    
  }
  
} # for (level in levels)
