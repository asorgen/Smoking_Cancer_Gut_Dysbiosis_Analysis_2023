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
# params <- c(params, "Smoke-FMT")

moduleRoot <- paste0("LinearModeling_CasevCtrl")

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

##### Linear Regression Modeling for Smoke v. Control (Pre- and Post-tumor) #####
for (level in levels) {
  outputLevel <- paste0(outputDir, level, "/")
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
    
    # Parse out metadata from counts
    taxa.df <- myT2[,taxaStart:ncol(myT2)]
    
    SmokeExposure <- myT2$SmokeExposure
    
    Taxa <- vector()
    pValue <- vector()
    R2 <- vector()
    direction <- vector()
    index <- 1
    
    for (i in 1:ncol(taxa.df)) {
      
      bug <- taxa.df[,i]
      bugName <- colnames(taxa.df)[i]
      
      if ( mean(bug > 0, na.rm = TRUE) > 0.1 ) {
        
        df <- data.frame( bug, SmokeExposure )
        
        myLm <- anova(lm( bug ~ SmokeExposure, data = df ))
        lmSum <- summary(lm( bug ~ SmokeExposure, data = df ))
        
        Taxa[index] <- bugName
        pValue[index] <- myLm$`Pr(>F)`[1]
        R2[index] <- lmSum$r.squared
        
        averages <- df %>%
          group_by(SmokeExposure) %>%
          get_summary_stats(bug, type = "mean")
        
        Ctrl <- averages$mean[which(averages$SmokeExposure == "Control")]
        Case <- averages$mean[which(averages$SmokeExposure == CSE)]
        
        direction[index] <- ifelse(Ctrl > Case, "Control",
                                   ifelse(Ctrl < Case, CSE,
                                          "same"))
        
        index <- index + 1
        
        
      } # if ( mean(bug > 0, na.rm = TRUE) > 0.1 )
      
    } # for (i in 1:ncol(taxa.df))
    
    dFrame <- data.frame( Taxa, pValue , R2, direction )
    dFrame$Adj_pValue <- p.adjust(dFrame$pValue, method = "BH")
    
    file.path <- paste0(outputLevel, level, "_by_SmokeExposure_", study, "_", t, "_LinearModelResults.tsv")
    write.table(dFrame, file=file.path, sep="\t", row.names=FALSE)
    
    file.path <- paste0(outputLevel, level,  "_by_SmokeExposure_", study, "_", t, "_pValue_Histogram.pdf")
    pdf(file.path)
    
    pValues <- dFrame[,which(colnames(dFrame) == "pValue")]
    n=length(pValues)
    title.lab <- paste0( "Linear Regression\n", t, " FMT Control v. ", CSE, " p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
    
    dev.off()
    
    

  } # for (t in tumorStatus)
  
} # for (level in levels)


##### Plotting significant results #####
for (level in levels) {
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  # tumorStatus <- unique(myT$TumorStatus)
  
  for (t in tumorStatus) {
    
    lmFileName <- paste0(outputLevel, level, "_by_SmokeExposure_", study, "_", t, "_LinearModelResults.tsv")
    dFrame = read.table(lmFileName,sep='\t',header = TRUE)
    
    dFrame2 <- dFrame[dFrame$Adj_pValue < 0.05,]
    dFrame2$logP <- ifelse(dFrame2$direction == CSE, -log10(dFrame2$Adj_pValue),
                           log(dFrame2$Adj_pValue))
    
    min <- ifelse(min(dFrame2$logP) < -32, floor(min(dFrame2$logP)), -32)
    max <- ifelse(max(dFrame2$logP) > 30, ceiling(max(dFrame2$logP)), 30)
    breakValues <- seq(min, max, 2)
    
    
    dFrame2$Taxa <- gsub("d__Bacteria.", "", dFrame2$Taxa)
    # TaxaNames <- strsplit(dFrame2$Taxa, "__")
    # TaxaNames <- do.call(rbind, TaxaNames)
    # dFrame2$TaxaNames <- TaxaNames[,ncol(TaxaNames)]
    
    dFrame2 <- dFrame2 %>%
      mutate(
        Taxa = factor(Taxa, levels = Taxa[order(logP, decreasing = FALSE)]),
        label_y = ifelse(logP < 0, 0.2, -0.2),
        label_hjust = ifelse(logP < 0, 0, 1)
      )
    
    plot <- ggplot(dFrame2, aes(x = Taxa, y = logP, fill = direction)) +
      geom_bar(stat = "identity", col = "black")+
      geom_text(aes(y = label_y, label = Taxa, hjust = label_hjust), size = 3)+
      coord_flip()+
      scale_fill_manual(values = c("red","forestgreen"))+
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 15),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 15),
            # legend.justification = 0,
            legend.title = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
            panel.grid.minor.x = element_blank())+
      scale_y_continuous(expression(log[10]("p value")),
                         breaks = c(breakValues), limits = c(min, max)); plot
    
    # file.path <- paste0(outputLevel, level, "_by_SmokeExposure_", study, "_", t, "-tumor_LinearModelResults.tsv")
    # write.table(dFrame, file=file.path, sep="\t", row.names=FALSE)
    if (level == "Genus" & study == "Smoke-FMT" & t == "Post-tumor") {
      file.path <- paste0(outputLevel, "Figure3B_", level,  "_by_SmokeExposure_", study, "_", t, "_LinearModel_plot.pdf")
      pdf(file.path, width = 16, height = 7)
    } else if (level == "Genus" & study == "NNK-FMT" & t == "Pre-tumor") {
      file.path <- paste0(outputLevel, "FigureS3D_", level,  "_by_SmokeExposure_", study, "_", t, "_LinearModel_plot.pdf")
      pdf(file.path, width = 17, height = 8)
    } else {
      file.path <- paste0(outputLevel, level,  "_by_SmokeExposure_", study, "_", t, "_LinearModel_plot.pdf")
      pdf(file.path, width = 14, height = 7)
    }
    print(plot)
    dev.off()
    
  } # for (t in tumorStatus)
  
} # for (level in levels)
