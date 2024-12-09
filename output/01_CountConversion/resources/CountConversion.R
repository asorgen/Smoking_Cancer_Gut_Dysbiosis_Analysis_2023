#Author: Alicia Sorgen
#Date: 08-18-21
#Description: Merge taxa count tables with metadata

##### Libraries #####
library(stringr); message("stringr:", packageVersion("stringr"))
library(data.table); message("data.table:", packageVersion("data.table"))

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "mouseSmoke"
params <- vector()
params <- c(params, "~/git/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023")
params <- c(params, "Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic")

moduleRoot <- paste0("CountConversion")

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
  
  # if (args[2] %in% c("Abx-Pancreatic", "Abx-Bladder", "Abx-Colon", "Neo-Pancreatic")) {
  #   mod1 <- sapply(strsplit(module, "_"), "[", 1)
  #   mod2 <- sapply(strsplit(module, "_"), "[", 2)
  #   module <- paste0(mod1, "_", args[2], "_", mod2)
  #   rm(mod1, mod2)
  # }
  
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
levels=c(levels, "Class")
levels=c(levels, "Order")
levels=c(levels, "Family")
levels=c(levels, "Genus")
levels=c(levels, "Species")

studies <- args[2:length(args)]

##### Set up input #####
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "TaxaMetaMerge"),"/output/")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Count Conversion #####
for (study in studies) {

  outputStudy = paste0(outputDir, study, "/")
  dir.create(outputStudy, showWarnings = FALSE)

  inputStudy = paste0(inputDir, "/", study, "/")

  levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")

  for (level in levels) {

    metaFile = paste0(level, "_metaCounts_", study, ".tsv")
    
    myT <- read.delim(paste0(inputStudy, metaFile), sep="\t",header = TRUE)
    startAbundanceIndex <- which(colnames(myT)=="Study")+1
    no_meta=myT[,startAbundanceIndex:ncol(myT)]
    meta <- myT[,1:(startAbundanceIndex-1)]
    no_meta=no_meta[,!(colnames(no_meta) %like% "d__Eukaryota.")]
    
    list=colnames(no_meta)
    
    num=grep(x = list,pattern = "mbiguous")
    num=c(num, grep(x = list, pattern = "Unassigned"))
    num=c(num, grep(x = list, pattern = "unculture"))
    list_other=list[num]
    reduced=no_meta[,-c(num)]
      
    if (length(num)>1) {
      reduced$Other=rowSums(no_meta[,num])
    } else {
      reduced$Other=no_meta[,num]
    }
    
    dFrame <- cbind(meta, reduced)
    write.table(dFrame, paste0(outputStudy, level, "_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    
    startTaxaIndex <- which(colnames(dFrame)=="Study")+1
    
    rel_abun=dFrame
    for (x in 1:nrow(dFrame)){
      rel_abun[x,startTaxaIndex:ncol(dFrame)]=((dFrame[x,startTaxaIndex:ncol(dFrame)])/(rowSums(dFrame[,startTaxaIndex:ncol(dFrame)])[x]))
    }
    write.table(rel_abun, paste0(outputStudy, level, "_relCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    
    Normalized=dFrame
    for (i in 1:nrow(dFrame)){
      Normalized[i,startTaxaIndex:ncol(dFrame)] = 
        log10(((dFrame[i,startTaxaIndex:ncol(dFrame)]) / (rowSums(dFrame[,startTaxaIndex:ncol(dFrame)])[i]) * 
                 rowMeans(dFrame[,startTaxaIndex:ncol(dFrame)])[i]) + 1)
    }
    write.table(Normalized, paste0(outputStudy, level, "_normCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    
      
  } # for (level in levels)

} # for (study in studies)

# Error
