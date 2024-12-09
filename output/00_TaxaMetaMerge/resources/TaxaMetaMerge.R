#Author: Alicia Sorgen
#Date: 07-13-21
#Description: Merge taxa count tables with metadata

## Libraries
library(stringr)

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Aug13"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = paste0("~/BioLockJ_pipelines/FMTmouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("TaxaMetaMerge"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("NNK-FMT", "Smoke-FMT")
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = file.path(pipeRoot, "input/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")

for (study in studies) {
  
  outputStudy = paste0(output, study, "/")
  dir.create(outputStudy, showWarnings = FALSE)

  inputStudy = paste0(inputDir, study, "/")
  metaFile = paste0(study, "_metadata.txt")
  
  metaTable <- read.delim(paste0(inputStudy, metaFile), sep="\t",header = TRUE)
  
  levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  x = 1
  
  for (level in levels) {
    
    x = x + 1
    countFile = paste0("level-", x, ".csv")
    countTable = read.delim(paste0(inputStudy, countFile), sep=",",header = TRUE)
    metaStart = which(colnames(countTable) == "MouseModel")
    names(countTable)[names(countTable) == "index"] <- "SampleID"
    countTable = countTable[,-(metaStart:ncol(countTable))]
    
    merged = merge(metaTable, countTable, by = "SampleID")
    write.table(merged, paste0(outputStudy,level,"_metaCounts_", study, ".tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    
  }
  
}

#Error