#Author: Alicia Sorgen
#Date: 10-15-21
#Description: Performs basic linear modeling and produces simple boxplots

R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))
library(tidyr); message("tidyr:", packageVersion("tidyr"))

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Oct06"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = paste0("~/BioLockJ_pipelines/FMTmouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("CombineStudies"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  # studies = c("Abx-Bladder", "Abx-Colon", "Abx-Pancreatic", "Neo-Pancreatic")
  studies = c("NNK-FMT", "Smoke-FMT")
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "TableGather"),"/output/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
# source(funcScript)

levels = c("Phylum", "Class", "Order", "Family", "Genus")


for (level in levels) {
  
  dFrame <-data.frame()
  
  for (study in studies) {
    
    inputPath <- paste0(inputDir, study, "/")
    fileName <- paste0(level, "_table.tsv")
    fullPath <- paste0(inputPath, fileName)
    myT <- read.table(file = fullPath, sep="\t", header = TRUE)
    dFrame <- rbind(dFrame, myT)
    
  }
  
  outputFile <- paste0(output, level, "_CombineStudies_ALL.tsv")
  write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
  
  # dFrame <- dFrame[dFrame$AntibioticExp == "None",]
  # dFrame <- dFrame[dFrame$TumorImp == "Post",]
  # 
  # outputFile <- paste0(output, level, "_CombineStudies_NoAbx_PostTumor.tsv")
  # write.table(dFrame, outputFile, sep="\t", quote = FALSE, row.names = FALSE)
  
}
