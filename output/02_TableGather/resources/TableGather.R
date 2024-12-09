#Author: Alicia Sorgen
#Date: 08-19-21
#Description: gather all metadata, and taxa abundances

## Libraries
library(stringr)
library(tidyr)

rm(list=ls())

date = Sys.Date()
date = format(date, "%Y%b%d")
# date = "2021Aug16"

root = paste0("~/BioLockJ_pipelines/mouseSmoke_analysis_", date)
root = paste0("~/BioLockJ_pipelines/FMTmouseSmoke_analysis_", date)
root = dir(root, pattern=paste0("TableGather"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  studies = commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
  studies = c("NNK-FMT", "Smoke-FMT")
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/")
moduleDir = dirname(getwd())
output = file.path(moduleDir,"output/")

for (study in studies) {
  
  outputStudy = paste0(output, study, "/")
  dir.create(outputStudy, showWarnings = FALSE)
  
  inputStudy = paste0(inputDir, "/", study, "/")
  
  levels = c("Phylum", "Class", "Order", "Family", "Genus")
  count <- 1
  
  for (level in levels) {
    
    count <- count + 1
    inputFile <- paste0(inputDir, study, "/", level, "_rawCounts.tsv")
    raw <- read.table(file = inputFile,sep="\t",header = TRUE)
    
    inputFile <- paste0(inputDir, study, "/", level, "_relCounts.tsv")
    rel <- read.table(file = inputFile, sep="\t", header = TRUE)
    
    inputFile <- paste0(inputDir, study, "/", level, "_normCounts.tsv")
    norm <- read.table(file = inputFile, sep="\t", header = TRUE)
    
    startTaxaIndex <- which(colnames(raw)=="Study")+1
    
    df <- gather(norm, "Taxon","NormalizedCounts", startTaxaIndex:ncol(norm))
    df$Taxon=gsub(pattern =  "d__", x = df$Taxon,replacement = "")
    df$Taxon=gsub(pattern =  "[.].__", x = df$Taxon,replacement = "_/_")
    other <- paste(replicate(count,"Other"), collapse = "_/_")
    df$Taxon=gsub(pattern ="^Other$",x = df$Taxon,replacement = other)
    string <- strsplit(as.character(df$Taxon),split = "_/_")

    temp_string=do.call(rbind,string)
    df <- cbind(temp_string,df)
    
    colnames(df)[colnames(df)=="1"] <- "Domain"
    colnames(df)[colnames(df)=="2"] <- "Phylum"
    colnames(df)[colnames(df)=="3"] <- "Class"
    colnames(df)[colnames(df)=="4"] <- "Order"
    colnames(df)[colnames(df)=="5"] <- "Family"
    colnames(df)[colnames(df)=="6"] <- "Genus"
    
    rel2 <- gather(rel, "Taxon","Abundance", startTaxaIndex:ncol(rel))
    RelativeAbundance <- rel2$Abundance
    df <- cbind(df,RelativeAbundance)
    df$PercentAbundance <- df$RelativeAbundance * 100
    
    raw2 <- gather(raw, "Taxon","Raw_Counts", startTaxaIndex:ncol(raw))
    RawCounts <- raw2$Raw_Counts
    df <- cbind(df,RawCounts)
    
    write.table(df, paste0(outputStudy, level, "_table.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
    
  } # for (level in levels)
  
} # for (study in studies)

# Error