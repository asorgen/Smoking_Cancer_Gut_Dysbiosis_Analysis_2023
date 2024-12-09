#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 05-15-23
#Description: Graph for the increase in Bacteroides to Firmicutes ratio with CSE


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
params <- c(params, "Abx-Pancreatic")

moduleRoot <- paste0("Bacteroides_Firmicutes_Ratio")

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
  
  if (args[2] %in% c("Abx-Pancreatic", "Abx-Bladder", "Abx-Colon", "Neo-Pancreatic")) {
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
level="Phylum"

study <- args[2]
cancer <- sapply(strsplit(study, "-"), "[", 2)

##### Set up input #####
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CountConversion"),"/output/", study, "/")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Linear Regression Modeling for Antibiotic-free Smoke v. Control (Pre- and Post-tumor) #####

myT=read.table(paste(inputDir,level,"_normCounts.tsv",sep=""),sep='\t',header = TRUE)

# t <- "Post"
# myT2 <- myT[myT$TumorImp == t,]
myT2 <- myT
myT2$Treatment <- factor(myT2$Treatment)
taxaStart <- which(colnames(myT2)=="Study")+1

myT2[, taxaStart:ncol(myT2)][is.na(myT2[, taxaStart:ncol(myT2)])] <- 0
myT2 <- myT2[myT2$AntibioticExp == "None",]

SampleID <- myT2$SampleID
Treatment <- myT2$Treatment
Treatment <- gsub("Smoke only", "CSE", Treatment)
Bacteroidota <- myT2$d__Bacteria.p__Bacteroidota
Firmicutes <- myT2$d__Bacteria.p__Firmicutes
Tumor <- paste0(myT2$TumorImp, "-tumor")
Tumor <- factor(Tumor, levels = c("Pre-tumor", "Post-tumor"))

df <- data.frame(SampleID, Tumor, Treatment, Bacteroidota, Firmicutes)
df$Ratio <- df$Bacteroidota / df$Firmicutes

df2 <- df %>%
  gather("Phylum", "Abundance", 4:5)

plot <- ggplot(df2, aes(fill = Phylum, x = Treatment, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Total Log Abundance") +
  facet_wrap(~Tumor); plot

file.path <- paste0(outputDir,  "Total_Bacteroidota_Firmicutes_by_Treatment_Tumor_", study, ".pdf")
pdf(file.path, width = 14, height = 7)
print(plot)
dev.off()

stat.t_test <- df2 %>%
  group_by(Tumor) %>%
  t_test(Ratio ~ Treatment) %>%
  add_xy_position(x = "Treatment")

stat.wilcox_test <- df2 %>%
  group_by(Tumor) %>%
  wilcox_test(Ratio ~ Treatment) %>%
  add_xy_position(x = "Treatment")

plot <- ggplot(df2, aes(x = Treatment, y = Ratio)) +
  geom_boxplot() +
  ylab("Ratio (Bacteroidota:Firmicutes)") +
  facet_wrap(~Tumor); plot

plot_t_test <- plot +
  stat_pvalue_manual(
    stat.t_test,
    bracket.nudge.y = 0.01,
    size = 4,
    color = "black",
    hide.ns = FALSE,
    label = "p",
    step.increase = 0.07) + labs(title = "t-test"); plot_t_test

plot_wilcox_test <- plot +
  stat_pvalue_manual(
    stat.wilcox_test,
    bracket.nudge.y = 0.01,
    size = 4,
    color = "black",
    hide.ns = FALSE,
    label = "p",
    step.increase = 0.07) + labs(title = "Wilcoxon"); plot_wilcox_test

file.path <- paste0(outputDir,  "Ratio_Bacteroidota_Firmicutes_by_Treatment_Tumor_", study, ".pdf")
pdf(file.path, width = 14, height = 7)
print(plot_t_test)
print(plot_wilcox_test)
dev.off()

# plot <- ggplot(df, aes(fill = Phylum, x = Tumor, y = Abundance)) +
#   geom_bar(position = "stack", stat = "identity") +
#   facet_wrap(~Treatment); plot


Ratio_Avg <- df %>% 
  group_by(Tumor, Treatment) %>%
  get_summary_stats(Ratio, type = "mean")
names(Ratio_Avg)[names(Ratio_Avg) == "mean"] <- "Ratio"

Bacteroidota_Avg <- df %>% 
  group_by(Tumor, Treatment) %>%
  get_summary_stats(Bacteroidota, type = "mean")
Bacteroidota <- Bacteroidota_Avg$mean

Firmicutes_Avg <- df %>% 
  group_by(Tumor, Treatment) %>%
  get_summary_stats(Firmicutes, type = "mean")
Firmicutes <- Firmicutes_Avg$mean

df3 <- cbind(Ratio_Avg, Bacteroidota)
df3 <- cbind(df3, Firmicutes)

df3 <- df3 %>%
  gather("Phylum", "Abundance", 6:7)

file.path <- paste0(outputDir,  "Mean_Bacteroidota_Firmicutes_Ratio_by_Treatment_Tumor_", study, ".tsv")
write.table(df3, file=file.path, sep="\t", row.names=FALSE)


plot <- ggplot(df3, aes(fill = Phylum, x = Treatment, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Mean Log Abundance") +
  facet_wrap(~Tumor); plot

file.path <- paste0(outputDir,  "Mean_Bacteroidota_Firmicutes_by_Treatment_Tumor_", study, ".pdf")
pdf(file.path, width = 14, height = 7)
print(plot)
dev.off()

