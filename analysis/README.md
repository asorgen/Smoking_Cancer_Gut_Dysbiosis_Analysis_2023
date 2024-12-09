# Steps to Reproduce Smoking Cancer Gut Dysbiosis Analysis

## Using BioLockJ

### 1. Ensure that BioLockJ (v1.3.13 or newer) is on local system
BioLockJ - https://biolockj-dev-team.github.io/BioLockJ/Getting-Started/

### 2. Download Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 directory
git clone https://github.com/asorgen/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023.git

### 3. Set up required software

#### Option A) using Docker

Install docker.
Docker Desktop - https://www.docker.com/products/docker-desktop

Make sure the ` docker run hello-world ` command runs successfully.

The docker images required for this pipeline will be automatically pulled from the docker hub as need.  The first time the pipeline runs, startup will be slow as images are downloaded. 

**_If_** the specified images cannot be retrieved, they can be built from the docker files.  See the BioLockJ/Analysis/dockerfiles folder.  Build instructions are included in each file.

#### Option B) not using Docker

Make sure R is installed.  See https://www.r-project.org/.  These scripts were written with 4.0.2.

Make sure all required R packages are installed                                

- data.table
- ggplot2
- ggpubr
- gridExtra
- nlme
- rstatix
- scales
- stringr
- tidyr

### 4. Run BioLockJ pipeline

Move to the analysis folder:            
`cd <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023/analysis/BLJ_config_files`

To run the pipeline using **locally installed software**:                 
`biolockj mouseSmoke_analysis.properties`

To run the pipeline using **docker images**, add the -d argument:                                    
`biolockj -d mouseSmoke_analysis.properties`


## Running pipeline locally with R

### R & package versions used for this project

**R version 4.0.2 (2020-06-22)**

- **vegan**: Version 
- **stringr**: Version 
- **rstatix**: Version 
- **ggplot2**: Version 
- **ggpubr**: Version 
- **gridExtra**: Version 
- **data.table**: Version 


# Module Output

#### Figure 3A
- 07_LinearModeling_Abx-Pancreatic_CSEvCtrl > Genus > Figure3A_Genus_by_Treatment_Abx-Pancreatic_Post-tumor_LinearModel_plot.pdf
- 03_PCoA_Abx-Pancreatic_CSEvCtrl > Genus > Figure3A_Abx-Pancreatic_Genus_Post-tumor_Treatment_PCoA.pdf

#### Figure 3B
- 15_LinearModeling_Smoke-FMT_CasevCtrl > Genus > Figure3B_Genus_by_SmokeExposure_Smoke-FMT_Post-tumor_LinearModel_plot.pdf
- 12_PCoA_Smoke-FMT_CasevCtrl > Genus > Figure3B_Smoke-FMT_Genus_Post-tumor_SmokeExposure_PCoA.pdf



#### Supplemental Figure 3A 
- 03_PCoA_Abx-Pancreatic_CSEvCtrl > Genus > FigureS3A_Abx-Pancreatic_Genus_Pre-tumor_Treatment_PCoA.pdf

#### Supplemental Figure 3B 
- 07_LinearModeling_Abx-Pancreatic_CSEvCtrl > Genus > FigureS3B_Genus_by_Treatment_Abx-Pancreatic_Pre-tumor_LinearModel_plot.pdf

#### Supplemental Figure 3C
- 11_PCoA_NNK-FMT_CasevCtrl > Genus > FigureS3C_NNK-FMT_Genus_Pre-tumor_SmokeExposure_PCoA.pdf

#### Supplemental Figure 3D
- 16_LinearModeling_NNK-FMT_CasevCtrl > Genus > FigureS3D_Genus_by_SmokeExposure_NNK-FMT_Pre-tumor_LinearModel_plot.pdf

#### Supplemental Figure 3E
- 04_PCoA_Abx-Colon_CSEvCtrl > Genus > FigureS3E_Abx-Colon_Genus_Post-tumor_Treatment_PCoA.pdf

#### Supplemental Figure 3F
- 08_LinearModeling_Abx-Colon_CSEvCtrl > Genus > FigureS3F_Genus_by_Treatment_Abx-Colon_Post-tumor_LinearModel_plot.pdf

#### Supplemental Figure 3G
- 05_PCoA_Abx-Bladder_CSEvCtrl > Genus > FigureS3G_Abx-Bladder_Genus_Post-tumor_Treatment_PCoA.pdf

#### Supplemental Figure 3H
- 09_LinearModeling_Abx-Bladder_CSEvCtrl > Genus > FigureS3H_Genus_by_Treatment_Abx-Bladder_Post-tumor_LinearModel_plot.pdf

#### Supplemental Figure 4A
- 13_PCoA_Smoke-FMT_DonorRecipient-CasevCtrl > Genus > FigureS4A_Smoke-FMT_Genus_Pre-tumor_SmokeExposure_DonorRecipient_PCoA.pdf






### Run scripts in order by module

Move to the analysis folder:            
`cd <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023/analysis/Rscripts/`

#### 00. TaxaMetaMerge
`Rscript TaxaMetaMerge.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Bladder Abx-Colon Abx-Pancreatic Neo-Pancreatic NNK-FMT Smoke-FMT`

#### 01. CountConversion
`Rscript CountConversion.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Bladder Abx-Colon Abx-Pancreatic Neo-Pancreatic NNK-FMT Smoke-FMT`

#### 02. TableGather
`Rscript TableGather.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Bladder Abx-Colon Abx-Pancreatic Neo-Pancreatic NNK-FMT Smoke-FMT`

#### 03. PCoA_Abx-Pancreatic_CSEvCtrl
`Rscript PCoA_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Pancreatic`

#### 04. PCoA_Abx-Colon_CSEvCtrl
`Rscript PCoA_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Colon`

#### 05. PCoA_Abx-Bladder_CSEvCtrl
`Rscript PCoA_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Bladder`

#### 06. PCoA_Neo-Pancreatic_CSEvCtrl
`Rscript PCoA_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Neo-Pancreatic`

#### 07. LinearModeling_Abx-Pancreatic_CSEvCtrl
`Rscript LinearModeling_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Pancreatic`

#### 08. LinearModeling_Abx-Colon_CSEvCtrl
`Rscript LinearModeling_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Colon`

#### 09. LinearModeling_Abx-Bladder_CSEvCtrl
`Rscript LinearModeling_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Abx-Bladder`

#### 10. LinearModeling_Neo-Pancreatic_CSEvCtrl
`Rscript LinearModeling_CSEvCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Neo-Pancreatic`

#### 11. PCoA_NNK-FMT_CasevCtrl
`Rscript PCoA_CasevCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 NNK-FMT`

#### 12. PCoA_Smoke-FMT_CasevCtrl
`Rscript PCoA_CasevCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Smoke-FMT`

#### 13. PCoA_Smoke-FMT_DonorRecipient-CasevCtrl
`Rscript PCoA_DonorRecipient-CasevCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Smoke-FMT`

#### 14. PCoA_NNK-FMT_DonorRecipient-CasevCtrl
`Rscript PCoA_DonorRecipient-CasevCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 NNK-FMT`

#### 15. LinearModeling_Smoke-FMT_CasevCtrl
`Rscript LinearModeling_CasevCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Smoke-FMT`

#### 16. LinearModeling_NNK-FMT_CasevCtrl
`Rscript LinearModeling_CasevCtrl.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 NNK-FMT`

#### 17. LinearModeling_Smoke-FMT_DonorvRecipient
`Rscript LinearModeling_DonorvRecipient.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 Smoke-FMT`

#### 18. LinearModeling_NNK-FMT_DonorvRecipient
`Rscript LinearModeling_DonorvRecipient.R <path/to>/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 NNK-FMT`


