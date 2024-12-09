# Smoking_Cancer_Gut_Dysbiosis_Analysis_2023 16S rRNA QIIME2 Sequence Processing

The folders `Dudeja_002_V4_analysis_qiime`, `Dudeja_004_V4_analysis`, `NNK_FMT_16S`, and `UMGC_Neo_16S` contain the sample manifest files with the SampleIDs and absolute file paths to the forward and reverse fastq files.

The `QIIME2` folder contains directories for the sub-studies within the analysis, each with individual metadata tables.

`slurmScripts/` contains the full pipeline and individual step-by-step scripts to perform the processing.
- Be sure to edit file paths for each script to reflect the correct file paths on your device.

This processing was originally performed on the UNCC HPC (Red Hat Enterprise Linux release 8.3 (Ootpa)) using the Slurm Workload Manager.

QIIME2 Version: 2019.1

## Processing pre-requisites

1. The manifest must reflect the absolute filepath for the sequences on your device.
2. A reference sequenece and taxonomy database must be downloaded. Pre-formatted reference sequence and taxonomy files (https://docs.qiime2.org/2021.2/data-resources/).

### All .qzv files can be downloaded and viewed at https://view.qiime2.org/

## To fully reproduce this processing pipeline:

Run slurmScripts in the following order:

1. `01_import.slurm`
2. `02_qzaSummary.slurm`
3. `03_D2_denoise.slurm`, `03_D4_denoise.slurm`, `03_NNK_denoise.slurm`, `03_UMGC_denoise.slurm` - the run order of these 4 scripts is not important.
4. `04_studyPartitioning.slurm`
5. `05_op[enRef_OTU.slurm`
6. `06_OTUfilter.slurm`
7. `07_readFilter.slurm`
8. `08_diversity.slurm`
9. `09_taxonomy.slurm`

The final outputs for this pipeline can be found in each of the `QIIME2` subdirectories under `visualizations/*_OR-taxonomy-barplot.qzv`
- Taxonomy count tables can be downloaded from these .qzv files  using https://view.qiime2.org/.
