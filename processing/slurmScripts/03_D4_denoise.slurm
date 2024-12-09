#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=D4_denoise                                    # Job name
#SBATCH --time=10:00:00                                      # Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         # Number of nodes requested
#SBATCH --ntasks-per-node=16                                 # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb                                          # Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=asorgen@uncc.edu                            # Specified email address
#SBATCH --output=/users/asorgen/mouseSmoke/outputLogs/%x.%j.out      # Standard output
#SBATCH --error=/users/asorgen/mouseSmoke/outputLogs/%x.%j.err          # Error log

logDir=/users/asorgen/mouseSmoke/outputLogs
jobName=$SLURM_JOB_NAME
outFile=$SLURM_JOB_NAME.$SLURM_JOB_ID

# echo Output file: $outFile

mv ${logDir}/${jobName}.* ${logDir}/oldLogs
mv ${logDir}/oldLogs/${outFile}.* ${logDir}

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load openmpi
module load qiime2/2019.1


study=Dudeja_004_V4_analysis
# Total samples: 269
# Total sequence counts: 17339389

inputDirPath=/scratch/asorgen/16Sdatasets/mouseSmoke

mkdir ${inputDirPath}/${study}/tables

seqs=${inputDirPath}/${study}/seqs/${study}
visualizations=${inputDirPath}/${study}/visualizations/${study}
tables=${inputDirPath}/${study}/tables/${study}


# Denoise and quality trim sequences
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ${seqs}_raw-seqs.qza \
--p-trim-left-f 10 \
--p-trim-left-r 10 \
--p-trunc-len-f 220 \
--p-trunc-len-r 220 \
--o-table ${tables}_ASVtable.qza \
--o-representative-sequences ${seqs}_ASVseqs.qza \
--o-denoising-stats ${tables}_denoising-stats.qza \
--verbose \
--p-n-threads 16


# Produce readable ASV table file
qiime feature-table summarize \
--i-table ${tables}_ASVtable.qza \
--o-visualization ${visualizations}_ASVtable.qzv


# Produce readable ASV sequence file
qiime feature-table tabulate-seqs \
--i-data ${seqs}_ASVseqs.qza \
--o-visualization ${visualizations}_ASVseqs.qzv


# Remove ASVs with frequency <10
qiime feature-table filter-features \
--i-table ${tables}_ASVtable.qza \
--p-min-frequency 10 \
--o-filtered-table ${tables}_filtered-ASVtable.qza


# Produce readable filtered ASV table file
qiime feature-table summarize \
--i-table ${tables}_filtered-ASVtable.qza \
--o-visualization ${visualizations}_filtered-ASVtable.qzv


# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${seqs}_ASVseqs.qza \
--i-table ${tables}_filtered-ASVtable.qza \
--o-filtered-data ${seqs}_filtered-ASVseqs.qza
