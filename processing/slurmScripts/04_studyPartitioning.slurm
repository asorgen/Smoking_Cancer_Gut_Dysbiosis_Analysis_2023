#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=studies                                    # Job name
#SBATCH --time=10:00:00                                      # Time limit (hrs:min:sec)
#SBATCH --nodes=1                              	           # Number of nodes requested
#SBATCH --ntasks-per-node=16                                 # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb                                          # Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=asorgen@uncc.edu                            # Specified email address
#SBATCH --output=/users/asorgen/mouseSmoke/outputLogs/%x.%j.out   	    # Standard output
#SBATCH --error=/users/asorgen/mouseSmoke/outputLogs/%x.%j.err          # Error log

logDir=/users/asorgen/mouseSmoke/outputLogs
jobName=$SLURM_JOB_NAME
outFile=$SLURM_JOB_NAME.$SLURM_JOB_ID

echo Output file: $outFile

mv ${logDir}/${jobName}.* ${logDir}/oldLogs
mv ${logDir}/oldLogs/${outFile}.* ${logDir}

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load openmpi
module load qiime2/2019.1

# studies=(UMGC_Neo_16S NNK_FMT_16S Dudeja_004_V4_analysis Dudeja_002_V4_analysis_qiime)

inputDirPath=/scratch/asorgen/16Sdatasets/mouseSmoke
outputDirPath=/users/asorgen/mouseSmoke/QIIME2


study=UMGC_Neo_16S
EXP=Neo-Pancreatic

IN=${inputDirPath}/${study}/tables/${study}
META=${outputDirPath}/${EXP}/${EXP}_metadata.txt
tables=${outputDirPath}/${EXP}/tables
seqs=${outputDirPath}/${EXP}/seqs
oldSeqs=${inputDirPath}/${study}/seqs/${study}

rm -r ${tables}
rm -r ${seqs}

mkdir ${tables}
mkdir ${seqs}

# Filter out samples based on individual experiments
qiime feature-table filter-samples \
--i-table ${IN}_ASVtable.qza \
--m-metadata-file ${META} \
--p-where "[Study]='Neo-Pancreatic'" \
--o-filtered-table ${tables}/${EXP}_ASVtable.qza

# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${oldSeqs}_ASVseqs.qza \
--i-table ${tables}/${EXP}_ASVtable.qza \
--o-filtered-data ${seqs}/${EXP}_ASVseqs.qza


#######################################################

study=NNK_FMT_16S
EXP=NNK-FMT

IN=${inputDirPath}/${study}/tables/${study}
META=${outputDirPath}/${EXP}/${EXP}_metadata.txt
tables=${outputDirPath}/${EXP}/tables
seqs=${outputDirPath}/${EXP}/seqs
oldSeqs=${inputDirPath}/${study}/seqs/${study}

rm -r ${tables}
rm -r ${seqs}

mkdir ${tables}
mkdir ${seqs}

# Filter out samples based on individual experiments
qiime feature-table filter-samples \
--i-table ${IN}_ASVtable.qza \
--m-metadata-file ${META} \
--p-where "[Study]='NNK-FMT'" \
--o-filtered-table ${tables}/${EXP}_ASVtable.qza

# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${oldSeqs}_ASVseqs.qza \
--i-table ${tables}/${EXP}_ASVtable.qza \
--o-filtered-data ${seqs}/${EXP}_ASVseqs.qza

#######################################################

study=Dudeja_004_V4_analysis
EXP=Smoke-FMT

IN=${inputDirPath}/${study}/tables/${study}
META=${outputDirPath}/${EXP}/${EXP}_metadata.txt
tables=${outputDirPath}/${EXP}/tables
seqs=${outputDirPath}/${EXP}/seqs
oldSeqs=${inputDirPath}/${study}/seqs/${study}

rm -r ${tables}
rm -r ${seqs}

mkdir ${tables}
mkdir ${seqs}

# Filter out samples based on individual experiments
qiime feature-table filter-samples \
--i-table ${IN}_ASVtable.qza \
--m-metadata-file ${META} \
--p-where "[Study]='Smoke-FMT'" \
--o-filtered-table ${tables}/${EXP}_ASVtable.qza

# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${oldSeqs}_ASVseqs.qza \
--i-table ${tables}/${EXP}_ASVtable.qza \
--o-filtered-data ${seqs}/${EXP}_ASVseqs.qza

#######################################################

study=Dudeja_004_V4_analysis
EXP=Abx-Colon

IN=${inputDirPath}/${study}/tables/${study}
META=${outputDirPath}/${EXP}/${EXP}_metadata.txt
tables=${outputDirPath}/${EXP}/tables
seqs=${outputDirPath}/${EXP}/seqs
oldSeqs=${inputDirPath}/${study}/seqs/${study}

rm -r ${tables}
rm -r ${seqs}

mkdir ${tables}
mkdir ${seqs}

# Filter out samples based on individual experiments
qiime feature-table filter-samples \
--i-table ${IN}_ASVtable.qza \
--m-metadata-file ${META} \
--p-where "[Study]='Abx-Colon'" \
--o-filtered-table ${tables}/${EXP}_ASVtable.qza

# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${oldSeqs}_ASVseqs.qza \
--i-table ${tables}/${EXP}_ASVtable.qza \
--o-filtered-data ${seqs}/${EXP}_ASVseqs.qza

#######################################################

study=Dudeja_002_V4_analysis_qiime
EXP=Abx-Pancreatic

IN=${inputDirPath}/${study}/tables/${study}
META=${outputDirPath}/${EXP}/${EXP}_metadata.txt
tables=${outputDirPath}/${EXP}/tables
seqs=${outputDirPath}/${EXP}/seqs
oldSeqs=${inputDirPath}/${study}/seqs/${study}

rm -r ${tables}
rm -r ${seqs}

mkdir ${tables}
mkdir ${seqs}

# Filter out samples based on individual experiments
qiime feature-table filter-samples \
--i-table ${IN}_ASVtable.qza \
--m-metadata-file ${META} \
--p-where "[Study]='Abx-Pancreatic'" \
--o-filtered-table ${tables}/${EXP}_ASVtable.qza

# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${oldSeqs}_ASVseqs.qza \
--i-table ${tables}/${EXP}_ASVtable.qza \
--o-filtered-data ${seqs}/${EXP}_ASVseqs.qza

#######################################################

study=Dudeja_002_V4_analysis_qiime
EXP=Abx-Bladder

IN=${inputDirPath}/${study}/tables/${study}
META=${outputDirPath}/${EXP}/${EXP}_metadata.txt
tables=${outputDirPath}/${EXP}/tables
seqs=${outputDirPath}/${EXP}/seqs
oldSeqs=${inputDirPath}/${study}/seqs/${study}

rm -r ${tables}
rm -r ${seqs}

mkdir ${tables}
mkdir ${seqs}

# Filter out samples based on individual experiments
qiime feature-table filter-samples \
--i-table ${IN}_ASVtable.qza \
--m-metadata-file ${META} \
--p-where "[Study]='Abx-Bladder'" \
--o-filtered-table ${tables}/${EXP}_ASVtable.qza

# Remove sequences filtered out of ASV table
qiime feature-table filter-seqs \
--i-data ${oldSeqs}_ASVseqs.qza \
--i-table ${tables}/${EXP}_ASVtable.qza \
--o-filtered-data ${seqs}/${EXP}_ASVseqs.qza

#######################################################
