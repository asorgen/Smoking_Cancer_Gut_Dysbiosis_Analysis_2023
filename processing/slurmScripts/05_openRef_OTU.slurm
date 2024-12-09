#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=or_OTU                                    # Job name
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

echo Output file: $outFile

mv ${logDir}/${jobName}.* ${logDir}/oldLogs
mv ${logDir}/oldLogs/${outFile}.* ${logDir}

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load openmpi
module load qiime2/2019.1

studies=(Abx-Pancreatic Abx-Bladder Abx-Colon Smoke-FMT NNK-FMT Neo-Pancreatic)

inputDirPath=/users/asorgen/mouseSmoke/QIIME2

for study in ${studies[@]}; do
	
	mkdir ${inputDirPath}/${study}/visualizations

	seqs=${inputDirPath}/${study}/seqs/${study}
	visualizations=${inputDirPath}/${study}/visualizations/${study}
	tables=${inputDirPath}/${study}/tables/${study}
	chimera=${inputDirPath}/${study}/${study}_OR-chimeraID
	
	## Open reference OTU clustering 
	qiime vsearch cluster-features-open-reference \
	--i-table ${tables}_ASVtable.qza \
	--i-sequences ${seqs}_ASVseqs.qza \
	--i-reference-sequences /users/asorgen/REFdatabase/16SrRNA/silva-138-99-seqs.qza \
	--p-perc-identity 0.99 \
	--o-clustered-table ${tables}_OR-OTUtable.qza \
	--o-clustered-sequences ${seqs}_OR-OTUseqs.qza \
	--o-new-reference-sequences ${seqs}_OR-new-ref-OTUseqs.qza \
	--p-threads 16 \
	--verbose	

	qiime feature-table summarize \
	--i-table ${tables}_OR-OTUtable.qza \
	--o-visualization ${visualizations}_OR-OTUtable.qzv

	## Identify chimeras
	qiime vsearch uchime-denovo \
	--i-table ${tables}_OR-OTUtable.qza \
	--i-sequences ${seqs}_OR-OTUseqs.qza \
	--output-dir ${chimera}

	## Filter out chimeras from OTU table
	qiime feature-table filter-features \
	--i-table ${tables}_OR-OTUtable.qza \
	--m-metadata-file ${chimera}/nonchimeras.qza \
	--o-filtered-table ${tables}_OR-OTUtable-nochim.qza

	qiime feature-table summarize \
	--i-table ${tables}_OR-OTUtable-nochim.qza \
	--o-visualization ${visualizations}_OR-OTUtable-nochim.qzv

	## Filter out chimeras from sequences
	qiime feature-table filter-seqs \
	--i-data ${seqs}_OR-OTUseqs.qza \
	--m-metadata-file ${chimera}/nonchimeras.qza \
	--o-filtered-data ${seqs}_OR-OTUseqs-nochim.qza

	## Remove singletons from OTU table
	qiime feature-table filter-features \
	--i-table ${tables}_OR-OTUtable-nochim.qza \
	--p-min-frequency 2 \
	--o-filtered-table ${tables}_OR-OTUtable-nochim-nosingle.qza

	qiime feature-table summarize \
	--i-table ${tables}_OR-OTUtable-nochim-nosingle.qza \
	--o-visualization ${visualizations}_OR-OTUtable-nochim-nosingle.qzv

	## Filter out singletons from sequences
	qiime feature-table filter-seqs \
	--i-data ${seqs}_OR-OTUseqs-nochim.qza \
	--i-table ${tables}_OR-OTUtable-nochim-nosingle.qza \
	--o-filtered-data ${seqs}_OR-OTUseqs-nochim-nosingle.qza


done
