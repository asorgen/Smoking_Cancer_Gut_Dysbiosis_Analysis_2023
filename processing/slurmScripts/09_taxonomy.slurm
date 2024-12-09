#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=taxonomy                                    # Job name
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

studies=(Abx-Bladder Abx-Colon Abx-Pancreatic Neo-Pancreatic NNK-FMT Smoke-FMT)

inputDirPath=/users/asorgen/mouseSmoke/QIIME2



for study in ${studies[@]}; do 

	seqs=${inputDirPath}/${study}/seqs/${study}
	visualizations=${inputDirPath}/${study}/visualizations/${study}
	tables=${inputDirPath}/${study}/tables/${study}
	metadata=${inputDirPath}/${study}/${study}_metadata.txt

	## Assign taxonomy to OTUs using Silva 
	qiime feature-classifier classify-consensus-vsearch \
	--i-query ${seqs}_final-filtered-OR-OTUseqs.qza \
	--i-reference-reads /users/asorgen/REFdatabase/16SrRNA/silva-138-99-seqs.qza \
	--i-reference-taxonomy /users/asorgen/REFdatabase/16SrRNA/silva-138-99-tax.qza \
	--o-classification ${tables}_OR-taxonomy.qza

	qiime metadata tabulate \
	--m-input-file ${tables}_OR-taxonomy.qza \
	--o-visualization ${visualizations}_OR-taxonomy.qzv


	## Remove all mitochondria/chloroplast identification from taxonomy
	qiime taxa filter-table \
	--i-table ${tables}_final-filtered-OR-OTUtable.qza \
	--i-taxonomy ${tables}_OR-taxonomy.qza \
	--p-exclude mitochondria,chloroplast \
	--o-filtered-table ${tables}_OR-OTUtable.qza

	qiime feature-table summarize \
	--i-table ${tables}_OR-OTUtable.qza \
	--o-visualization ${visualizations}_OR-OTUtable.qzv


	## Make taxonomy barplots
	qiime taxa barplot \
	--i-table ${tables}_OR-OTUtable.qza \
	--m-metadata-file ${metadata} \
	--i-taxonomy ${tables}_OR-taxonomy.qza \
	--o-visualization ${visualizations}_OR-taxonomy-barplot.qzv

done

