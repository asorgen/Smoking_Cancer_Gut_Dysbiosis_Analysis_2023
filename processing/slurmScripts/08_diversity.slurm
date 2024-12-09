#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=diversity                                    # Job name
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
	diversity=${inputDirPath}/${study}/${study}_diversity

	## Generate tree for phylogenetic diversity analyses
	qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences ${seqs}_final-filtered-OR-OTUseqs.qza \
	--o-alignment ${seqs}_mafft-aligned-OR-seqs.qza \
	--o-masked-alignment ${seqs}_masked-mafft-aligned-OR-seqs.qza \
	--o-tree ${seqs}_OR-unrooted-tree.qza \
	--o-rooted-tree ${seqs}_OR-rooted-tree.qza



	 ## Calculate alpha and beta diversity based on phylogeny
	 qiime diversity core-metrics-phylogenetic \
	 --i-phylogeny ${seqs}_OR-rooted-tree.qza \
	 --i-table ${tables}_final-filtered-OR-OTUtable.qza \
	 --p-sampling-depth 11500 \
	 --m-metadata-file ${metadata} \
	 --output-dir ${diversity}


	 ## 18.) Calculate Shannon alpha diversity and observed OTUs (maybe not necessary)
	 qiime diversity alpha \
	 --i-table ${tables}_final-filtered-OR-OTUtable.qza \
	 --p-metric 'shannon' \
	 --o-alpha-diversity ${tables}_shannon-diversity.qza

	 qiime metadata tabulate \
	 --m-input-file ${tables}_shannon-diversity.qza \
	 --o-visualization ${visualizations}_shannon-diversity.qzv

	 qiime diversity alpha \
	 --i-table ${tables}_final-filtered-OR-OTUtable.qza \
	 --p-metric 'observed_otus' \
	 --o-alpha-diversity ${tables}_observed-otus.qza

	 qiime metadata tabulate \
	 --m-input-file ${tables}_observed-otus.qza \
	 --o-visualization ${visualizations}_observed-otus.qzv


done

