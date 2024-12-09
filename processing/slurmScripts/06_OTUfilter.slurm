#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=OTUfilt                                    # Job name
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
cutoffs=(77 59 170 77 584 231)

inputDirPath=/users/asorgen/mouseSmoke/QIIME2


x=0

for study in ${studies[@]}; do 


	echo -e "\n${study} 0.01% cutoff is ${cutoffs[x]} \n"

	seqs=${inputDirPath}/${study}/seqs/${study}
	visualizations=${inputDirPath}/${study}/visualizations/${study}
	tables=${inputDirPath}/${study}/tables/${study}

	qiime feature-table filter-features \
	--i-table ${tables}_OR-OTUtable-nochim-nosingle.qza \
	--p-min-frequency  ${cutoffs[x]} \
	--o-filtered-table ${tables}_OR-OTUtable-nochim-nosingle-abun.qza

	qiime feature-table summarize \
	--i-table ${tables}_OR-OTUtable-nochim-nosingle-abun.qza \
	--o-visualization ${visualizations}_OR-OTUtable-nochim-nosingle-abun.qzv

	x=$((x+1))

done

