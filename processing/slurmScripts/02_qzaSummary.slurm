#!/bin/bash

#SBATCH --partition=Orion                                    # Partition/queue requested on server
#SBATCH --job-name=qzaSum                                    # Job name
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

studies=(UMGC_Neo_16S NNK_FMT_16S Dudeja_004_V4_analysis Dudeja_002_V4_analysis_qiime)

inputDirPath=/scratch/asorgen/16Sdatasets/mouseSmoke


for study in ${studies[@]}; do
	
	echo \n Starting ${study} import \n

	IN=${inputDirPath}/${study}
	OUT=${IN}/visualizations
	
	mkdir ${OUT}


	qiime demux summarize \
	--i-data ${IN}/${study}_raw-seqs.qza \
	--o-visualization ${OUT}/${study}_raw-seqs.qzv

	echo \n ${study} import is complete \n

done
