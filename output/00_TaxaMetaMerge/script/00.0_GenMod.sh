#!/bin/bash

# BioLockJ.v1.4.2: ${scriptDir}/00.0_GenMod.sh

export BLJ=/app/biolockj
export HOSTNAME=4d28528f15d3

pipeDir="/mnt/efs/pipelines/mouseSmoke_analysis_2024Dec09"
modDir="${pipeDir}/00_TaxaMetaMerge"
scriptDir="${modDir}/script"
tempDir="${modDir}/temp"
logDir="${modDir}/log"
outputDir="${modDir}/output"

touch "${scriptDir}/00.0_GenMod.sh_Started"

exec 1>${logDir}/00.0_GenMod.log
exec 2>&1

cd ${scriptDir}

function scriptFailed() {
    echo "Line #${2} failure status code [ ${3} ]:  ${1}" >> "${scriptDir}/00.0_GenMod.sh_Failures"
    exit ${3}
}

function executeLine() {
    ${1}
    statusCode=$?
    [ ${statusCode} -ne 0 ] && scriptFailed "${1}" ${2} ${statusCode}
}

executeLine "Rscript ${modDir}/resources/TaxaMetaMerge.R Abx-Bladder Abx-Colon Abx-Pancreatic Neo-Pancreatic NNK-FMT Smoke-FMT" ${LINENO}
touch "${scriptDir}/00.0_GenMod.sh_Success"
echo 'Created Success flag.'
