#!/bin/bash

# BioLockJ.v1.4.2: ${scriptDir}/15.0_GenMod.sh

export BLJ=/app/biolockj
export HOSTNAME=c6698da27901

pipeDir="/mnt/efs/pipelines/mouseSmoke_analysis_2024Dec09"
modDir="${pipeDir}/15_LinearModeling_Smoke-FMT_CasevCtrl"
scriptDir="${modDir}/script"
tempDir="${modDir}/temp"
logDir="${modDir}/log"
outputDir="${modDir}/output"

touch "${scriptDir}/15.0_GenMod.sh_Started"

exec 1>${logDir}/15.0_GenMod.log
exec 2>&1

cd ${scriptDir}

function scriptFailed() {
    echo "Line #${2} failure status code [ ${3} ]:  ${1}" >> "${scriptDir}/15.0_GenMod.sh_Failures"
    exit ${3}
}

function executeLine() {
    ${1}
    statusCode=$?
    [ ${statusCode} -ne 0 ] && scriptFailed "${1}" ${2} ${statusCode}
}

executeLine "Rscript ${modDir}/resources/LinearModeling_CasevCtrl.R BLJ Smoke-FMT" ${LINENO}
touch "${scriptDir}/15.0_GenMod.sh_Success"
echo 'Created Success flag.'
