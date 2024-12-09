#!/bin/bash

# BioLockJ v1.4.2: ${scriptDir}/MAIN_06_PCoA_Neo-Pancreatic_CSEvCtrl.sh

export BLJ=/app/biolockj
export HOSTNAME=c6698da27901

pipeDir="/mnt/efs/pipelines/mouseSmoke_analysis_2024Dec09"
modDir="${pipeDir}/06_PCoA_Neo-Pancreatic_CSEvCtrl"
scriptDir="${modDir}/script"
tempDir="${modDir}/temp"
logDir="${modDir}/log"
outputDir="${modDir}/output"

touch "${scriptDir}/MAIN_06_PCoA_Neo-Pancreatic_CSEvCtrl.sh_Started"

exec 1>${logDir}/MAIN.log
exec 2>&1
cd ${scriptDir}

# Spawn Docker container
function spawnDockerContainer() {
    SCRIPT_ID=$(basename $1)
    containerId=$(docker run --detach --rm \
     -v /Users/aliciasorgen/BioLockJ_pipelines:/mnt/efs/pipelines:delegated \
     -v /Users/aliciasorgen/git/Dudeja_Projects/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023/analysis/BLJ_config_files:/mnt/efs/vol_1:ro \
     -v /Users/aliciasorgen/git/Dudeja_Projects/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023/analysis/RScripts:/mnt/efs/vol_2:ro \
     -v /Users/aliciasorgen/git/Dudeja_Projects/Smoking_Cancer_Gut_Dysbiosis_Analysis_2023/analysis/data:/mnt/efs/vol_3:ro \
     -e BLJ=/app/biolockj  \
     -e HOSTNAME=c6698da27901  \
     asorgen/gg-tidyr:v5 \
    /bin/bash -c "$1" )
    echo "Launched docker image: asorgen/gg-tidyr:v5"
    echo "To execute module: GenMod"
    echo "Docker container id: $containerId"
    echo "${SCRIPT_ID}:docker:${containerId}" >> ${scriptDir}/MAIN_06_PCoA_Neo-Pancreatic_CSEvCtrl.sh_Started
    docker inspect ${containerId}
}

function scriptFailed() {
    echo "Line #${2} failure status code [ ${3} ]:  ${1}" >> "${scriptDir}/MAIN_06_PCoA_Neo-Pancreatic_CSEvCtrl.sh_Failures"
    exit ${3}
}

function executeLine() {
    ${1}
    statusCode=$?
    [ ${statusCode} -ne 0 ] && scriptFailed "${1}" ${2} ${statusCode}
}

executeLine "spawnDockerContainer ${scriptDir}/06.0_GenMod.sh" ${LINENO}

touch "${scriptDir}/MAIN_06_PCoA_Neo-Pancreatic_CSEvCtrl.sh_Success"
