#!/bin/bash

NUMBER_OF_CORES=40

# add additional paths to run snakemake on here
PATHS_TO_RUN=("software" "annotation" "cuffdiff2_analysis" "bottomly" "geuvadis" ".")

ARGUMENTS="-p -j ${NUMBER_OF_CORES} --dryrun"

########################################################################

BASE_PATH=$(pwd)

for path in "${PATHS_TO_RUN[@]}"
do
  cd ${path} || exit
  echo "In path: ${path}"
  snakemake ${ARGUMENTS}
  cd ${BASE_PATH} || exit
done
