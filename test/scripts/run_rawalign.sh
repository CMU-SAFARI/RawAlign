#!/bin/bash

#Runs the indexing and mapping step of rawalign, generates the output file as well as the corresponding .time files

OUTDIR=$1 #Path to the output directory to store all the files to be generated
PREFIX=$2 #A string prefix that you want to attach to the file names to use as an identifier for your current run (e.g., mytestrun)
SIGNALS=$3 #Path to the directory that contains the fast5 files
REF=$4 #Path to the reference genome
PORE=$5 #Path to the k-mer model file
PRESETX=$6 #Default preset of rawalign for the run (e.g., viral)
THREAD=$7 #Number of threads to use
PARAMS=$8 #(optional -- you can keep it empty) custom parameters to set on top of the default parameters

#convert params to a sanitized string
#spaces are replaced by underscores
#hyphens are removed
SANITIZED_PARAMS=$(echo ${PARAMS} | sed 's/ /_/g' | sed 's/-//g')

#add leading underscore if params are not empty
if [ ! -z "${SANITIZED_PARAMS}" ]; then
    SANITIZED_PARAMS="_${SANITIZED_PARAMS}"
fi

/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawalign_index_${PRESETX}${SANITIZED_PARAMS}.time" rawalign -x ${PRESETX} -t ${THREAD} ${PARAMS} -p "${PORE}" -d "${OUTDIR}/${PREFIX}_rawalign_${PRESETX}${SANITIZED_PARAMS}.ind" ${REF}
/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawalign_map_${PRESETX}${SANITIZED_PARAMS}.time" rawalign -x ${PRESETX} -t ${THREAD} ${PARAMS} -o "${OUTDIR}/${PREFIX}_rawalign_${PRESETX}${SANITIZED_PARAMS}.paf" "${OUTDIR}/${PREFIX}_rawalign_${PRESETX}${SANITIZED_PARAMS}.ind" ${SIGNALS}
