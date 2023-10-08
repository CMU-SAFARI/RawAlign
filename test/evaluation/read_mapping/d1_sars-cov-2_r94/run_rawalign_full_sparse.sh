#!/bin/bash

THREAD=$1

#d1_sars-cov-2_r94
OUTDIR="./outdir/"
FAST5="../../../data/d1_sars-cov-2_r94/fast5_files/"
REF="../../../data/d1_sars-cov-2_r94/ref.fa"
PORE="../../../../extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="d1_sars-cov-2_r94"
mkdir -p ${OUTDIR}

#Viral preset (default for viral genomes)
PRESET="viral"
PARAMS="--dtw-evaluate-chains --dtw-border-constraint sparse"

#convert params to a sanitized string
#spaces are replaced by underscores
#hyphens are removed
SANITIZED_PARAMS=$(echo ${PARAMS} | sed 's/ /_/g' | sed 's/-//g')

#add leading underscore if params are not empty
if [ ! -z "${SANITIZED_PARAMS}" ]; then
    SANITIZED_PARAMS="_${SANITIZED_PARAMS}"
fi
bash ../../../scripts/run_rawalign.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawalign_${PRESET}${SANITIZED_PARAMS}.out" 2> "${OUTDIR}/${PREFIX}_rawalign_${PRESET}${SANITIZED_PARAMS}.err"
