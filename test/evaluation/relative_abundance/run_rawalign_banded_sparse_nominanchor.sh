#!/bin/bash

THREAD=$1

#if $2 is passed, use it as the band radius fraction, otherwise default to 15%
if [ ! -z "$2" ]; then
    BAND_RADIUS_FRACTION=$2
else
    BAND_RADIUS_FRACTION=0.10
fi

#if $3 is passed, use it as the match bonus, otherwise default to 0.5
if [ ! -z "$3" ]; then
    MATCH_BONUS=$3
else
    MATCH_BONUS=0.4
fi

#if $4 is passed, use it as the dtw min score, otherwise default to 0.0
if [ ! -z "$4" ]; then
    DTW_MIN_SCORE=$4
else
    DTW_MIN_SCORE=20.0
fi

#relative_abundance
OUTDIR="./outdir/"
FAST5="../../data/relative_abundance/fast5_files/"
REF="../../data/relative_abundance/ref.fa"
PORE="../../../extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="relative_abundance"
mkdir -p ${OUTDIR}

#Viral preset (default for viral genomes)
PRESET="fast"
PARAMS="--dtw-evaluate-chains \
        --dtw-border-constraint sparse \
        --dtw-fill-method banded=${BAND_RADIUS_FRACTION} \
        --dtw-match-bonus ${MATCH_BONUS} \
        --dtw-min-score ${DTW_MIN_SCORE} \
        --stop-min-anchor 2"

#convert params to a sanitized string
#spaces are replaced by underscores
#hyphens are removed
SANITIZED_PARAMS=$(echo ${PARAMS} | sed 's/ /_/g' | sed 's/-//g')

#add leading underscore if params are not empty
if [ ! -z "${SANITIZED_PARAMS}" ]; then
    SANITIZED_PARAMS="_${SANITIZED_PARAMS}"
fi
bash ../../scripts/run_rawalign.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawalign_${PRESET}${SANITIZED_PARAMS}.out" 2> "${OUTDIR}/${PREFIX}_rawalign_${PRESET}${SANITIZED_PARAMS}.err"
