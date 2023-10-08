#!/bin/bash

THREAD=$1

#contamination
OUTDIR="./outdir/"
FAST5="../../data/contamination/fast5_files/"
REF="../../data/contamination/ref.fa"
PORE="../../../extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="contamination"
mkdir -p ${OUTDIR}

#Viral preset (default for viral genomes)
PRESET="viral"
bash ../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
