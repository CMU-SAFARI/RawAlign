#!/bin/bash

THREAD=$1

#d5_human_na12878_r94
OUTDIR="./outdir/"
FAST5="../../../data/d5_human_na12878_r94/fast5_files/"
REF="../../../data/d5_human_na12878_r94/ref.fa"
PREFIX="d5_human_na12878_r94"
mkdir -p ${OUTDIR}

#Default
bash ../../../scripts/run_uncalled.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${THREAD} > "${OUTDIR}/${PREFIX}_uncalled.out" 2> "${OUTDIR}/${PREFIX}_uncalled.err"
