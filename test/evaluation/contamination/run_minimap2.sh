#!/bin/bash

THREAD=$1

minimap2 -x map-ont -t ${THREAD} -o true_mappings.paf ../../data/contamination/ref.fa ../../data/d1_sars-cov-2_r94/reads.fasta ../../data/d5_human_na12878_r94/reads.fasta
