#!/bin/bash

mkdir -p ../annotated
mkdir -p ../throughput
for i in ../outdir/*.paf; do
    FILENAME=$(basename "${i%.paf}")
    ANNOTATED=../annotated/${FILENAME}_ann.paf
    THROUGHPUT=../throughput/${FILENAME}.throughput
    
    #only annotate if i is newer than ANNOTATED
    if [ -f "$ANNOTATED" ] && [ "$ANNOTATED" -nt "$i" ]; then
        echo "Skipping $i, because it was already annotated"
        continue
    else
        echo "Annotating $i"
    fi
    
    uncalled pafstats -r ../true_mappings.paf --annotate $i > ${ANNOTATED} 2> ${THROUGHPUT}
done

paf_categories=""
#iterate over file and prepend Uncalled= if the filename contains uncalled, and add to paf_categories
for i in ../annotated/*.paf; do
    FILENAME=$(basename "${i%.paf}")
    if [[ $FILENAME == *"uncalled"* ]]; then
        categorized="Uncalled=$i"
    elif [[ $FILENAME == *"sigmap"* ]]; then
        categorized="Sigmap=$i"
    elif [[ $FILENAME == *"rawhash"* ]]; then
        categorized="RawHash=$i"
    elif [[ $FILENAME == *"rawalign"* ]]; then
        categorized="RawAlign=$i"
    else
        echo "Error: could not detect paf catogory of $i"
        exit 1
    fi
    paf_categories="${paf_categories}${categorized} "
done

python ../../../../scripts/compare_pafs.py \
    ${paf_categories} \
    > d1_sars-cov-2_r94.comparison
