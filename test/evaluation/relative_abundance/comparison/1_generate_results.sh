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

python ../../../scripts/compare_pafs.py \
    ${paf_categories} \
    > relative_abundance.comparison

rm -f relative_abundance.abundance

get_numreads_numbases_meanlength(){
    DATASET_PATH="$1"
    SEQKIT_RESULT=$(seqkit stat "$DATASET_PATH" --tabular)
    echo "$SEQKIT_RESULT" | awk 'NR==2 {print $4, $5, $7}'
}

read -r COVID_NUMREADS COVID_NUMBASES COVID_MEANLENGTH <<< $(get_numreads_numbases_meanlength ../../../data/d1_sars-cov-2_r94/reads.fasta)
read -r ECOLI_NUMREADS ECOLI_NUMBASES ECOLI_MEANLENGTH <<< $(get_numreads_numbases_meanlength ../../../data/d2_ecoli_r94/reads.fasta)
read -r YEAST_NUMREADS YEAST_NUMBASES YEAST_MEANLENGTH <<< $(get_numreads_numbases_meanlength ../../../data/d3_yeast_r94/reads.fasta)
read -r ALGAE_NUMREADS ALGAE_NUMBASES ALGAE_MEANLENGTH <<< $(get_numreads_numbases_meanlength ../../../data/d4_green_algae_r94/reads.fasta)
read -r HUMAN_NUMREADS HUMAN_NUMBASES HUMAN_MEANLENGTH <<< $(get_numreads_numbases_meanlength ../../../data/d5_human_na12878_r94/reads.fasta)
TOTALREADS=$(($COVID_NUMREADS+$ECOLI_NUMREADS+$YEAST_NUMREADS+$ALGAE_NUMREADS+$HUMAN_NUMREADS))
TOTALBASES=$(($COVID_NUMBASES+$ECOLI_NUMBASES+$YEAST_NUMBASES+$ALGAE_NUMBASES+$HUMAN_NUMBASES))
COVID_READFRAC=$(echo "scale=5; $COVID_NUMREADS / $TOTALREADS" | bc)
ECOLI_READFRAC=$(echo "scale=5; $ECOLI_NUMREADS / $TOTALREADS" | bc)
YEAST_READFRAC=$(echo "scale=5; $YEAST_NUMREADS / $TOTALREADS" | bc)
ALGAE_READFRAC=$(echo "scale=5; $ALGAE_NUMREADS / $TOTALREADS" | bc)
HUMAN_READFRAC=$(echo "scale=5; $HUMAN_NUMREADS / $TOTALREADS" | bc)
COVID_BASEFRAC=$(echo "scale=5; $COVID_NUMBASES / $TOTALBASES" | bc)
ECOLI_BASEFRAC=$(echo "scale=5; $ECOLI_NUMBASES / $TOTALBASES" | bc)
YEAST_BASEFRAC=$(echo "scale=5; $YEAST_NUMBASES / $TOTALBASES" | bc)
ALGAE_BASEFRAC=$(echo "scale=5; $ALGAE_NUMBASES / $TOTALBASES" | bc)
HUMAN_BASEFRAC=$(echo "scale=5; $HUMAN_NUMBASES / $TOTALBASES" | bc)
echo "fastq Ratio of bases: covid: $COVID_BASEFRAC ecoli: $ECOLI_BASEFRAC yeast: $YEAST_BASEFRAC green_algae: $ALGAE_BASEFRAC human: $HUMAN_BASEFRAC" >> relative_abundance.abundance
echo "fastq Mapping lengths: covid: $COVID_MEANLENGTH ecoli: $ECOLI_MEANLENGTH yeast: $YEAST_MEANLENGTH green_algae: $ALGAE_MEANLENGTH human: $HUMAN_MEANLENGTH" >> relative_abundance.abundance
echo "fastq Ratio of reads: covid: $COVID_READFRAC ecoli: $ECOLI_READFRAC yeast: $YEAST_READFRAC green_algae: $ALGAE_READFRAC human: $HUMAN_READFRAC" >> relative_abundance.abundance

#extract 
for i in '../true_mappings.paf' '../outdir/'*.paf ; do
    if test -f $i; then
        FILENAME=$(basename "${i%.paf}");
        TOOLNAME="${FILENAME#relative_abundance_}";
        #echo -n "${TOOLNAME} " >> relative_abundance.abundance;
        awk -v toolname="${TOOLNAME}" \
        'BEGIN{tot=0; covid=0; human=0; ecoli=0; yeast=0; algae=0; totsum=0; covidsum=0; humansum=0; ecolisum=0; yeastsum=0; algaesum=0;}
        {
            if(substr($6,1,5) == "ecoli" && $3 != "*" && $4 != "*" && $4 >= $3){ecoli++; count++; ecolisum+=$4-$3}
            else if(substr($6,1,5) == "yeast" && $3 != "*" && $4 != "*" && $4 >= $3){yeast++; count++; yeastsum+=$4-$3}
            else if(substr($6,1,11) == "green_algae" && $3 != "*" && $4 != "*" && $4 >= $3){algae++; count++; algaesum+=$4-$3}
            else if(substr($6,1,5) == "human" && $3 != "*" && $4 != "*" && $4 >= $3){human++; count++; humansum+=$4-$3}
            else if(substr($6,1,5) == "covid" && $3 != "*" && $4 != "*" && $4 >= $3){covid++; count++; covidsum+=$4-$3}
        }
        END{
            totsum = ecolisum + yeastsum + algaesum + humansum + covidsum;
            if(totsum == 0){
                print toolname " Ratio of bases: covid: 0 ecoli: 0 yeast: 0 green_algae: 0 human: 0";
                print toolname " Mapping lengths: covid: 0 ecoli: 0 yeast: 0 green_algae: 0 human: 0";
                print toolname " Ratio of reads: covid: 0 ecoli: 0 yeast: 0 green_algae: 0 human: 0";
            }
            else{
                print toolname " Ratio of bases: covid: " covidsum/totsum " ecoli: " ecolisum/totsum " yeast: " yeastsum/totsum " green_algae: " algaesum/totsum " human: " humansum/totsum;
                print toolname " Mapping lengths: covid: " covidsum/covid " ecoli: " ecolisum/ecoli " yeast: " yeastsum/yeast " green_algae: " algaesum/algae " human: " humansum/human;
                print toolname " Ratio of reads: covid: " covid/count " ecoli: " ecoli/count " yeast: " yeast/count " green_algae: " algae/count " human: " human/count;
            }
        }' $i >> relative_abundance.abundance
    fi
done;
