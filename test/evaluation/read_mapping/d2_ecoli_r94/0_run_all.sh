pushd ../../../..
source ensure_full_environment.sh
popd

#if $1 is passed, use that as the number of threads, otherwise default to 64
THREADS=64	
if [ ! -z "$1" ]; then
    THREADS=$1
fi

bash run_minimap2.sh $THREADS
bash run_sigmap.sh $THREADS
bash run_uncalled.sh $THREADS
bash run_rawhash.sh $THREADS
bash run_rawalign_full_global.sh $THREADS
bash run_rawalign_full_sparse.sh $THREADS
bash run_rawalign_banded_sparse.sh $THREADS
bash run_rawalign_banded_sparse_nominanchor.sh $THREADS
bash run_rawalign_banded_sparse_nominanchor_anchorchainalignment.sh $THREADS

#list of band-radius-fractions to test
BAND_RADIUS_FRACTIONS=""
BAND_RADIUS_FRACTIONS=$BAND_RADIUS_FRACTIONS" 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09"
BAND_RADIUS_FRACTIONS=$BAND_RADIUS_FRACTIONS" 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19"
BAND_RADIUS_FRACTIONS=$BAND_RADIUS_FRACTIONS" 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00"

#run rawalign with different band radius fractions
for BAND_RADIUS_FRACTION in $BAND_RADIUS_FRACTIONS; do
    bash run_rawalign_banded_sparse_nominanchor.sh $THREADS $BAND_RADIUS_FRACTION
done

BAND_RADIUS_FRACTION=0.10
MATCH_BONUSES="0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.50"
DTW_MIN_SCORES="0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0"

#run rawalign with different band radius fractions
#for BAND_RADIUS_FRACTION in $BAND_RADIUS_FRACTIONS; do
for MATCH_BONUS in $MATCH_BONUSES; do
    for DTW_MIN_SCORE in $DTW_MIN_SCORES; do
        bash run_rawalign_banded_sparse_nominanchor.sh $THREADS $BAND_RADIUS_FRACTION $MATCH_BONUS $DTW_MIN_SCORE
    done
done
