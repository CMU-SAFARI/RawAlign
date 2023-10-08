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
bash run_rawalign_banded_sparse.sh $THREADS
bash run_rawalign_banded_sparse_nominanchor.sh $THREADS
bash run_rawalign_banded_global.sh $THREADS
bash run_rawalign_full_global.sh $THREADS
bash run_rawalign_full_sparse.sh $THREADS

BAND_RADIUS_FRACTION=0.10
MATCH_BONUSES="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.8 0.9 1.0"
DTW_MIN_SCORES="0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0"

#run rawalign with different band radius fractions
#for BAND_RADIUS_FRACTION in $BAND_RADIUS_FRACTIONS; do
for MATCH_BONUS in $MATCH_BONUSES; do
    for DTW_MIN_SCORE in $DTW_MIN_SCORES; do
        bash run_rawalign_banded_sparse.sh $THREADS $BAND_RADIUS_FRACTION $MATCH_BONUS $DTW_MIN_SCORE
    done
done
