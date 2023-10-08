pushd ../../..
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
bash run_rawalign_banded_sparse_nominanchor.sh $THREADS
