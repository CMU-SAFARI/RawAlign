# Read Mapping

The scripts and a [README](./read_mapping/README.md) can be found in the [read mapping directory](./read_mapping/) to perform read mapping using UNCALLED, Sigmap and RawHash.

# Relative Abundance Estimation

The scripts and a [README](./relative_abundance/README.md) can be found in the [relative abundance directory](./relative_abundance/) to perform the relative abundance estimation using UNCALLED, Sigmap and RawHash.

# Contamination Analysis

The scripts and a [README](./contamination/README.md) can be found in the [contamination directory](./contamination/) to perform the contamination analysis using UNCALLED, Sigmap and RawHash.

# Running Everything in One Go
The following commands will run all tools and configurations for all datasets and compare them to each other. Note that this will take a long time (days) to finish - you may want to run only a subset of the commands below, or split them across multiple machines. Some of these datasets (e.g., human, relative abundance) have large reference databases, which may require a large amount of memory to run.

```bash
THREADS=64

cd read_mapping
cd d1_sars-cov-2_r94       && bash 0_run_all.sh $THREADS && cd comparison &&  bash 0_run.sh && cd ../..
cd d2_ecoli_r94            && bash 0_run_all.sh $THREADS && cd comparison &&  bash 0_run.sh && cd ../..
cd d3_yeast_r94            && bash 0_run_all.sh $THREADS && cd comparison &&  bash 0_run.sh && cd ../..
cd d4_green_algae_r94      && bash 0_run_all.sh $THREADS && cd comparison &&  bash 0_run.sh && cd ../..
cd d5_human_na12878_r94    && bash 0_run_all.sh $THREADS && cd comparison &&  bash 0_run.sh && cd ../..
cd ..

cd relative_abundance && bash 0_run_all.sh $THREADS && cd ..
cd contamination && bash 0_run_all.sh $THREADS && cd ..
```