# Evaluation - Contamination
We assume your current directory points to this directory ([`contamination`](./))

## Organization

This directory includes the scripts to run RawAlign, UNCALLED, Sigmap, RawHash, and minimap2 for the corresponding dataset where the dataset name is specified in the directory name. To run all tools and configurations, you can simply run the `0_run_all.sh`.

After running minimap2, each directory should include `true_mappings.paf` that can be used as the ground truth mapping information for the corresponding dataset.

After running RawAlign, UNCALLED, Sigmap, and RawHash, each directory should include a subdirectory called `outdir` which includes the mapping output of each tool in `PAF` format.

### Comparing RawAlign to UNCALLED, Sigmap, and RawHash

This directory includes a subdirectory called `comparison` where there are three scripts: `0_run.sh`, `1_generate_results.sh`, and `2_output_results.sh`. These scripts use the output generated by each tool for a certain dataset to compare them each other.

## Running and comparing the tools
Use the following commands to run all tools and compare them to each other:

```bash
bash 0_run_all.sh && cd comparison &&  bash 0_run.sh && cd ..
```