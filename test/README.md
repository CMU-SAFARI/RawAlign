# Reproducing the results

## Prerequisites

We compare RawAlign with [UNCALLED](https://github.com/skovaka/UNCALLED), [Sigmap](https://github.com/haowenz/sigmap), and [RawHash](https://github.com/CMU-SAFARI/RawHash). We specify the versions we use for each tool below.

We list the links to download and compile each tool for comparison below:

* [UNCALLED](https://github.com/skovaka/UNCALLED/tree/4c0584ab60a811e74664b0b0d241257d39b967ae)
* [Sigmap](https://github.com/haowenz/sigmap/commit/c9a40483264c9514587a36555b5af48d3f054f6f)
* [RawHash](https://github.com/CMU-SAFARI/RawHash/commit/e9a56fec18007341231af89b32eec9578b1ba622)

We use minimap2 to generate the ground truth mapping information by mapping basecalled seqeunces to their corresponding reference genomes. We use the following minimap2 version:

* [Minimap2 v2.24](https://github.com/lh3/minimap2/releases/tag/v2.24)

We use various tools to process and analyze the data we generate using each tool. The following tools must also be installed in your machine. We suggest using conda to install these tools with their specified versions as almost all of them are included in the conda repository.

* Python 3.6+ (We worked with Python 3.10.6)
* pip3 (We worked with 22.0.2)
* ONT VBZ HDF Plugin - We provide the [`ensure_hdf5_vbz_plugin.sh`](../ensure_hdf5_vbz_plugin.sh) script to install this plugin and add it to the `PATH`.

### Shortcut
We recommend the [`ensure_full_environment.sh`](../ensure_full_environment.sh) script to install all baselines, tools, dependencies, and environment variables.

# Datasets

Please follow the instructions in the [Datasets README](./data/README.md).

# Evaluation

Please follow the instructions in the [Evaluation README](./evaluation/README.md).

# Reproducing the Figures

Please follow the instructions in the [Paper Plot Scripts README](../paperplotscripts/README.md).
