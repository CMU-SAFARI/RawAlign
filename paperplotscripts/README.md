# Reproducing the Figures

After all experimental data is generated as described in the [Test README](../test/README.md), the Python scripts in this directory can be used to reproduce the figures in the paper.

The plotting scripts have the following pip3 dependencies:
* regex
* matplotlib
* numpy

You can install them as follows:
```bash
pip3 install regex matplotlib numpy
```

Then you can run the plotting scripts as follows:
```bash
python3 plot_accuracy_throughput_tradeoff.py
python3 plot_band_radius_parameter_sweep.py
python3 plot_matchbonus_parameter_sweep.py
python3 plot_seeding_chaining_alignment.py
python3 plot_spider_tradeoffs.py
python3 table_full_results.py
python3 table_numeric_results.py
python3 table_relative_abundance.py
```
The resulting PDF figures and LaTeX tables will be placed in [paperfigures](../paperfigures/).
