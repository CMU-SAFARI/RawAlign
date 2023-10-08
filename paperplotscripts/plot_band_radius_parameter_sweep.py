import argparse
import subprocess
import re
import numpy as np
import pathlib
from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


def parse_throughputs(output_results):
    throughput_pattern = r"^(\S+)\s+BP per sec:\s+(\d+\.\d+)\s+(\d+\.\d+)$"
    throughput_matches = re.findall(throughput_pattern, output_results, re.MULTILINE)

    res = {}
    for match in throughput_matches:
        tool_name = match[0]
        mean_throughput = float(match[1])
        median_throughput = float(match[2])
        res[tool_name] = (mean_throughput, median_throughput)
    
    return res


def parse_analysis_latency(output_results):
    latency_pattern = r"^(\S+)\s+Mean time per read :\s+(\d+\.\d+)$"
    latency_matches = re.findall(latency_pattern, output_results, re.MULTILINE)

    res = {}
    for match in latency_matches:
        tool_name = match[0]
        mean_latency = float(match[1])
        res[tool_name] = mean_latency
    
    return res


def parse_accuracy(output_results):
    precision_pattern = r"^(\S+)\s+precision:\s+(\d+\.\d+)$"
    recall_pattern = r"^(\S+)\s+recall:\s+(\d+\.\d+)$"
    f1_pattern = r"^(\S+)\s+F-1 score:\s+(\d+\.\d+)$"

    precision_matches = re.findall(precision_pattern, output_results, re.MULTILINE)
    recall_matches = re.findall(recall_pattern, output_results, re.MULTILINE)
    f1_matches = re.findall(f1_pattern, output_results, re.MULTILINE)

    precisions = {}
    recalls = {}
    f1s = {}
    for dict, matches in [(precisions, precision_matches), (recalls, recall_matches), (f1s, f1_matches)]:
        for match in matches:
            tool_name = match[0]
            score = float(match[1])
            dict[tool_name] = score
    
    #set intersection
    tools = set(precisions.keys()) & set(recalls.keys()) & set(f1s.keys())
    res = {}
    for tool in tools:
        res[tool] = (precisions[tool], recalls[tool], f1s[tool])

    return res


def parse_sequencing_latencies(output_results):
    bp_latency_pattern = r"^(\S+)\s+Mean # of sequenced bases per read :\s+((?:\d+\.\d+)|inf)$"
    bp_latency_matches = re.findall(bp_latency_pattern, output_results, re.MULTILINE)

    chunk_latency_pattern = r"^(\S+)\s+Mean # of sequenced chunks per read :\s+((?:\d+\.\d+)|inf)$"
    chunk_latency_matches = re.findall(chunk_latency_pattern, output_results, re.MULTILINE)

    bp_latencies = {}
    for match in bp_latency_matches:
        tool_name = match[0]
        mean_bp_latency = float(match[1])
        bp_latencies[tool_name] = mean_bp_latency
    
    chunk_latencies = {}
    for match in chunk_latency_matches:
        tool_name = match[0]
        mean_chunk_latency = float(match[1])
        chunk_latencies[tool_name] = mean_chunk_latency

    #merge
    tools = set(bp_latencies.keys()) & set(chunk_latencies.keys())
    res = {}
    for tool in tools:
        res[tool] = (bp_latencies[tool], chunk_latencies[tool])

    return res


def parse_memory_footprints(output_results):
    memory_footprint_pattern = r"^(\S+)\s+Memory \(GB\):\s+(\d+\.\d+)$"
    memory_footprint_matches = re.findall(memory_footprint_pattern, output_results, re.MULTILINE)

    res = {}
    for match in memory_footprint_matches:
        tool_name = match[0]
        if '_map' in tool_name:
            tool_name = ''.join(tool_name.split('_map'))
            footprint_type = 'map'
        elif '_index' in tool_name:
            tool_name = ''.join(tool_name.split('_index'))
            footprint_type = 'index'
        else:
            raise ValueError(f"Could not parse memory footprint type (map/index) from {tool_name}")
        memory_footprint = float(match[1])
        #res[tool_name] = (mean_memory_footprint, median_memory_footprint)
        new_entry = res.get(tool_name, {})
        new_entry[footprint_type] = memory_footprint
        res[tool_name] = new_entry

    return res


def parse_rawalign_config(config_string):
    config = {}
    config['preset'] = config_string.split('_')[0]

    if 'dtwfillmethod_' in config_string:
        config['dtwfillmethod'] = config_string.split('dtwfillmethod_')[-1].split('_')[0]
    if 'dtwfillmethod_banded=' in config_string:
        config['band_radius'] = float(config_string.split('dtwfillmethod_banded=')[-1].split('_')[0])
    if 'dtwborderconstraint_' in config_string:
        config['dtwborderconstraint'] = config_string.split('dtwborderconstraint_')[-1].split('_')[0]
    if 'dtwmatchbonus_' in config_string:
        config['dtwmatchbonus'] = float(config_string.split('dtwmatchbonus_')[-1].split('_')[0])
    if 'dtwminscore' in config_string:
        config['dtwminscore'] = float(config_string.split('dtwminscore_')[-1].split('_')[0])
    if 'stopminanchor' in config_string:
        config['stopminanchor'] = int(config_string.split('stopminanchor_')[-1].split('_')[0])

    return config


def tool_name_to_pretty_name(tool_name):
    if 'rawalign' in tool_name:
        return 'RawAlign'
    elif 'uncalled' in tool_name:
        return 'Uncalled'
    elif 'rawhash' in tool_name:
        return 'RawHash'
    elif 'sigmap' in tool_name:
        return 'Sigmap'
    else:
        return tool_name


def tool_name_to_color(tool_name):
    #use matplotlib color map to get a color for each tool
    #cm = colormaps['Dark2']
    #colors = colormaps['Dark2'](np.linspace(0, 1, 8))
    colors = colormaps['tab10'](np.linspace(0, 1, 10))
    keys_to_look_for = ['rawalign', 'uncalled', 'rawhash', 'sigmap']
    for i, key in enumerate(keys_to_look_for):
        if key in tool_name.lower():
            return colors[i]
    return (0.0, 0.0, 0.0)


def metric_name_to_color(metric_name):
    colors = colormaps['tab10'](np.linspace(0, 1, 10))
    keys_to_look_for = ['throughput', 'f1', 'recall', 'precision']

    for i, key in enumerate(keys_to_look_for):
        if key in metric_name.lower():
            return colors[i]

    return (0.0, 0.0, 0.0)


def plot(tputs, latency, accs, sequencing_latencies, memory_footprints, plottitle="Read Mapping Spider Plot"):
    toolnames = list(tputs.keys())
    for toolname in toolnames:
        assert(toolname in tputs)
        assert(toolname in latency)
        assert(toolname in accs)
        assert(toolname in sequencing_latencies)
        assert(toolname in memory_footprints)

    xs = []
    f1_ys = []
    precision_ys = []
    recall_ys = []
    tput_ys = []

    for toolname in toolnames:
        assert('rawalign_' in toolname.lower())
        rawalign_config = parse_rawalign_config(toolname.replace('rawalign_', ''))
        if 'dtwmatchbonus' not in rawalign_config or rawalign_config['dtwmatchbonus'] != 0.4: continue
        if 'dtwminscore' not in rawalign_config or rawalign_config['dtwminscore'] != 20.0: continue
        if 'stopminanchor' not in rawalign_config or rawalign_config['stopminanchor'] != 2: continue
        if 'dtwborderconstraint' not in rawalign_config or rawalign_config['dtwborderconstraint'] != 'sparse': continue
        if 'band_radius' not in rawalign_config: continue

        xs.append(rawalign_config['band_radius']*2) #*2 to convert from radius to band width
        tput_ys.append(tputs[toolname][0])
        f1_ys.append(accs[toolname][2])
        precision_ys.append(accs[toolname][0])
        recall_ys.append(accs[toolname][1])


    fig, ax = plt.subplots(figsize=(6, 3))
    ax2 = ax.twinx()

    #ax.set_title(plottitle)
    ax.set_xlabel("Band Width (%)",     weight='bold', fontsize=12)
    ax.set_ylabel("Accuracy (%)",       weight='bold', fontsize=12)
    ax2.set_ylabel("Throughput (BP/s)", weight='bold', fontsize=12)

    #set label paddings to 1
    ax.xaxis.labelpad = 0
    ax.yaxis.labelpad = 0
    ax2.yaxis.labelpad = 3

    major_linewidth = 3.0
    minor_linewidth = 1.0

    ax.grid(True, linewidth=0.5)

    l0 = ax2.plot(xs, tput_ys,       color=metric_name_to_color('throughput'),   linestyle='-',  linewidth=major_linewidth, label='Throughput')
    l1 = ax.plot(xs, f1_ys,          color=metric_name_to_color('f1'),           linestyle='--', linewidth=major_linewidth, label='F-1 Score')
    l2 = ax.plot(xs, precision_ys,   color=metric_name_to_color('precision'),    linestyle='--', linewidth=minor_linewidth, label='Precision')
    l3 = ax.plot(xs, recall_ys,      color=metric_name_to_color('recall'),       linestyle='--', linewidth=minor_linewidth, label='Recall')

    ax.set_xticklabels( [f"{tick:.2f}" for tick in ax.get_xticks() ], weight='bold'    , fontsize=10)
    ax.set_yticklabels( [f"{tick:.2f}" for tick in ax.get_yticks() ], weight='bold'    , fontsize=10)
    ax2.set_yticklabels([f"{int(tick):,}"  for tick in ax2.get_yticks()], weight='bold', fontsize=10)

    ax.tick_params(axis='x', pad=1)
    ax.tick_params(axis='y', pad=0)
    ax2.tick_params(axis='y', pad=1)

    lines = l0 + l1 + l2 + l3
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, loc='lower right',
              edgecolor=(0, 0, 0, 1.),
              facecolor=(1, 1, 1, 1.),
              framealpha=1.0,
    )

    # Save the figure
    figspath = pathlib.Path("paperfigures")
    figspath.mkdir(exist_ok=True)
    plottitle_slug = plottitle.lower().replace(' ', '_').replace('.', '')
    figpath = figspath / f"{plottitle_slug}.pdf"
    fig.savefig(figpath, dpi=300, bbox_inches='tight')


def get_output_results(comparison_directory:Path):
    bash_script_path = "2_output_results.sh"
    res = subprocess.run(["bash", bash_script_path], check=True, cwd=comparison_directory, capture_output=True)
    return res.stdout.decode("utf-8")


def filter_metrics(metrics, filters={}):
    #filters is a dict of shorttoolname => list of fulltoolname_with_config
    #if any shorttoolname is found in a metric, it is kept only if fulltoolname_with_config contains one of the values in the list
    #if no shorttoolname matches the metric, it is kept

    res = {fulltoolname_with_config:metric
            for (fulltoolname_with_config, metric) in metrics.items()
            if not any([ #check if the toolname is in the filter list. If it is, keep anyway
                shorttoolname in fulltoolname_with_config
                for shorttoolname in filters.keys()
            ])
            or any([ #if the toolname is in the filter list, check if fulltoolname_with_config should be kept
                shorttoolname in fulltoolname_with_config
                and any([f in fulltoolname_with_config for f in filters[shorttoolname]])
                for shorttoolname in filters.keys()
            ])
    }
    return res


if __name__ == "__main__":
    main_rawalign_configs = [
        'dtwevaluatechains_dtwborderconstraint_sparse_dtwfillmethod_banded',
    ]
    main_rawalign_filter = {
        'uncalled': [],
        'sigmap': [],
        'rawhash': [],
        'rawalign': main_rawalign_configs,
    }

#    #d1
#    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d1_sars-cov-2_r94/comparison/'))
#    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
#    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
#    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
#    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
#    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
#    plot(tputs, latency, accs, bp_latencies, memory_footprints,
#         plottitle='d1 SARS-CoV-2',
#    )

    #d2
    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d2_ecoli_r94/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='Band Width Parameter Sweep on d2 E.coli',
    )

#    #d3
#    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d3_yeast_r94/comparison/'))
#    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
#    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
#    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
#    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
#    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
#    plot(tputs, latency, accs, bp_latencies, memory_footprints,
#         plottitle='d3 Yeast',
#    )
#
#    #d4
#    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d4_green_algae_r94/comparison/'))
#    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
#    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
#    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
#    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
#    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
#    plot(tputs, latency, accs, bp_latencies, memory_footprints,
#         plottitle='d4 Green Algae',
#    )
#
#    #d5
#    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d5_human_na12878_r94/comparison/'))
#    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
#    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
#    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
#    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
#    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
#    plot(tputs, latency, accs, bp_latencies, memory_footprints,
#         plottitle='d5 Human',
#    )
#
#    #contamination
#    results = get_output_results(pathlib.Path('test/evaluation/contamination/comparison/'))
#    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
#    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
#    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
#    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
#    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
#    plot(tputs, latency, accs, bp_latencies, memory_footprints,
#         plottitle='Contamination',
#    )
#
#    #relative abundance
#    results = get_output_results(pathlib.Path('test/evaluation/relative_abundance/comparison/'))
#    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
#    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
#    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
#    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
#    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
#    plot(tputs, latency, accs, bp_latencies, memory_footprints,
#         plottitle='Relative Abundance',
#    )

    plt.show()
