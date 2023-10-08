import argparse
import subprocess
import re
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import colormaps


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


def parse_rawalign_config(config_string):
    config = {}
    config['preset'] = config_string.split('_')[0]

    if 'dtwfillmethod_' in config_string:
        config['dtwfillmethod'] = config_string.split('dtwfillmethod_')[-1].split('_')[0]
    if 'dtwfillmethod_banded=' in config_string:
        config['band_radius'] = float(config_string.split('dtwfillmethod_banded=')[-1].split('_')[0])
    if 'dtwborderconstraint_' in config_string:
        config['dtwborderconstraint'] = config_string.split('dtwborderconstraint_')[-1].split('_')[0]

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


def plot(accs, tputs, plottitle="Accuracy-Performance Tradeoff Space", annotation_offsets = {}, better_offset = (0, 0)):
    fig, ax = plt.subplots()
    ax.set_title(plottitle)
    #ax.set_xlabel("Median Throughput (BP/s)", weight='bold', fontsize=15)
    ax.set_xlabel("Mean Throughput (BP/s)", weight='bold', fontsize=15)
    ax.set_ylabel("F-1 Score (%)", weight='bold', fontsize=15)
    ax.grid(True, linewidth=0.5)
    ax.set_axisbelow(True)

    rawalign_accs = {
        toolname:
        accs[toolname] for toolname in accs if "rawalign" in toolname
    }
    rawalign_tputs = {
        toolname:
        tputs[toolname] for toolname in tputs if "rawalign" in toolname
    }

    non_rawalign_accs = {
        toolname:
        accs[toolname] for toolname in accs if "rawalign" not in toolname
    }

    rawalign_xs = []
    rawalign_ys = []
    for toolname in rawalign_accs:
        precision, recall, f1 = rawalign_accs[toolname]
        mean_tput, median_tput = rawalign_tputs[toolname]
        #rawalign_xs.append(median_tput)
        rawalign_xs.append(mean_tput)
        rawalign_ys.append(100*f1)
    
    other_tools = {}
    for toolname in non_rawalign_accs:
        precision, recall, f1 = non_rawalign_accs[toolname]
        mean_tput, median_tput = tputs[toolname]
        #other_tools[toolname] = (median_tput, 100*f1)
        other_tools[toolname] = (mean_tput, 100*f1)

    # Plot the rawalign points
    ax.scatter(rawalign_xs, rawalign_ys, s=8, label="RawAlign", zorder=10, color=tool_name_to_color('RawAlign'))
    color = ax.collections[-1].get_facecolor()[0]
    ann = ax.annotate(
        f"RawAlign",
        (min(rawalign_xs), max(rawalign_ys)),
        fontsize=12,
        xytext=annotation_offsets.get("RawAlign", (0, 2)),
        textcoords="offset points",
        ha='center',
        va='bottom',
        rotation=0,
        color=color,
        weight='bold',
        bbox=dict(fc='white', ec='none', pad=0.5)
    )
    ann.zorder = 3

    for toolname in other_tools:
        median_tput, f1 = other_tools[toolname]
        ax.scatter(median_tput, f1, s=8, label=toolname, zorder=10, color=tool_name_to_color(toolname))
        #get color from scatter plot
        color = ax.collections[-1].get_facecolor()[0]
        ann = ax.annotate(
            f"{tool_name_to_pretty_name(toolname)}",
            (median_tput, f1),
            fontsize=12,
            xytext=annotation_offsets.get(toolname, (0, 2)),
            textcoords="offset points",
            ha='center',
            va='bottom',
            rotation=0,
            color=color,
            weight='bold',
            bbox=dict(fc='white', ec='none', pad=0.5)
        )
        ann.zorder = 2

    max_y = max(max([f1 for _, f1 in other_tools.values()]), max(rawalign_ys))
    min_y_other_tools = min([f1 for _, f1 in other_tools.values()])
    y_range = max_y - min_y_other_tools
    ylower = min_y_other_tools - 0.1*y_range
    yupper = max_y + 0.1*y_range
    ax.set_ylim(ylower, yupper)

    # Fat arrow towards the top right corner saying "better"
    #better_color = (1.0, 0.0, 0.6)
    better_color = (0.0, 0.0, 0.0)
    ann = ax.annotate(
        "better",
        #(max(rawalign_xs), max(rawalign_ys)),
        #top-right corner based on everything plotted so far
        (ax.get_xlim()[1]+better_offset[0], ax.get_ylim()[1]+better_offset[1]),
        fontsize=15,
        xytext=(-50, -50),
        textcoords="offset points",
        ha='center',
        va='bottom',
        rotation=0,
        color=better_color,
        weight='bold',
        arrowprops=dict(arrowstyle="-|>", mutation_scale=30, color=better_color, linewidth=5),
        bbox=dict(fc='white', ec='none', pad=0.5)
    )
    ann.zorder = 1

    # Add the legend
    #ax.legend()

    # Save the figure
    figspath = Path("paperfigures")
    figspath.mkdir(exist_ok=True)
    plottitle_slug = plottitle.lower().replace(' ', '_')
    figpath = figspath / f"{plottitle_slug}.pdf"
    fig.savefig(figpath, dpi=300, bbox_inches='tight')


def get_output_results(comparison_directory:Path):
    bash_script_path = "2_output_results.sh"
    res = subprocess.run(["bash", bash_script_path], check=True, cwd=comparison_directory, capture_output=True)
    return res.stdout.decode("utf-8")


if __name__ == "__main__":
    #d1
    results = get_output_results(Path('test/evaluation/read_mapping/d1_sars-cov-2_r94/comparison/'))
    tputs = parse_throughputs(results)
    accs = parse_accuracy(results)
    plot(accs, tputs,
         plottitle='d1 SARS-CoV-2 R9.4 Accuracy-Performance Tradeoff Space',
         annotation_offsets={
             'rawhash_viral': (-20, -18),
             'uncalled': (20, 2)
         },
        better_offset=(-200000, 0)
    )

    #d2
    results = get_output_results(Path('test/evaluation/read_mapping/d2_ecoli_r94/comparison/'))
    tputs = parse_throughputs(results)
    accs = parse_accuracy(results)
    plot(accs, tputs,
         plottitle='d2 Ecoli R9.4 Accuracy-Performance Tradeoff Space',
         annotation_offsets={
             'rawhash_sensitive': (-20, -18),
             'uncalled': (20, 2),
             'RawAlign': (-20, 2),
        },
        better_offset=(-15000, 0)
    )
    
    #d3
    results = get_output_results(Path('test/evaluation/read_mapping/d3_yeast_r94/comparison/'))
    tputs = parse_throughputs(results)
    accs = parse_accuracy(results)
    plot(accs, tputs,
         plottitle='d3 Yeast R9.4 Accuracy-Performance Tradeoff Space',
         annotation_offsets={
             'rawhash_sensitive': (-20, 2),
             'uncalled': (20, 2),
             'RawAlign': (-20, 2),
        },
        better_offset=(-5000, 0)
    )

    #d4
    results = get_output_results(Path('test/evaluation/read_mapping/d4_green_algae_r94/comparison/'))
    tputs = parse_throughputs(results)
    accs = parse_accuracy(results)
    plot(accs, tputs,
         plottitle='d4 Green Algae R9.4 Accuracy-Performance Tradeoff Space',
         annotation_offsets={
             'uncalled': (-20, 2),
             'sigmap': (30, -8)
        }
    )    

    #d5
    results = get_output_results(Path('test/evaluation/read_mapping/d5_human_na12878_r94/comparison/'))
    tputs = parse_throughputs(results)
    accs = parse_accuracy(results)
    plot(accs, tputs,
         plottitle='d5 Human R9.4 Accuracy-Performance Tradeoff Space',
        annotation_offsets={
            'rawhash_fast': (20, 2),
            'uncalled': (-20, 2),
            'RawAlign': (20, 2),
        }
    )

    #plt.show()
