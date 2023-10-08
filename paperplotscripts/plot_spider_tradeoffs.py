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


#from https://matplotlib.org/stable/gallery/specialty_plots/radar_chart.html
def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            return self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


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


def plot(tputs, latency, accs, sequencing_latencies, memory_footprints, plottitle="Read Mapping Spider Plot"):
    toolnames = list(tputs.keys())
    for toolname in toolnames:
        assert(toolname in tputs)
        assert(toolname in latency)
        assert(toolname in accs)
        assert(toolname in sequencing_latencies)
        assert(toolname in memory_footprints)
    metric_names = set(['Throughput', 'Analysis Latency', 'Accuracy', 'Sequencing Latency', 'Memory Footprint'])
    ordered_metric_names = ['Memory Footprint', 'Throughput', 'Analysis Latency', 'Sequencing Latency', 'Accuracy']

    opt_strings = {
        'Throughput': 'High.',
        'Analysis Latency': 'Low.',
        'Accuracy': 'High.',
        'Sequencing Latency': 'Low.',
        'Memory Footprint': 'Low.',
    }

    N = len(metric_names)
    theta = radar_factory(N, frame='polygon')
    fig, ax = plt.subplots(subplot_kw=dict(projection='radar'))
    #fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    ylim = (0, 1.05)
    ax.set_ylim(*ylim)
    ax.set_rgrids([0.2, 0.4, 0.6, 0.8, 1.0], ['', '', '', '', ''])
    ax.set_title(plottitle)
    ax.spines['polar'].set_visible(False)

    #lines, labels = ax.set_varlabels(ordered_metric_names)
    lines, labels = ax.set_varlabels([''] * len(ordered_metric_names))
    #ax.xaxis.grid(True, color='black')
    for metric_name, angle in zip(ordered_metric_names, theta):
        labelangle = angle + theta[1]/2
        upsidedown = labelangle > 1/2*np.pi and labelangle < 3/2*np.pi
        ha='right'
        va='bottom'
        if upsidedown:
            labelangle += np.pi
            ha='left'
            va='top'
        ax.annotate(f'{opt_strings[metric_name]} {metric_name}',
                xy=(angle+np.pi*0.008, 1.0),
                ha=ha, va=va,
                rotation=np.rad2deg(labelangle), rotation_mode='anchor',
                #bold
                weight='bold',
                fontsize=9,
                )
    ax.scatter(theta, [1.0] * len(theta), [20]*len(theta), marker='o', color='black', zorder=10000)

    tool_datas = {}
    for toolname in toolnames:
        tool_data = {
            'Throughput': tputs[toolname][0], #mean throughput, bp/s
            'Analysis Latency': -latency[toolname], #mean latency, ms, negated so that lower is better
            'Accuracy': accs[toolname][2], #f1 score (frac)
            'Sequencing Latency': -sequencing_latencies[toolname][1], #mean chunk latency, negated so that lower is better
            'Memory Footprint': -memory_footprints[toolname]['map'] #map memory footprint in GB, negated so that lower is better
        }
        if 'uncalled' in toolname.lower():
            tool_data['Sequencing Latency'] = -sequencing_latencies[toolname][0]/450
        tool_datas[toolname] = tool_data

    maxs = {metric_name: -float('inf') for metric_name in metric_names}
    mins = {metric_name: float('inf') for metric_name in metric_names}
    for toolname in toolnames:
        for metric_name in metric_names:
            maxs[metric_name] = max(maxs[metric_name], tool_datas[toolname][metric_name])
            mins[metric_name] = min(mins[metric_name], tool_datas[toolname][metric_name])


    for toolname in toolnames:
        color = tool_name_to_color(toolname)
        normalized_tool_data = [
            (tool_datas[toolname][metric_name] - mins[metric_name]) / (maxs[metric_name] - mins[metric_name])
            for metric_name in ordered_metric_names
        ]

        value_start = 0.0
        value_end = 1.0
        usable_value_range = value_end - value_start
        adjusted_tool_data = [
            value_start + (value * usable_value_range)
            for value in normalized_tool_data
        ]

        pretty_tool_name = tool_name_to_pretty_name(toolname)
        zorder = 3 if pretty_tool_name == 'RawAlign' else 1
        alpha = 0.3 if pretty_tool_name == 'RawAlign' else 0.25
        linewidth = 3 if pretty_tool_name == 'RawAlign' else 1.5
        linestyle = '-' if pretty_tool_name == 'RawAlign' else '--'
        ax.plot(theta, adjusted_tool_data, color=color, label=pretty_tool_name, linewidth=linewidth, zorder=zorder+1, linestyle=linestyle)
        ax.fill(theta, adjusted_tool_data, facecolor=color, alpha=alpha, zorder=zorder)

    #fig.legend()

    # Save the figure
    figspath = pathlib.Path("paperfigures")
    figspath.mkdir(exist_ok=True)
    plottitle_slug = plottitle.lower().replace(' ', '_').replace('.', '')
    figpath = figspath / f"{plottitle_slug}.pdf"
    fig.savefig(figpath, dpi=300, bbox_inches='tight')


def plot_legend():
    #generate a fake plot using the same colors as the real plot
    dummy_toolnames = [
        'sigmap',
        'uncalled',
        'rawhash',
        'rawalign',
    ]
    legendfig = plt.figure()

    patches = []
    toolnames = []
    for toolname in dummy_toolnames:
        pretty_tool_name = tool_name_to_pretty_name(toolname)
        alpha = 0.3 if pretty_tool_name == 'RawAlign' else 0.25
        linewidth = 3 if pretty_tool_name == 'RawAlign' else 1.5
        linestyle = '-' if pretty_tool_name == 'RawAlign' else '--'
        color = tool_name_to_color(toolname)
        facecolor = np.array(color)
        facecolor[3] = alpha
        patch = Rectangle((0,0), 1, 1, facecolor=facecolor, edgecolor=color, linewidth=1.5, linestyle=linestyle)
        patches.append(patch)
    legendfig.legend(
        patches,
        [tool_name_to_pretty_name(t) for t in dummy_toolnames], loc='center', ncol=1, fontsize=12,
        #modern look
        frameon=False,
    )

    # Save the figure
    figspath = pathlib.Path("paperfigures")
    figspath.mkdir(exist_ok=True)
    figpath = figspath / f"tools_legend.pdf"
    legendfig.savefig(figpath, dpi=300, bbox_inches='tight', pad_inches=0.1)


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
        'dtwevaluatechains_dtwborderconstraint_sparse_dtwfillmethod_banded=0.10_dtwmatchbonus_0.4_dtwminscore_20.0_stopminanchor_2',
    ]
    main_rawalign_filter = {'rawalign': main_rawalign_configs}

    #d1
    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d1_sars-cov-2_r94/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='d1 SARS-CoV-2',
    )

    #d2
    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d2_ecoli_r94/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='d2 E.coli',
    )

    #d3
    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d3_yeast_r94/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='d3 Yeast',
    )

    #d4
    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d4_green_algae_r94/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='d4 Green Algae',
    )

    #d5
    results = get_output_results(pathlib.Path('test/evaluation/read_mapping/d5_human_na12878_r94/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='d5 Human',
    )

    #contamination
    results = get_output_results(pathlib.Path('test/evaluation/contamination/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='Contamination',
    )

    #relative abundance
    results = get_output_results(pathlib.Path('test/evaluation/relative_abundance/comparison/'))
    tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter)
    latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter)
    accs = filter_metrics(parse_accuracy(results), main_rawalign_filter)
    bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter)
    memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter)
    plot(tputs, latency, accs, bp_latencies, memory_footprints,
         plottitle='Relative Abundance',
    )

    plot_legend()

    plt.show()
