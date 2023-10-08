import argparse
import subprocess
import re
import numpy as np
import math
from pathlib import Path


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


def dataset_name_to_pretty_name(dataset_name):
    if 'covid' in dataset_name:
        return 'd1 SARS-CoV-2'
    elif 'ecoli' in dataset_name:
        return 'd2 E.coli'
    elif 'yeast' in dataset_name:
        return 'd3 Yeast'
    elif 'green_algae' in dataset_name:
        return 'd4 Green Algae'
    elif 'human' in dataset_name:
        return 'd5 Human'
    elif 'contamination' in dataset_name:
        return 'Contamination'
    elif 'relative_abundance' in dataset_name:
        return 'Relative Abundance'
    else:
        return dataset_name

def latex_escape(string):
    return string.replace('_', '\_')


def get_latex_table_single_tabular_contents(dataset_results, dataset_name = '', omit_metric_names=False, omit_first_line=False):
    tputs = dataset_results['tputs']
    latency = dataset_results['latency']
    accs = dataset_results['accs']
    sequencing_latencies = dataset_results['bp_latencies']
    memory_footprints = dataset_results['memory_footprints']

    toolnames = list(tputs.keys())
    for toolname in toolnames:
        assert(toolname in tputs)
        assert(toolname in latency)
        assert(toolname in accs)
        assert(toolname in sequencing_latencies)
        assert(toolname in memory_footprints)
    metric_names = set(['Indexing Memory\nFootprint (GB)', 'Mapping Memory\nFootprint (GB)', 'Mean Throughput\n(bp/s)', 'Median Throughput\n(bp/s)', 'Analysis\nLatency (ms)', 'Mean Sequencing\nLatency (Bases)', 'Mean Sequencing\nLatency (Chunks)', 'Precision', 'Recall', 'F-1'])
    ordered_metric_names = ['Indexing Memory\nFootprint (GB)', 'Mapping Memory\nFootprint (GB)', 'Mean Throughput\n(bp/s)', 'Median Throughput\n(bp/s)', 'Analysis\nLatency (ms)', 'Mean Sequencing\nLatency (Bases)', 'Mean Sequencing\nLatency (Chunks)', 'Precision', 'Recall', 'F-1']
    ordered_tools = ['Uncalled', 'Sigmap', 'RawHash', 'RawAlign']

    metric_optimization_functions = {
        'Indexing Memory\nFootprint (GB)': min,
        'Mapping Memory\nFootprint (GB)': min,
        'Mean Throughput\n(bp/s)': max,
        'Median Throughput\n(bp/s)': max,
        'Analysis\nLatency (ms)': min,
        'Mean Sequencing\nLatency (Bases)': min,
        'Mean Sequencing\nLatency (Chunks)': min,
        'Precision': max,
        'Recall': max,
        'F-1': max
    }

    tool_datas = {}
    for toolname in toolnames:
        tool_data = {
            'Mean Throughput\n(bp/s)': tputs[toolname][0], #mean throughput, bp/s
            'Median Throughput\n(bp/s)': tputs[toolname][1], #median throughput, bp/s
            'Analysis\nLatency (ms)': latency[toolname], #mean latency, ms, negated so that lower is better
            'Precision': accs[toolname][0], #precision (frac)
            'Recall': accs[toolname][1], #recall (frac)
            'F-1': accs[toolname][2], #f1 score (frac)
            'Mean Sequencing\nLatency (Bases)': sequencing_latencies[toolname][0], #mean base latency, negated so that lower is better
            'Mean Sequencing\nLatency (Chunks)': sequencing_latencies[toolname][1], #mean chunk latency, negated so that lower is better
            'Indexing Memory\nFootprint (GB)': memory_footprints[toolname]['index'], #map memory footprint in GB, negated so that lower is better
            'Mapping Memory\nFootprint (GB)': memory_footprints[toolname]['map'] #map memory footprint in GB, negated so that lower is better
        }
        if 'uncalled' in toolname.lower():
            tool_data['Mean Sequencing\nLatency (Chunks)'] = sequencing_latencies[toolname][0]/450
        tool_datas[tool_name_to_pretty_name(toolname)] = tool_data

    annotated_tool_datas = {toolname:{} for toolname in tool_datas.keys()}
    for metric_name in metric_names:
        metric_values = [tool_datas[toolname][metric_name] for toolname in tool_datas.keys()]
        best_metric_value = metric_optimization_functions[metric_name](metric_values)
        for toolname in tool_datas.keys():
            annotated_tool_datas[toolname][metric_name] = (tool_datas[toolname][metric_name], tool_datas[toolname][metric_name] == best_metric_value)

    res = []

    n = '\n'
    if not omit_first_line:
        res.append(r'\hline' + n)
    
    headers = [dataset_name] + (ordered_metric_names if not omit_metric_names else [])#['']*len(metric_names))

    num_header_rows = max([len(header.split('\n')) for header in headers])
    padded_headers = []
    for header in headers:
        current_header_lines = header.split('\n')
        num_padding_lines = num_header_rows - len(current_header_lines)
        padded_header = '\n'*num_padding_lines + header
        padded_headers.append(padded_header)

    for row_i in range(num_header_rows):
        header_row_cells = []
        for column_i, padded_header in enumerate(padded_headers):
            header_component = padded_header.split('\n')[row_i]
            alignment = 'c' if column_i > 0 else 'l'
            if len(headers) > 1:
                header_row_cells.append(r'\multicolumn{1}{' + alignment + r'}{\textbf{' + header_component + r'}}')
            else:
                header_row_cells.append(r'\multicolumn{2}{' + alignment + r'}{\textbf{' + header_component + r'}}')
        res.append(' & '.join(header_row_cells))
        res.append(r'\\' + n)
        
    res.append(r'\hline' + n)
    for tool in ordered_tools:
        annotated_metrics = annotated_tool_datas.get(tool)
        if annotated_metrics is None:
            continue
        res.append(tool + ' & ')
        formatted_numbers = [
            f'{annotated_metrics[metric_name][0]:,.3f}'
            if not math.isinf(annotated_metrics[metric_name][0]) else
            '-'
            for metric_name in ordered_metric_names
        ]
        res.append(' & '.join([
            r'\textbf{' + formatted_number + r'}'
            if annotated_metrics[metric_name][1] else
            formatted_number
            for (metric_name, formatted_number) in zip(ordered_metric_names, formatted_numbers)
        ]))
        res.append(r'\\' + n)
    res.append(r'\hline' + n)

    return ''.join(res)


def write_latex_table(full_results, output_path):
    metric_count = 10
    n = '\n'

    with open(output_path, 'w') as f:
        f.write(r'\begin{tabular}{l' + ' r'*metric_count + '}' + n)

        is_first = True
        for dataset_name, tool_results in full_results.items():
            #if not is_first:
            #    f.write(r'\\' + '\n')
            f.write(get_latex_table_single_tabular_contents(tool_results, dataset_name=dataset_name_to_pretty_name(dataset_name), omit_metric_names=not is_first, omit_first_line=not is_first))
            is_first = False

        f.write(r'\end{tabular}' + n)
    

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

    total_results = {}

    #d1
    results = get_output_results(Path('test/evaluation/read_mapping/d1_sars-cov-2_r94/comparison/'))
    total_results['covid'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    #d2
    results = get_output_results(Path('test/evaluation/read_mapping/d2_ecoli_r94/comparison/'))
    total_results['ecoli'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    #d3
    results = get_output_results(Path('test/evaluation/read_mapping/d3_yeast_r94/comparison/'))
    total_results['yeast'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    #d4
    results = get_output_results(Path('test/evaluation/read_mapping/d4_green_algae_r94/comparison/'))
    total_results['green_algae'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    #d5
    results = get_output_results(Path('test/evaluation/read_mapping/d5_human_na12878_r94/comparison/'))
    total_results['human'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    #contamination
    results = get_output_results(Path('test/evaluation/contamination/comparison/'))
    total_results['contamination'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    #relative abundance
    results = get_output_results(Path('test/evaluation/relative_abundance/comparison/'))
    total_results['relative_abundance'] = dict(
        tputs = filter_metrics(parse_throughputs(results), main_rawalign_filter),
        latency = filter_metrics(parse_analysis_latency(results), main_rawalign_filter),
        accs = filter_metrics(parse_accuracy(results), main_rawalign_filter),
        bp_latencies = filter_metrics(parse_sequencing_latencies(results), main_rawalign_filter),
        memory_footprints = filter_metrics(parse_memory_footprints(results), main_rawalign_filter),
    )

    write_latex_table(total_results, Path('paperfigures/full_results_table.tex'))

