import subprocess
import re
import pathlib


def get_output_results(comparison_directory:pathlib.Path):
    bash_script_path = "2_output_results.sh"
    res = subprocess.run(["bash", bash_script_path], check=True, cwd=comparison_directory, capture_output=True)
    return res.stdout.decode("utf-8")


def parse_relative_abundances(output_results):
    read_ratio_pattern = r"^(\S+)\s+Ratio of reads:(.*?)$"
    bases_ratio_pattern = r"^(\S+)\s+Ratio of bases:(.*?)$"
    label_value_pattern = r"(\S+):\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)+(?:\s|$)"

    read_ratios = {}
    for line in output_results.split("\n"):
        read_ratio_match = re.search(read_ratio_pattern, line)
        if not read_ratio_match:
            continue
        tool_name = read_ratio_match.group(1)
        label_value_pairs = read_ratio_match.group(2)
        label_value_matches = re.findall(label_value_pattern, label_value_pairs)
        tool_read_ratios = {label: float(value) for label, value in label_value_matches}
        read_ratios[tool_name] = tool_read_ratios

    bases_ratios = {}
    for line in output_results.split("\n"):
        bases_ratio_match = re.search(bases_ratio_pattern, line)
        if not bases_ratio_match:
            continue
        tool_name = bases_ratio_match.group(1)
        label_value_pairs = bases_ratio_match.group(2)
        label_value_matches = re.findall(label_value_pattern, label_value_pairs)
        tool_bases_ratios = {label: float(value) for label, value in label_value_matches}
        bases_ratios[tool_name] = tool_bases_ratios

    #merge
    tools = set(read_ratios.keys()).intersection(set(bases_ratios.keys()))
    res = {}
    for tool in tools:
        res[tool] = {}
        for label in read_ratios[tool].keys():
            res[tool][label] = {"read_ratio": read_ratios[tool][label], "bases_ratio": bases_ratios[tool][label]}

    return res


def get_euclidean_distance(first_tool_relative_abundances, second_tool_relative_abundances):
    read_ratio_diff_squared_sum = 0
    base_ratio_diff_squared_sum = 0
    for oranism, values in first_tool_relative_abundances.items():
        read_ratio_diff = values["read_ratio"] - second_tool_relative_abundances[oranism]["read_ratio"]
        bases_ratio_diff = values["bases_ratio"] - second_tool_relative_abundances[oranism]["bases_ratio"]
        read_ratio_diff_squared_sum += read_ratio_diff**2
        base_ratio_diff_squared_sum += bases_ratio_diff**2
    return {
        'read_ratio': read_ratio_diff_squared_sum**0.5,
        'bases_ratio': base_ratio_diff_squared_sum**0.5
    }


def latex_escape(string):
    return string.replace('_', '\_')


def write_latex_table(relative_abundances, output_path, data_source='read_ratio'):
    #relative_abundances has toolnames as keys
    #table will have tools in rows and datasets in columns
    organism_count = len(relative_abundances['true_mappings'])
    for tool, organisms in relative_abundances.items():
        assert len(organisms) == organism_count
    
    renamed_relative_abundances = {}
    for tool, organisms in relative_abundances.items():
        escaped_organisms = {}
        for organism in organisms:
            escaped_organisms[latex_escape(organism_name_to_prety_name(organism))] = organisms[organism]
        renamed_relative_abundances[tool_name_to_pretty_name(tool)] = escaped_organisms

    ordered_tools = [
        'Ground Truth',
        'minimap2',
        'Uncalled',
        'Sigmap',
        'RawHash',
        'RawAlign'
    ]

    ordered_organisms = [
        'SARS-CoV-2',
        'E.coli',
        'Yeast',
        'Green Algae',
        'Human',
        'Distance'
    ]

    with open(output_path, 'w') as f:
        n = '\n'
        f.write(r'\begin{tabular}{l' + ' r'*organism_count + '}' + n)
        f.write(r'\hline' + n)
        f.write(
            ' & '.join(
                f'\\textbf{{{header}}}'
                for header in ['Tool'] + ordered_organisms
            )
        )
        f.write(r'\\' + n)
        #f.write(r'Tool & ' + ' & '.join(ordered_organisms) + r'\\' + n)

        best_distance = min([renamed_relative_abundances[tool]['Distance'][data_source] for tool in ordered_tools if tool != 'Ground Truth'])

        f.write(r'\hline' + n)
        for tool in ordered_tools:
            organisms = renamed_relative_abundances[tool]
            f.write(tool + ' & ')
            formatted_numbers = [
                (
                    f'{organisms[organism][data_source]:,.3f}'
                    if organisms[organism][data_source] != best_distance else
                    f'\\textbf{{{organisms[organism][data_source]:,.3f}}}'
                )
                if organism != 'Distance' or tool != 'Ground Truth'
                else '-'
                for organism in ordered_organisms
            ]
            f.write(' & '.join(formatted_numbers))
            f.write(r'\\' + n)
            #f.write(r'\hline' + n)
        f.write(r'\hline' + n)
        f.write(r'\end{tabular}' + n)


def tool_name_to_pretty_name(tool_name):
    if 'rawalign' in tool_name:
        return 'RawAlign'
    elif 'uncalled' in tool_name:
        return 'Uncalled'
    elif 'rawhash' in tool_name:
        return 'RawHash'
    elif 'sigmap' in tool_name:
        return 'Sigmap'
    elif 'true_mappings' in tool_name:
        return 'minimap2'
    elif 'fastq' in tool_name:
        return 'Ground Truth'
    else:
        return tool_name


def organism_name_to_prety_name(organism_name):
    if 'covid' in organism_name:
        return 'SARS-CoV-2'
    elif 'ecoli' in organism_name:
        return 'E.coli'
    elif 'yeast' in organism_name:
        return 'Yeast'
    elif 'green_algae' in organism_name:
        return 'Green Algae'
    elif 'human' in organism_name:
        return 'Human'
    elif 'distance' in organism_name:
        return 'Distance'
    else:
        return organism_name


if __name__ == '__main__':
    results = get_output_results(pathlib.Path('test/evaluation/relative_abundance/comparison/'))
    relative_abundances = parse_relative_abundances(results)

    ground_truth = relative_abundances['fastq']
    for tool, organisms in relative_abundances.items():
        organisms['euclidean_distance'] = get_euclidean_distance(organisms, ground_truth)

    write_latex_table(relative_abundances, pathlib.Path('paperfigures/read_ratio_relative_abundances.tex'), data_source='read_ratio')
    write_latex_table(relative_abundances, pathlib.Path('paperfigures/bases_ratio_relative_abundances.tex'), data_source='bases_ratio')
