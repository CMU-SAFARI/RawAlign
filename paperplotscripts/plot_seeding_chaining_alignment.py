import re
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def parse_tags(line):
    tags = {}
    tag_pattern = r'(\S+):(\S+):(\S+)'  # Regex pattern to match tag_key:type

    # Find all occurrences of tag_key:type using regex
    matches = re.findall(tag_pattern, line)
    for key, valuetype, value in matches:
        if key == 'aln':
            # Parse the tuples into arrays
            tuples = re.findall(r'\(([^)]+)\)', value)
            parsed_tuples = []
            for t in tuples:
                t_values = t.split(',')
                x, y, z = int(t_values[0]), int(t_values[1]), float(t_values[2])
                parsed_tuples.append((x, y, z))
            tags[key] = parsed_tuples
        elif key == 'anchors':
            # Parse the tuples into arrays
            tuples = re.findall(r'\(([^)]+)\)', value)
            parsed_tuples = []
            for t in tuples:
                t_values = t.split(',')
                x, y = int(t_values[0]), int(t_values[1])
                parsed_tuples.append((x, y))
            tags[key] = parsed_tuples
        else:
            tags[key] = value
    return tags

def parse_paf(input_file, limit=None):
    res = {}
    with open(input_file, 'r') as file:
        i = 0
        for line in file:
            line = line.strip()
            if line:
                default_fields = list(line.split('\t')[:12])
                tag_data = parse_tags(line)
                readname = default_fields[0]
                strand = default_fields[4]
                refname = default_fields[5]
                if 'aln' not in tag_data or 'anchors' not in tag_data:
                    continue
                res[(readname, refname, strand)] = {
                    'alignment': tag_data['aln'],
                    'chain': tag_data['anchors'],
                }

                i += 1
                if i % 1000 == 0:
                    print(f'Parsed {i} lines from PAF')

                if limit is not None and i >= limit:
                    break
    return res

def plot_aln(ax, aln, label, color='red'):
    xs = []
    ys = []
    for x, y, z in aln:
        xs.append(x)
        ys.append(y)
    ax.plot(xs, ys, label=label, color=color, linewidth=2)

def plot_anchors(ax, anchors, label, color='red'):
    xs = []
    ys = []
    for x, y in anchors:
        xs.append(x)
        ys.append(y)

    ax.scatter(xs, ys, label=label, color=color)

def plot_chain(ax, anchors, label, color='red'):
    xs = []
    ys = []
    for x, y in anchors:
        xs.append(x)
        ys.append(y)

    #ax.scatter(xs, ys, label=label, color=color)
    ax.plot(xs, ys, label=label, color=color, linestyle='--', marker='o', linewidth=2)

def parse_anchor_log(input_file, line_limit=None, filter_list=None):
    anchor_line_pattern = r'readname=([\S]+)\srefname=(\S+)\sstrand=(\d)\sanchors=(\(\d+,\d+\)(?:\(\d+,\d+\))*)'
    readref_pattern = r'readname=(\S+)\srefname=(\S+)\sstrand=(\d)'

    res = {}
    with open(input_file, 'r') as file:
        i = 0
        for line in file:
            if line_limit is not None and i >= line_limit:
                break

            parsed_log_frequency = 100
            if (i) % parsed_log_frequency == 0:
                print(f'Parsed {i} lines from log')

            i += 1

            if filter_list is not None:
                match = re.match(readref_pattern, line[:1000])

                if not match:
                    continue
                readname = match.group(1)
                refname = match.group(2)
                strand = '-' if int(match.group(3))==1 else '+'

                if (readname, refname, strand) not in filter_list:
                    continue

            match = re.match(anchor_line_pattern, line)
            if not match:
                #print(line)
                continue

            readname = match.group(1)
            refname = match.group(2)
            strand = '-' if int(match.group(3))==1 else '+'
            anchors_str = match.group(4)

            # Split the anchors string and parse into a list of tuples
            anchors_pairs = re.findall(r'\((\d+,\d+)\)', anchors_str)
            anchors = [tuple(map(int, pair.split(','))) for pair in anchors_pairs]

            if (readname, refname, strand) not in res:
                res[(readname, refname, strand)] = []
            res[(readname, refname, strand)]+=anchors

    return res

def plot(readname, refname, strand, anchors, chain, alignment, figurespath):
    fig, axs = plt.subplots(1, 3, layout="constrained")
    fig.set_size_inches(8, 4)

    #fig.suptitle(f'Read {readname}')

    for ax in axs:
        ax.set_aspect('equal')
        ax.grid(linewidth=1)
        ax.set_axisbelow(True)

    plot_anchors(axs[0], anchors, 'Anchors')

    plot_anchors(axs[1], set(anchors) - set(chain), 'Anchors', color=(0.5, 0.5, 0.5))
    plot_chain(axs[1], chain, 'Chain')

    plot_chain(axs[2], chain, 'Chain', color=(0.5, 0.5, 0.5))
    plot_aln(axs[2], alignment, 'Alignment')
    axs[0].set_title('1. Seeding',     fontweight='bold', size=14, pad=0)
    axs[1].set_title('2. Chaining',    fontweight='bold', size=14, pad=0)
    axs[2].set_title('3. Alignment',   fontweight='bold', size=14, pad=0)

    min_plotted_x = min([x for x, y, z in alignment])
    max_plotted_x = max([x for x, y, z in alignment])
    min_plotted_y = min([y for x, y, z in alignment])
    max_plotted_y = max([y for x, y, z in alignment])
    plot_width = max(max_plotted_x - min_plotted_x, max_plotted_y - min_plotted_y)
    for ax in axs:
        ax.set_xlim(min_plotted_x - plot_width*0.05, min_plotted_x + plot_width + plot_width*0.05)
        ax.set_ylim(min_plotted_y - plot_width*0.05, min_plotted_y + plot_width + plot_width*0.05)
        ax.set_xlabel('Read Pos.', fontsize=12, fontweight='bold', labelpad=0)
        ax.legend(fontsize=10, framealpha=1, edgecolor='black', facecolor='white',
                  loc='upper right', fancybox=False,
                  bbox_to_anchor=(1.0, 1.0), borderaxespad=0,
                  labelspacing=0.2,
                  handletextpad=0.2, handlelength=1.0,
                  )#borderpad=0.2,)
        ax.xaxis.set_tick_params(labelsize=10, pad=1)
        ax.yaxis.set_tick_params(labelsize=10, pad=1)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=100))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=100))
        ax.invert_yaxis()
        for ticklabel in ax.xaxis.get_ticklabels():
            ticklabel.set_fontweight('bold')

    for ax in axs[1:]:
        ax.set_yticklabels([])
    axs[0].set_ylabel('Reference Pos.', fontsize=12, fontweight='bold', labelpad=0)
    y_offset = min_plotted_y - (min_plotted_y % 1000)
    axs[0].yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=y_offset))
    axs[0].yaxis.get_offset_text().set_fontsize(10)
    axs[0].yaxis.get_offset_text().set_fontweight('bold')
    axs[0].yaxis.get_offset_text().set_x(-0.1)
    for ticklabel in axs[0].yaxis.get_ticklabels():
        ticklabel.set_fontweight('bold')

    fig.savefig(figurespath / f'{readname}.pdf', dpi=480, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

if __name__ == "__main__":
    figspath = pathlib.Path('paperfigures/seeding_chaining_alignment')
    figspath.mkdir(exist_ok=True)

    paf_file_path = pathlib.Path('test/evaluation/read_mapping/d2_ecoli_r94/outdir/d2_ecoli_r94_rawalign_sensitive_dtwevaluatechains_loganchors_outputchains_dtwoutputcigar_dtwborderconstraint_sparse_dtwfillmethod_banded=0.10_dtwmatchbonus_0.4_dtwminscore_20.0_stopminanchor_2.paf')    
    paf = parse_paf(paf_file_path, limit=100)

    anchor_log_path = pathlib.Path('test/evaluation/read_mapping/d2_ecoli_r94/outdir/d2_ecoli_r94_rawalign_sensitive_dtwevaluatechains_loganchors_outputchains_dtwoutputcigar_dtwborderconstraint_sparse_dtwfillmethod_banded=0.10_dtwmatchbonus_0.4_dtwminscore_20.0_stopminanchor_2.err')
    anchors = parse_anchor_log(anchor_log_path, line_limit=10000, filter_list=list(paf.keys()))

    keys_in_both = set(paf.keys()) & set(anchors.keys())

    merged = {}
    for key in keys_in_both:
        merged[key] = {
            'alignment': paf[key]['alignment'],
            'chain': paf[key]['chain'],
            'anchors': anchors[key],
        }

    for key, value in merged.items():
        readname, refname, strand = key
        alignment = value['alignment']
        chain = value['chain']
        anchors = value['anchors']

        if max(chain, key=lambda x: x[0]) not in anchors: #skip pairs for which the log was not parsed fully
            print(f'Skipping {readname}')
            continue
        plot(readname, refname, strand, anchors, chain, alignment, figspath)

    #plt.show()
