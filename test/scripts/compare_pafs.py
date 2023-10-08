import os
import sys
import argparse
import subprocess
import fileinput
from pathlib import Path

from statistics import median
from statistics import mean

def analyze_paf_file(toolname, paf_type, paf_path):
	if 'Uncalled' in paf_type:
		tp, fp, fn, tn, time_per_read, maplast_pos, umaplast_pos \
			= analyze_uncalled(paf_path)
		maplast_chunk = []
		umaplast_chunk = []
	elif 'Sigmap' in paf_type:
		tp, fp, fn, tn, time_per_mapped_read, time_per_unmapped_read, maplast_chunk, umaplast_chunk \
			= analyze_sigmap(paf_path)
		maplast_pos = []
		umaplast_pos = []
		time_per_read = time_per_mapped_read + time_per_unmapped_read
	elif 'RawHash' in paf_type:
		tp, fp, fn, tn, time_per_mapped_read, time_per_unmapped_read, maplast_pos, umaplast_pos, maplast_chunk, umaplast_chunk \
			= analyze_rawhash(paf_path)
		time_per_read = time_per_mapped_read + time_per_unmapped_read
	elif 'RawAlign' in paf_type:
		tp, fp, fn, tn, time_per_mapped_read, time_per_unmapped_read, maplast_pos, umaplast_pos, maplast_chunk, umaplast_chunk \
			= analyze_rawalign(paf_path)
		time_per_read = time_per_mapped_read + time_per_unmapped_read
		maplast_pos = [] #temporarily disabled since the last round of experiments had an erroneous read length in the rawalign binary
		umaplast_pos = []
	else:
		print("Unknown paf type " + paf_type)
		sys.exit(1)

	print(f"{toolname} TP: " + str(tp))
	print(f"{toolname} FP: " + str(fp))
	print(f"{toolname} FN: " + str(fn))
	print(f"{toolname} TN: " + str(tn))
	precision = tp / (tp + fp) if (tp + fp) > 0 else 0
	print(f"{toolname} precision: " + str(precision))
	recall = tp / (tp + fn) if (tp + fn) > 0 else 0
	print(f"{toolname} recall: " + str(recall))
	f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
	print(f"{toolname} F-1 score: " + str(f1))
	#initie if no time_per_read
	mean_time = mean(time_per_read) if len(time_per_read) > 0 else float('inf')
	print(f"{toolname} Mean time per read : " + str(mean_time))
	median_time = median(time_per_read) if len(time_per_read) > 0 else float('inf')
	print(f"{toolname} Median time per read : " + str(median_time))
	mean_maplast_pos = mean(maplast_pos) if len(maplast_pos) > 0 else float('inf')
	print(f"{toolname} Mean (only mapped) # of sequenced bases per read : " + str(mean_maplast_pos))
	mean_umaplast_pos = mean(umaplast_pos) if len(umaplast_pos) > 0 else float('inf')
	print(f"{toolname} Mean (only unmapped) # of sequenced bases per read : " + str(mean_umaplast_pos))
	mean_last_pos = mean(maplast_pos + umaplast_pos) if len(maplast_pos + umaplast_pos) > 0 else float('inf')
	print(f"{toolname} Mean # of sequenced bases per read : " + str(mean_last_pos))
	mean_maplast_chunk = mean(maplast_chunk) if len(maplast_chunk) > 0 else float('inf')
	print(f"{toolname} Mean (only mapped) # of sequenced chunks per read : " + str(mean_maplast_chunk))
	mean_umaplast_chunk = mean(umaplast_chunk) if len(umaplast_chunk) > 0 else float('inf')
	print(f"{toolname} Mean (only unmapped) # of sequenced chunks per read : " + str(mean_umaplast_chunk))
	mean_last_chunk = mean(maplast_chunk + umaplast_chunk) if len(maplast_chunk + umaplast_chunk) > 0 else float('inf')
	print(f"{toolname} Mean # of sequenced chunks per read : " + str(mean_last_chunk))
	print(f"#Done with {toolname}\n")


def analyze_uncalled(paf_path):
	tp = 0
	fp = 0
	fn = 0
	tn = 0
	time_per_read = []
	maplast_pos = []
	umaplast_pos = []

	for line in paf_path.open():
		cols = line.rstrip().split()
		mt = float(cols[14].split(":")[2])
		lastpos = int(cols[1])
		is_mapped = cols[2] != '*'

		if (is_mapped):
			maplast_pos.append(lastpos)
		else:
			umaplast_pos.append(lastpos)

		if (cols[15].split(":")[2] != 'na'):
			time_per_read.append(mt)
		if (cols[15].split(":")[2] == 'tp'):
			tp += 1
		if (cols[15].split(":")[2] == 'fp' or cols[15].split(":")[2] == 'na'):
			fp += 1
		if (cols[15].split(":")[2] == 'fn'):
			fn += 1
		if (cols[15].split(":")[2] == 'tn'):
			tn += 1

	return tp, fp, fn, tn, time_per_read, maplast_pos, umaplast_pos


def analyze_sigmap(paf_path):
	tp = 0
	fp = 0
	fn = 0
	tn = 0

	time_per_mapped_read = []
	time_per_unmapped_read = []
	maplast_chunk = []
	umaplast_chunk = []

	for line in paf_path.open():
		cols = line.rstrip().split()
		if (len(cols) == 24):
			mt = float(cols[12].split(":")[2])
			if (cols[23].split(":")[2] != 'na'):
				time_per_unmapped_read.append(mt)
			chunk = int(cols[13].split(":")[2])
			is_mapped = cols[2] != '*'
			if (is_mapped):
				maplast_chunk.append(chunk)
			else:
				umaplast_chunk.append(chunk)
			cm = int(cols[15].split(":")[2])
			nc = int(cols[16].split(":")[2])
			s1 = float(cols[17].split(":")[2])
			s2 = float(cols[18].split(":")[2])
			sm = float(cols[19].split(":")[2])
			ad = float(cols[20].split(":")[2])
			at = float(cols[21].split(":")[2])
			aq = float(cols[22].split(":")[2])
			if (cols[23].split(":")[2] == 'tp'):
				tp += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[23].split(":")[2] == 'fp' or cols[23].split(":")[2] == 'na'):
				fp += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[23].split(":")[2] == 'fn'):
				fn += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[23].split(":")[2] == 'tn'):
				tn += 1
				time_per_mapped_read.append(mt / chunk)
		if (len(cols) == 15):
			mt = float(cols[12].split(":")[2])
			if (cols[14].split(":")[2] != 'na'):
				time_per_unmapped_read.append(mt)
			if (cols[14].split(":")[2] == 'fn'):
				fn += 1
			if (cols[14].split(":")[2] == 'tn'):
				tn += 1

	return tp, fp, fn, tn, time_per_mapped_read, time_per_unmapped_read, maplast_chunk, umaplast_chunk


def analyze_rawhash(paf_path):
	tp = 0
	fp = 0
	fn = 0
	tn = 0

	time_per_mapped_read = []
	time_per_unmapped_read = []
	maplast_pos = []
	umaplast_pos = []
	maplast_chunk = []
	umaplast_chunk = []

	for line in paf_path.open():
		cols = line.rstrip().split()
		if (len(cols) == 23):
			mt = float(cols[12].split(":")[2])
			if (cols[22].split(":")[2] != 'na'):
				time_per_unmapped_read.append(mt)
			lastpos = int(cols[1])
			chunk = int(cols[13].split(":")[2])
			is_mapped = cols[2] != '*'
			if (is_mapped):
				maplast_pos.append(lastpos)
				maplast_chunk.append(chunk)
			else:
				umaplast_pos.append(lastpos)
				umaplast_chunk.append(chunk)
			cm = int(cols[15].split(":")[2])
			nc = int(cols[16].split(":")[2])
			s1 = float(cols[17].split(":")[2])
			s2 = float(cols[18].split(":")[2])
			sm = float(cols[19].split(":")[2])
			at = float(cols[20].split(":")[2])
			aq = float(cols[21].split(":")[2])
			if (cols[22].split(":")[2] == 'tp'):
				tp += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[22].split(":")[2] == 'fp' or cols[22].split(":")[2] == 'na'):
				fp += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[22].split(":")[2] == 'fn'):
				fn += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[22].split(":")[2] == 'tn'):
				tn += 1
				time_per_mapped_read.append(mt / chunk)
		if (len(cols) == 15):
			mt = float(cols[12].split(":")[2])
			if (cols[14].split(":")[2] != 'na'):
				time_per_unmapped_read.append(mt)
			if (cols[14].split(":")[2] == 'fn'):
				fn += 1
			if (cols[14].split(":")[2] == 'tn'):
				tn += 1

	return tp, fp, fn, tn, time_per_mapped_read, time_per_unmapped_read, maplast_pos, umaplast_pos, maplast_chunk, umaplast_chunk


def analyze_rawalign(paf_path):
	tp = 0
	fp = 0
	fn = 0
	tn = 0

	time_per_mapped_read = []
	time_per_unmapped_read = []
	maplast_pos = []
	umaplast_pos = []
	maplast_chunk = []
	umaplast_chunk = []
	for line in paf_path.open():
		cols = line.rstrip().split()
		if (len(cols) == 23):
			mt = float(cols[12].split(":")[2])
			if (cols[22].split(":")[2] != 'na'):
				time_per_unmapped_read.append(mt)
			lastpos = int(cols[1])
			chunk = int(cols[13].split(":")[2])
			is_mapped = cols[2] != '*'
			if (is_mapped):
				maplast_pos.append(lastpos)
				maplast_chunk.append(chunk)
			else:
				umaplast_pos.append(lastpos)
				umaplast_chunk.append(chunk)
			cm = int(cols[15].split(":")[2])
			nc = int(cols[16].split(":")[2])
			s1 = float(cols[17].split(":")[2])
			s2 = float(cols[18].split(":")[2])
			sm = float(cols[19].split(":")[2])
			at = float(cols[20].split(":")[2])
			aq = float(cols[21].split(":")[2])
			if (cols[22].split(":")[2] == 'tp'):
				tp += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[22].split(":")[2] == 'fp' or cols[22].split(":")[2] == 'na'):
				fp += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[22].split(":")[2] == 'fn'):
				fn += 1
				time_per_mapped_read.append(mt / chunk)
			if (cols[22].split(":")[2] == 'tn'):
				tn += 1
				time_per_mapped_read.append(mt / chunk)
		if (len(cols) == 15):
			mt = float(cols[12].split(":")[2])
			if (cols[14].split(":")[2] != 'na'):
				time_per_unmapped_read.append(mt)
			if (cols[14].split(":")[2] == 'fn'):
				fn += 1
			if (cols[14].split(":")[2] == 'tn'):
				tn += 1

	return tp, fp, fn, tn, time_per_mapped_read, time_per_unmapped_read, maplast_pos, umaplast_pos, maplast_chunk, umaplast_chunk


def longest_common_prefix(strings):
    if not strings:
        return ""

    prefix = strings[0]

    for string in strings[1:]:
        while string[:len(prefix)] != prefix:
            prefix = prefix[:-1]
            if not prefix:
                return ""
    
    return prefix


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Compare PAFs')
	parser.add_argument('pafs', metavar='PAF', type=str, nargs='+', help='annotated PAF files to compare in the format NAME=PATH')
	args = parser.parse_args()

	#ensure that all files exist
	pafs = []
	for namedpaf in args.pafs:
		paf_type, *paf = namedpaf.split('=')
		paf = '='.join(paf) #rejoin in case there were '=' in the paf path

		paf = Path(paf)
		if not paf.is_file():	
			print("PAF file " + str(paf) + " does not exist")
			sys.exit(1)
		pafs.append((paf_type, paf))

	#find shared prefix between paf filenames
	basenames = [paf.name for (paf_type, paf) in pafs]
	prefix = longest_common_prefix(basenames)

	for paf_type, paf in pafs:
		toolname = paf.name.replace(prefix, "").replace("_ann.paf", "")
		print(f"Analyzing {paf} (aka {toolname}) as {paf_type}")
		analyze_paf_file(toolname, paf_type, paf)
