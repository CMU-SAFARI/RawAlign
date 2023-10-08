#!/bin/bash


# This function takes a list of filenames as arguments and returns the
# longest common prefix of all of them.
# Only the basename is considered, so the directory path is ignored.
get_shared_filename_prefix() {
    local files=("$@")

    # Extract basenames and sort them
    basenames=()
    for file in "${files[@]}"; do
        basenames+=("$(basename "$file")")
    done

    IFS=$'\n' sorted_basenames=($(sort <<<"${basenames[*]}"))
    unset IFS

    # Find the common prefix by comparing the first and last basenames
    local first_basename="${sorted_basenames[0]}"
    local last_basename="${sorted_basenames[-1]}"

    local common_prefix=""

    for (( i=0; i<${#first_basename}; i++ )); do
        if [ "${first_basename:$i:1}" != "${last_basename:$i:1}" ]; then
            break
        fi
        common_prefix+="${first_basename:$i:1}"
    done

    echo "$common_prefix"
}

files=(../*/*.throughput)
prefix=$(get_shared_filename_prefix "${files[@]}")

echo "Throughput (mean/median)"
for i in `echo ../*/*.throughput`; do
    FILENAME=$(basename "${i%.throughput}");
	TOOLNAME="${FILENAME#${prefix}}";
	GREPRESULT=$(grep "BP per sec:" $i);
	if [ ! -z "$GREPRESULT" ]; then
		echo -n "${TOOLNAME} ";
		echo $GREPRESULT;
	fi;
done

echo;
echo "Mean time per read"
grep "Mean time per read" *.comparison

echo;
echo '(Indexing) Timing results:'
for i in `echo ../*/*index*.time`; do
    FILENAME=$(basename "${i%.time}");
	TOOLNAME="${FILENAME#${prefix}}";
	echo -n "${TOOLNAME} ";
	awk '{
		if(NR == 2){
			time = $NF
		}else if(NR == 3){
			time += $NF
			printf("CPU Time: %.2f\n", time)
		}
	}' $i;
done

echo;
echo '(Indexing) Memory usage results:'
for i in `echo ../*/*index*.time`; do
    FILENAME=$(basename "${i%.time}");
	TOOLNAME="${FILENAME#${prefix}}";
	echo -n "${TOOLNAME} ";
	awk '{
		if(NR == 10){
			printf("Memory (GB): %.2f\n", $NF/1000000)
		}
	}' $i;
done

echo;
echo '(Mapping) Timing results:'
for i in `echo ../*/*_map*.time`; do
    FILENAME=$(basename "${i%.time}");
	TOOLNAME="${FILENAME#${prefix}}";
	echo -n "${TOOLNAME} ";
	awk '{
		if(NR == 2){
			time = $NF
		}else if(NR == 3){
			time += $NF
			printf("CPU Time: %.2f\n", time)
		}
	}' $i;
done

echo;
echo '(Mapping) Memory usage results:'
for i in `echo ../*/*_map*.time`; do
    FILENAME=$(basename "${i%.time}");
	TOOLNAME="${FILENAME#${prefix}}";
	echo -n "${TOOLNAME} ";
	awk '{
		if(NR == 10){
			printf("Memory (GB): %.2f\n", $NF/1000000)
		}
	}' $i;
done

echo;
echo "Mean sequenced bases:"
grep "Mean # of sequenced bases per read :" *.comparison

echo;
echo "Mean sequenced chunks:"
grep "Mean # of sequenced chunks per read :" *.comparison

echo;
echo "Precision:"
grep "precision:" *.comparison

echo;
echo "Recall:"
grep "recall:" *.comparison

echo;
echo "F-1 score:"
grep "F-1 score:" *.comparison
