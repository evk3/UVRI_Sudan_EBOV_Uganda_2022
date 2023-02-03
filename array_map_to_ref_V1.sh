#!/bin/bash -l

#The purpose of this wrapper script is to submit
#the map_and_dedup script using all files in file_names.txt


files=$(awk 'END{print NR}' ./file_names.txt)
echo $files

#Submit kraken script as an array
#Submit script 1 to length of patient_array number of times.
qsub -t 1-$files:2 \
	-N EBOV \
	./EBOV_iterative_map.sh \

echo "Wrapper script done"
