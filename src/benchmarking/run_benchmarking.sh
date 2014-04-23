#!/bin/bash

TIMEFORMAT=%R

NUM_ITERATIONS=$1
UNCOMPRESSED_SUFFIX=".txt"
COMPRESSED_SUFFIX="_out.txt"

echo "USize,CSize,pop_size,num_variants,c1_time,u1_time,c2_time,u2_time"
for i in {2..9}
do
	filesize=$((10**$i))
	file_prefix="benchfile_$filesize"
	compressed_filename="$file_prefix$COMPRESSED_SUFFIX"
	uncompressed_filename="$file_prefix$UNCOMPRESSED_SUFFIX"
	../GenerateData -f $filesize 1 $uncompressed_filename
	./remove_nulls.sh $uncompressed_filename
	../NavarroCompress $uncompressed_filename $compressed_filename
	compressed_filesize=$(wc -c $compressed_filename | awk '{ print $1 }')
	uncompressed_filesize=$(wc -c $uncompressed_filename | awk '{ print $1 }')
	max_dimension=`expr $i-1`
	for (( j=1; j<=$max_dimension; j++ ))
	do
		population_size=$((10**$j))
		num_variants=$((10**`expr $i - $j`))
		c1_time=$( { time ./run_trials.sh $NUM_ITERATIONS -c1 $compressed_filename 9 $population_size; } 2>&1 )
		u1_time=$( { time ./run_trials.sh $NUM_ITERATIONS -u1 $uncompressed_filename 9 $population_size; } 2>&1 )
		c2_time=$( { time ./run_trials.sh $NUM_ITERATIONS -c2 $compressed_filename $num_variants $population_size; } 2>&1 )
		u2_time=$( { time ./run_trials.sh $NUM_ITERATIONS -u2 $uncompressed_filename $num_variants $population_size; } 2>&1 )
		echo "$uncompressed_filesize,$compressed_filesize,$population_size,$num_variants,$c1_time,$u1_time,$c2_time,$u2_time"
	done
done