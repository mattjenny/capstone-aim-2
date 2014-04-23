#!/bin/bash

NUM_ITERATIONS=$1
FLAG=$2
FILENAME=$3
VARIANT_ARG=$4
POPULATION_SIZE=$5
FREQUENCY_CUTOFF=0.5

COUNT=0
while [ $COUNT -lt $NUM_ITERATIONS ]; do
	../CapstoneQuery $FLAG $FILENAME $VARIANT_ARG $POPULATION_SIZE $FREQUENCY_CUTOFF >temp.out
	COUNT=$[$COUNT+1]
done