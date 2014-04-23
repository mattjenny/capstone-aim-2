#!/bin/bash

for i in {2..8}
do
    Num=$((10**$i))
    for j in {1..9}
    do
	echo "For $Num variants and $j patients:"
	./run_test.sh $Num $j
    done
done

Num=$((10**9))
for j in {1..4}
do
    echo "For $Num variants and $j patients:"
    ./run_test.sh $Num $j
done