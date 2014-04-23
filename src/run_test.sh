#!/bin/bash

./generate_data -f $1 $2 testfile.txt
./remove_nulls.sh testfile.txt
./NavSeq -compress testfile.txt testfile-compressed.txt
./NavSeq -access testfile-compressed.txt 142
./NavSeq -rank testfile-compressed.txt 981 0