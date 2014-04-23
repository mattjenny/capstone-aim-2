#!/bin/bash

tr < $1 -d '\000' > nonulls.txt
rm $1
mv nonulls.txt $1