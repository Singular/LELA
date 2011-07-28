#!/bin/bash

rm -rf output
mkdir output

echo "Running $* 100 times, putting output in ./output/test.output-0-99"

i=0

while test $i -lt 100; do
    $* output/test.output-$i && rm -f output/test.output-$i
    i=$[$i+1]
done
