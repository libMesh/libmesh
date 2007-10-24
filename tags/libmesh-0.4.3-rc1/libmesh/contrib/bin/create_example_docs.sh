#!/bin/bash
#
# This shell script is a wrapper for the script ex2html.sh
# also located in this directory.  This script calls ex2html.sh
# for all of the example programs.

# How many examples are there?
n_examples=`ls -l ../../examples/ | grep ex | wc -l`

for i in $(seq 1 $n_examples); do
    echo "Processing example $i";
    ./ex2html.sh ../../examples/ex$i/ex$i.C "make --no-print-directory -C ../../examples/ex$i run"
done