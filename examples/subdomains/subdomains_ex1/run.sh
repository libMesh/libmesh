#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=subdomains_ex1

options="-d 1 -n 40"
run_example "$example_name" "$options"

options="-d 2 -n 30"
run_example "$example_name" "$options"

options="-d 3 -n 15"
run_example "$example_name" "$options"

benchmark_example 1 "$example_name" "-d 1 -n 1000000"
benchmark_example 1 "$example_name" "-d 2 -n 1000"
benchmark_example 1 "$example_name" "-d 3 -n 120"
