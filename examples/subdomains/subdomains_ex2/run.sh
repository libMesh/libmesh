#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=subdomains_ex2

options="-d 1 -n 20"
run_example "$example_name" "$options"

options="-d 2 -n 15"
run_example "$example_name" "$options"

options="-d 3 -n 6"
run_example "$example_name" "$options"

# Only benchmarking 3D because 1-D and 2-D spend 95% of their time in
# the linear solver
benchmark_example 1 "$example_name" "-d 3 -n 50"
