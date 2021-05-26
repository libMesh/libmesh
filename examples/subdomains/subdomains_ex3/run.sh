#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=subdomains_ex3

options="-d 2"
run_example "$example_name" "$options"

options="-d 3"
run_example "$example_name" "$options"

# No benchmarks here - if I bump n_refinements above 4 I get "must
# have 3 input vertices" errors
# benchmark_example 1 "$example_name" "-d 2 -n_refinements 6"
# benchmark_example 1 "$example_name" "-d 3 -n_refinements 6"
