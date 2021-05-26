#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=reduced_basis_ex1

example_dir=examples/reduced_basis/$example_name

options="-online_mode 0"
run_example "$example_name" "$options"

options="-online_mode 1"
run_example "$example_name" "$options"

# Not benchmarking - 95% of our time is in the offline linear solve
