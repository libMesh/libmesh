#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=reduced_basis_ex3

example_dir=examples/reduced_basis/$example_name

options="-online_mode 0"
run_example "$example_name" "$options"
benchmark_example 1 "$example_name" "$options n_elem=150"

options="-online_mode 1"
run_example "$example_name" "$options"
benchmark_example 1 "$example_name" "$options n_elem=150"
