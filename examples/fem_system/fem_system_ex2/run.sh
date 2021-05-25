#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=fem_system_ex2

example_dir=examples/fem_system/$example_name

run_example "$example_name"

# Skipping benchmarking until we fix up our bash scripts and get this
# syntax working:
# benchmark_example 1 "$example_name" mesh/generation/num_elem='1 1 1' testvar=1
