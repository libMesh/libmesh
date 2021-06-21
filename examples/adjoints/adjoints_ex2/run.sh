#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex2

example_dir=examples/adjoints/$example_name

run_example "$example_name"

# Benchmark parameters
benchmark_example 1 "$example_name" coarserefinements=4
