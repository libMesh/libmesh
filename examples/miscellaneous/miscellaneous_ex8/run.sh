#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex8

options=""

run_example "$example_name" "$options"

benchmark_example 1 "$example_name" "n_refinements=1 n_source_points=10000 n_target_points=100"
