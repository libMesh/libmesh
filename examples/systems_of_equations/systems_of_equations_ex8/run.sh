#!/bin/sh

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=systems_of_equations_ex8

options=""

run_example "$example_name" "$options"

benchmark_example 1 "$example_name" "n_refinements=1"
