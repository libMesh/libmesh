#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex9

options=""

run_example "$example_name" "$options"

# 81% in the linear solve, kind of borderline
benchmark_example 1 "$example_name" "n_refinements=1"
