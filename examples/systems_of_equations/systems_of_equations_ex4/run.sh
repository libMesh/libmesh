#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex4

options="-ksp_type cg"

run_example "$example_name" "$options"

# No benchmark here - this spends 98% of time in the linear solve
# benchmark_example 1 "$example_name" "$options -nx 500 -ny 100"
