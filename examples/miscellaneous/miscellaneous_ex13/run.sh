#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex13

options="-pc_type jacobi -ksp_type cg"

run_example "$example_name" "$options"

options="-distributed_load 1 -pc_type jacobi -ksp_type cg"

run_example "$example_name" "$options"

# No benchmarks here - I assume it's as efficient as its brother
# miscellaneous_ex12
