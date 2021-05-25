#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex12

options="-pc_type jacobi -ksp_type cg"

run_example "$example_name" "$options"

options="-distributed_load 1 -pc_type jacobi -ksp_type cg"

run_example "$example_name" "$options"

# No benchmarks here - if we scale up a little we spend 90%+ of our
# time in the solver
