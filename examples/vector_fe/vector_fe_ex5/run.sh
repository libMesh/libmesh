#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex5

run_example "$example_name"

# No benchmark here - this spends 88% of time in PETSc?
# benchmark_example 1 "$example_name" "nx=250 ny=250"
