#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=fem_system_ex3

example_dir=examples/fem_system/$example_name

# Note: Too much ILU fails badly on this problem in single precision.
# We force simple Jacobi to be safe.

options="-ksp_type cg -pc_type jacobi"

run_example_no_extra_options "$example_name" "$options"
