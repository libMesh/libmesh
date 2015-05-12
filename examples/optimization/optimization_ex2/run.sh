#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=optimization_ex2

example_dir=examples/optimization/$example_name

options="-tao_monitor -tao_view -tao_type ipm -pc_type jacobi -ksp_max_it 100"

# Prevent the environment from running this example in parallel or
# changing the command line arguments.
LIBMESH_OPTIONS=
LIBMESH_RUN=
run_example "$example_name" "$options"
