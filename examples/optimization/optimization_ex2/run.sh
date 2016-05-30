#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=optimization_ex2

example_dir=examples/optimization/$example_name

options="-tao_monitor -tao_view -tao_type ipm -pc_type jacobi -ksp_max_it 100"

# This example is very sensitive to solver options
run_example_no_extra_options "$example_name" "$options"
