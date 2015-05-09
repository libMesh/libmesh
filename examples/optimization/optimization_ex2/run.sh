#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=optimization_ex2

example_dir=examples/optimization/$example_name

options="-tao_monitor -tao_view -tao_type ipm -pc_type jacobi -ksp_max_it 100"
run_example "$example_name" "$options"
