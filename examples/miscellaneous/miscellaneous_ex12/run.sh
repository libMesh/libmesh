#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex12

options="-pc_type lu -pc_factor_mat_solver_package mumps"

run_example "$example_name" "$options"
