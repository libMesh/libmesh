#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex6

options="-ksp_type cg -pc_type bjacobi -sub_pc_type icc"

run_example "$example_name" "$options"
