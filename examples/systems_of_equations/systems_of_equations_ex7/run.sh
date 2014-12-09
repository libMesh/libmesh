#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex7

options="-ksp_type cg -pc_type bjacobi"

run_example "$example_name" "$options"
