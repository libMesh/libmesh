#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex12

options="-pc_type jacobi -ksp_type cg"

run_example "$example_name" "$options"
