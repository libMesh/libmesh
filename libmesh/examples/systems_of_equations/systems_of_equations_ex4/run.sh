#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex4

options="-ksp_type cg"

run_example "$example_name" "$options"
