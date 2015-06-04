#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=fem_system_ex3

example_dir=examples/fem_system/$example_name

options="-ksp_type cg -pc_type bjacobi"

run_example "$example_name" "$options"
