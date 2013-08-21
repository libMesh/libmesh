#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex4

options="element_type=HEX20 -pc_type jacobi"
run_example "$example_name" "$options"

options="element_type=HEX27 -pc_type jacobi"
run_example "$example_name" "$options"
