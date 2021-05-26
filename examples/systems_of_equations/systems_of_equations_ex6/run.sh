#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex6

options=""

run_example "$example_name" "$options"

benchmark_example 1 "$example_name" "$options -nx 160 -ny 40 -nz 20"
