#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex1

run_example "$example_name"

benchmark_example 1 "$example_name" "nx=100 ny=100 nz=100"
