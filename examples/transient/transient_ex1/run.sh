#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=transient_ex1

run_example "$example_name" 

benchmark_example 1 "$example_name" "n_refinements=8"
