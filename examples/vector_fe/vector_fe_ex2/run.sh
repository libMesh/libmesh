#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex2

run_example "$example_name" 

# compute_shape_functions is way too expensive here?
benchmark_example 1 "$example_name" "grid_size=40"
