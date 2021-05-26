#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex1

run_example "$example_name" 

# No benchmark here - this spends 95+% of time in the linear solve
# benchmark_example 1 "$example_name" "-nx 200 -ny 200"
