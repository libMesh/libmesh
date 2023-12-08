#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex1

run_example "$example_name" -o FIRST

run_example "$example_name" -o SECOND

run_example "$example_name" -o SECOND -f HIERARCHIC_VEC

run_example "$example_name" -o THIRD -f HIERARCHIC_VEC

# Higher orders run into roundoff error limitations so fast we might
# as well not bother to refine them far
run_example "$example_name" -o FOURTH -nx 13 -ny 13 -f HIERARCHIC_VEC

run_example "$example_name" -o FIFTH -nx 6 -ny 6 -f HIERARCHIC_VEC

# No benchmark here - this spends 95+% of time in the linear solve
# benchmark_example 1 "$example_name" "-nx 200 -ny 200"
