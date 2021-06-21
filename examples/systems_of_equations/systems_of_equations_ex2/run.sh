#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex2

run_example "$example_name" 

# When we allow more complete linear solves, we finish faster ... but
# we also end up with 93% of our time in the linear solver.
#benchmark_example 1 "$example_name" "-n_elem 40 -max_iter 2500"
