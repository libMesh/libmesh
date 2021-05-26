#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex1

run_example "$example_name" 

# This is borderline, already ~85% in the linear solve.
benchmark_example 1 "$example_name" "-n_elem 300"
