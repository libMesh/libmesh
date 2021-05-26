#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex3

run_example "$example_name" 

benchmark_example 1 "$example_name" "-n_elem 80"
