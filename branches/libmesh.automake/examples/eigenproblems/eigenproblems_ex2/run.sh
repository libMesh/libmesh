#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex2

options="-n 5 -eps_type lapack"
run_example "$example_name" "$options"
