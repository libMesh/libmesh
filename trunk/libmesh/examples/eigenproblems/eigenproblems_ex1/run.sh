#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex1

options="-n 5"
run_example "$example_name" "$options"
