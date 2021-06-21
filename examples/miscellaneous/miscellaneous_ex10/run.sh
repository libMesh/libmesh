#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex10

options="-n 2"

run_example "$example_name" "$options"

benchmark_example 1 "$example_name" "-n 7"
