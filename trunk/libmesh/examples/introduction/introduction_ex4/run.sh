#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex4

options="-d 1 -n 20"
run_example "$example_name" "$options"

options="-d 2 -n 15"
run_example "$example_name" "$options"

options="-d 3 -n  6"
run_example "$example_name" "$options"
