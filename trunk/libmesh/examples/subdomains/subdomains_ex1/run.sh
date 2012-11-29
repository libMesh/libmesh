#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=subdomains_ex1

options="-d 1 -n 40"
run_example "$example_name" "$options"

options="-d 2 -n 30"
run_example "$example_name" "$options"

options="-d 3 -n 15"
run_example "$example_name" "$options"
