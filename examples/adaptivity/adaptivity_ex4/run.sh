#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex4

run_example "$example_name" "$options" approx_type=HERMITE
run_example "$example_name" "$options" approx_type=CLOUGH
#run_example "$example_name" "$options" approx_type=SECOND  # Broken?
