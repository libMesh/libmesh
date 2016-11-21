#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex4

run_example "$example_name" approx_type=HERMITE

# CLOUGH elements don't currently support threads
run_example_no_extra_options "$example_name" approx_type=CLOUGH $LIBMESH_OPTIONS --n_threads=1
#run_example "$example_name" approx_type=SECOND  # Broken?
