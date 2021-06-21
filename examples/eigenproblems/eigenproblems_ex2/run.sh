#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex2

options="-n 5 -eps_type lapack"
run_example "$example_name" "$options"

# No benchmark_example here - it spends 99.8+% of its time in the
# SLEPc solve at even small scale, with that eps_type
