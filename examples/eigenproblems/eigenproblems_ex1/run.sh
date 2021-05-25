#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex1

options="-n 5"
run_example "$example_name" "$options"

# No benchmark_example here - it spends 98+% of its time in the SLEPc
# solve at scale, no room to improve on our end.
