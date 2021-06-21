#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex7

options="--verbose dim=1 N=1024 initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-8"

run_example "$example_name" "$options"

# No benchmark here - I can't seem to pick parameters that aren't
# either way too easy or prone to solver convergence failure
