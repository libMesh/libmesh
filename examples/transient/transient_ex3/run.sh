#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=transient_ex3
example_dir=examples/transient/$example_name

options="-i advection_2D.in"
run_example "$example_name" "$options"
