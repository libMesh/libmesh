#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex3

# Note: these problems are particularly ill-conditioned, and common
# preconditioners like ILU don't all work very well.  We currently
# force Jacobi in these examples.

options="element_type=TRI6 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"

options="element_type=QUAD8 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"

options="element_type=QUAD9 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"

# No benchmark until we decide what to do about no_extra_options there
