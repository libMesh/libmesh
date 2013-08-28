#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex4

# Note: these problems are particularly ill-conditioned, and common
# preconditioners like ILU don't all work very well.  We currently
# force Jacobi in these examples.

options="element_type=HEX20 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"

options="element_type=HEX27 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"
