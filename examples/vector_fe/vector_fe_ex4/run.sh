#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex4

# Note: these problems are particularly ill-conditioned, and standard
# preconditioners like jacobi and ILU don't seem to work very well.
# We've therefore instead used LU in these examples.

options="element_type=HEX20 -pc_type lu"
run_example "$example_name" "$options"

options="element_type=HEX27 -pc_type lu"
run_example "$example_name" "$options"
