#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=reduced_basis_ex7

example_dir=examples/reduced_basis/$example_name

# This example requires a direct solver for -online_mode=0, e.g.
#
# -ksp_type preonly -pc_type lu
#
# or, using an external direct solver like mumps:
#
# -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps
options="-online_mode 0"
run_example "$example_name" "$options"

options="-online_mode 1"
run_example "$example_name" "$options"
