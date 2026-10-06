#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=fem_system_ex6

run_example "$example_name"

# Check the dual-number Jacobian against finite differences.
run_example "$example_name" verify_analytic_jacobians=1.e-6 grid_size=4
