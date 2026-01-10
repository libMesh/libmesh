#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=vector_fe_ex4

# Note: these problems are particularly ill-conditioned, so we use a robust
# preconditioner, hypre AMS. We use jacobi in two configs for coverage only.

options="element_type=TET10 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"

options="element_type=TET14 -pc_type hypre -pc_hypre_type ams"
run_example_no_extra_options "$example_name" "$options"

options="element_type=HEX20 -pc_type jacobi"
run_example_no_extra_options "$example_name" "$options"

options="element_type=HEX27 -pc_type hypre -pc_hypre_type ams"
run_example_no_extra_options "$example_name" "$options"
