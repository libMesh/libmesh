#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex3

# Note: these problems are particularly ill-conditioned, so we currently
# resort to a direct solver.

options="element_type=TRI6 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="element_type=TRI7 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="element_type=QUAD8 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="element_type=QUAD9 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="order=2 element_type=TRI6 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="order=2 element_type=TRI7 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="order=2 element_type=QUAD8 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"

options="order=2 element_type=QUAD9 -pc_type lu"
run_example_no_extra_options "$example_name" "$options"
