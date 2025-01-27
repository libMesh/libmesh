#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=vector_fe_ex10

# Note: these problems are particularly ill-conditioned, so we currently
# resort to a direct solver.

options="dim=2 element_type=TRI6 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 element_type=TRI7 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 element_type=QUAD8 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 element_type=QUAD9 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 order=3 element_type=TRI6 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 order=3 element_type=TRI7 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 order=3 element_type=QUAD8 -pc_type lu"
run_example "$example_name" "$options"

options="dim=2 order=3 element_type=QUAD9 -pc_type lu"
run_example "$example_name" "$options"

# Subdividing each hex into 24 tets gets expensive in dbg...
options="dim=3 element_type=TET14 grid_size=6 -pc_type lu"
run_example "$example_name" "$options"

options="dim=3 element_type=HEX27 -pc_type lu"
run_example "$example_name" "$options"
