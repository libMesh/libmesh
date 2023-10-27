#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex7

options="dim=2 element_type=TRI6"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=TRI7"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=QUAD8"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=QUAD9"
run_example_no_extra_options "$example_name" "$options"

# Subdividing each hex into 24 tets gets expensive in dbg...
options="dim=3 element_type=TET14 grid_size=4"
run_example_no_extra_options "$example_name" "$options"

options="dim=3 element_type=HEX27 grid_size=4"
run_example_no_extra_options "$example_name" "$options"
