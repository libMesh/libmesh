#!/bin/sh

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=systems_of_equations_ex7

options="-ksp_type cg -pc_type bjacobi -snes_linesearch_type basic"

run_example "$example_name" "$options"

benchmark_example 1 "$example_name" "$options n_elem_x=75 n_elem_y=15 n_elem_z=15"
