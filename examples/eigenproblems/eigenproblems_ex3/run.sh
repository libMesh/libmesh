#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=eigenproblems_ex3

example_dir=examples/eigenproblems/$example_name

slepc_options="-st_ksp_type gmres -st_pc_type bjacobi -st_sub_pc_type ilu -st_type sinvert"
options="-n_evals 5 -mesh_name drum1 -plotting_index 2 $slepc_options"
run_example_no_extra_options "$example_name" "$options"

options="-n_evals 5 -mesh_name drum2 -plotting_index 2 $slepc_options"
run_example_no_extra_options "$example_name" "$options"

# No benchmark problems here - these are run on mesh files without
# scaling built into the examples.
