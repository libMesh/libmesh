#!/bin/sh

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=reduced_basis_ex5

example_dir=examples/reduced_basis/$example_name

options="-online_mode 0 -ksp_type cg"
run_example_no_extra_options "$example_name" "$options"

options="-online_mode 1"
run_example_no_extra_options "$example_name" "$options"

options="-online_mode 0 -ksp_type cg mesh_filename=Twisted_Beam.bez.gz order=2"
run_example_no_extra_options "$example_name" "$options"

options="-online_mode 1 mesh_filename=Twisted_Beam.bez.gz order=2"
run_example_no_extra_options "$example_name" "$options"
