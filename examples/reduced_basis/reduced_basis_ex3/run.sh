#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=reduced_basis_ex3

example_dir=examples/reduced_basis/$example_name

# PETSc default Block-Jacobi+ILU tends to fail to converge for us on
# OSX on this problem.
petsc_options="-pc_type jacobi"

options="-online_mode 0"
run_example "$example_name" "$options" "$petsc_options"
benchmark_example 1 "$example_name" "$options n_elem=150"

options="-online_mode 1"
run_example "$example_name" "$options" "$petsc_options"
benchmark_example 1 "$example_name" "$options n_elem=150"
