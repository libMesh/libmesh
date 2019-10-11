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
#
# for PETSc versions > 3.9.5, the selection of solver has changed, so we loop over
# both options and skip the runs with incompatible arguments.

options="-online_mode 0 -ksp_type preonly -pc_type lu"
for solver in "-pc_factor_mat_solver_type mumps" "-pc_factor_mat_solver_package mumps" #"-pc_factor_mat_solver_type superlu" "-pc_factor_mat_solver_package superlu"
do
   run_example_no_extra_options "$example_name" "$options $solver"

   options="-online_mode 1"
   run_example "$example_name" "$options $solver"
done
