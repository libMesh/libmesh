#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex3

example_dir=examples/adjoints/$example_name

run_example "$example_name" $options

# This example needs better solvers on larger meshes + processor counts
#benchmark_example 1 "$example_name" "max_adaptivesteps=18 -pc_type bjacobi -sub_pc_factor_levels 4"

# This example needs even better solvers on larger meshes + much larger processor counts

# This is the best I can do that works on any processor count, but on
# one core it takes 6 minutes...
#benchmark_example 1 "$example_name" "max_adaptivesteps=15 --solver_variable_names --solver_group_u=0 --solver_group_v=0 --solver_group_C=0 --solver_group_p=1 -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_pc_type lu -fieldsplit_1_pc_type jacobi -fieldsplit_1_ksp_type tfqmr -ksp_monitor -fieldsplit_1_ksp_converged_reason"

# This only takes me 17s on one proc, 10s on two ... but 20s on three,
# then hits 10k iterations and fails line search on four
#benchmark_example 1 "$example_name" "max_adaptivesteps=15 --solver_variable_names --solver_group_u=0 --solver_group_v=0 --solver_group_C=0 --solver_group_p=1 -pc_type bjacobi -sub_pc_type ilu"

# Even more factor levels plus a shift type and we don't fail until
# 13 ... but this is starting to make us as slow as Schur
#benchmark_example 1 "$example_name" "max_adaptivesteps=15 --solver_variable_names --solver_group_u=0 --solver_group_v=0 --solver_group_C=0 --solver_group_p=1 -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_shift_type nonzero"
