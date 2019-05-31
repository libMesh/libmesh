#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=fem_system_ex1

example_dir=examples/fem_system/$example_name

# First, run with standard options from input files.
run_example "$example_name"

# Next, lets run with pure fieldsplit without gmg
fs_lu_options="--use_petsc_dm --node-major-dofs \
-snes_view -snes_monitor -snes_converged_reason -snes_rtol 1.0e-4 \
-ksp_type fgmres -ksp_converged_reason -ksp_monitor_true_residual \
-ksp_rtol 1.0e-5 \
-pc_type fieldsplit -pc_fieldsplit_0_fields 0,1 -pc_fieldsplit_1_fields 2 \
-pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition a11 \
-fieldsplit_0_pc_type lu \
-fieldsplit_0_ksp_type preonly \
-fieldsplit_0_ksp_converged_reason \
-fieldsplit_p_pc_type jacobi \
-fieldsplit_p_ksp_type gmres \
-fieldsplit_p_ksp_converged_reason \
-fieldsplit_p_inner_pc_type lu -fieldsplit_p_inner_ksp_type preonly \
-fieldsplit_p_upper_pc_type lu -fieldsplit_p_upper_ksp_type preonly"

run_example "$example_name" "solver_type=petscdiff coarsegridsize=6 transient=false n_timesteps=1 $fs_lu_options"

# And again with some timestepping coverage
run_example "$example_name" "solver_type=petscdiff coarsegridsize=6 coarserefinements=1 transient=true n_timesteps=3 $fs_lu_options"

# Rerun with options demonstrating fieldpslit+gmg on the velocity block.
fs_gmg_options="--use_petsc_dm --node-major-dofs \
 -snes_view -snes_monitor -snes_converged_reason -snes_rtol 1.0e-4 \
 -ksp_type fgmres -ksp_rtol 1.0e-5 \
-pc_type fieldsplit -pc_fieldsplit_0_fields 0,1 -pc_fieldsplit_1_fields 2 \
-pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition a11 \
-fieldsplit_0_pc_type mg -fieldsplit_0_pc_mg_galerkin both -fieldsplit_0_pc_mg_type full \
-pc_mg_levels 3 -fieldsplit_0_pc_mg_levels 3 \
-fieldsplit_0_mg_levels_pc_type sor -fieldsplit_0_mg_levels_pc_sor_its 5 \
-fieldsplit_0_mg_levels_ksp_type richardson -fieldsplit_0_mg_levels_ksp_richardson_self_scale \
-fieldsplit_0_mg_levels_ksp_max_it 5 \
-fieldsplit_0_mg_coarse_pc_type lu -fieldsplit_0_mg_coarse_ksp_type preonly \
-fieldsplit_0_ksp_max_it 1 -fieldsplit_0_ksp_type gmres \
-fieldsplit_p_ksp_max_it 5 -fieldsplit_p_ksp_type gmres \
-fieldsplit_p_pc_type none -fieldsplit_p_ksp_rtol 1.0e-4 \
-fieldsplit_p_inner_pc_type lu -fieldsplit_p_inner_ksp_type preonly \
-fieldsplit_p_upper_pc_type lu -fieldsplit_p_upper_ksp_type preonly"

run_example "$example_name" "solver_type=petscdiff coarsegridsize=6 coarserefinements=2 transient=false n_timesteps=1 $fs_gmg_options"

# Now rerun same thing with over a distributed mesh
run_example "$example_name" "mesh_type=distributed solver_type=petscdiff coarsegridsize=6 coarserefinements=2 transient=false n_timesteps=1 $fs_gmg_options"
