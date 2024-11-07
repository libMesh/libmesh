#!/bin/bash

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=fem_system_ex3

example_dir=examples/fem_system/$example_name

common_options="write_interval=1 solver_quiet=false relative_step_tolerance=1e-3 relative_residual_tolerance=1.e-3"

# Note: Too much ILU fails badly on this problem in single precision.
# We force simple Jacobi to be safe.
petsc_options="-ksp_type cg -pc_type jacobi"

# Note: Use 25 timesteps to simulate approximately three periods of oscillation.
options="deltat=0.25 n_timesteps=5 time_solver=newmark $common_options $petsc_options"
run_example_no_extra_options "$example_name" "$options"
benchmark_example 1 "$example_name" "$options nx=64 ny=16 nz=8"

options="time_solver=steady n_timesteps=1 $common_options $petsc_options"
run_example_no_extra_options "$example_name" "$options"
benchmark_example 1 "$example_name" "$options nx=48 ny=12 nz=6"

# With first order solvers, the Jacobian is no longer symmetric
petsc_options="-ksp_type gmres -pc_type bjacobi -sub_pc_type ilu"

options="deltat=0.25 n_timesteps=5 time_solver=euler theta=0.5 $common_options $petsc_options"
run_example_no_extra_options "$example_name" "$options"
benchmark_example 1 "$example_name" "$options nx=48 ny=12 nz=6"

options="deltat=0.25 n_timesteps=5 time_solver=euler2 theta=0.5 $common_options $petsc_options"
run_example_no_extra_options "$example_name" "$options"
benchmark_example 1 "$example_name" "$options nx=48 ny=12 nz=6"
