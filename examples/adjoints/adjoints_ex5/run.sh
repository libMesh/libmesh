#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex5

example_dir=examples/adjoints/$example_name

# Uniform time stepping
# Save previous timestep solution in memory
run_example "$example_name" n_timesteps=10 timesolver_tolerance=0.0 timesolver_upper_tolerance=0.0 solution_history_type=memory
# Save previous timestep solution on disk
run_example "$example_name" n_timesteps=10 timesolver_tolerance=0.0 timesolver_upper_tolerance=0.0 solution_history_type=file

# Adaptive time stepping
# Save previous timestep solution in memory
run_example "$example_name" n_timesteps=7 timesolver_tolerance=1.0 timesolver_upper_tolerance=1.2 solution_history_type=memory
# Save previous timestep solution on disk
run_example "$example_name" n_timesteps=7 timesolver_tolerance=1.0 timesolver_upper_tolerance=1.2 solution_history_type=file
