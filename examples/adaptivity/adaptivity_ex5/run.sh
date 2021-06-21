#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex5

# Specify the number of timesteps to do
n_timesteps=25

# Specify the amount of initial uniform refinement to do
n_refinements=5

# Specify how often to write output files
output_freq=10

# Lack of '' around the fparser function is intentional, as is lack of
# whitespace in the function.  Yes, we need a better way to pass
# these.
options="-n_timesteps $n_timesteps -n_refinements $n_refinements \
    -output_freq $output_freq -init_timestep 0 \
    -exact_solution 10*exp(-(pow(x-0.8*t-0.2,2)+pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)"

run_example "$example_name" $options

echo " "
echo "***** Finished first" $n_timesteps "steps, now read in" \
    "saved solution and continue *****"
echo " "

options="-read_solution -n_timesteps $n_timesteps \
    -output_freq $output_freq -init_timestep $n_timesteps"
run_example "$example_name" "$options"
