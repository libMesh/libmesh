#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex5
example_dir=examples/adaptivity/$example_name

# Specify the number of timesteps to do
n_timesteps=25

# Specify the amount of initial uniform refinement to do
n_refinements=5

# Specify how often to write output files
output_freq=10

options="-n_timesteps $n_timesteps -n_refinements $n_refinements \
    -output_freq $output_freq -init_timestep 0 \
    -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2)+pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)'"


# This is frustrating - I cannot get the fparser expression to survive as a shell variable

# run_example "$example_name" $options
# run_example "$example_name" -n_timesteps $n_timesteps -n_refinements $n_refinements \
#         -output_freq $output_freq -init_timestep 0 \
#         -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)'

executable=$example_name
if (test ! -x ${executable}); then
    echo "ERROR: cannot find ${executable}!"
    exit 1
fi

message_running $example_name $executable $options
    
$LIBMESH_RUN ./$executable -n_timesteps $n_timesteps -n_refinements $n_refinements \
    -output_freq $output_freq -init_timestep 0 \
    -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)' \
    $LIBMESH_OPTIONS || exit 1
    
message_done_running $example_name $executable $options


echo " "
echo "***** Finished first" $n_timesteps "steps, now read in" \
    "saved solution and continue *****"
echo " "

options="-read_solution -n_timesteps $n_timesteps \
    -output_freq $output_freq -init_timestep $n_timesteps"
run_example "$example_name" "$options"

options="\[-read_solution\] -n_timesteps $n_timesteps -init_timestep \[0\|$n_timesteps\]"
