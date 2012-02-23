#!/bin/sh

set -x

example_name=adaptivity_ex5
example_dir=examples/adaptivity/$example_name

# Specify the number of timesteps to do
n_timesteps=25

# Specify the amount of initial uniform refinement to do
n_refinements=5

# Specify how often to write output files
output_freq=10

echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name \[-read_solution\] -n_timesteps $n_timesteps -n_refinements $n_refinements -init_timestep \[0\|$n_timesteps\] $LIBMESH_OPTIONS
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name -n_timesteps $n_timesteps -n_refinements $n_refinements \
    -output_freq $output_freq -init_timestep 0 $LIBMESH_OPTIONS || exit 1
echo " "
echo "***** Finished first" $n_timesteps "steps, now read in" \
    "saved solution and continue *****"
echo " "
$LIBMESH_RUN ./$example_name -read_solution -n_timesteps $n_timesteps \
    -output_freq $output_freq -init_timestep $n_timesteps $LIBMESH_OPTIONS || exit 1
echo " "
$LIBMESH_RUN ./$example_name -n_timesteps $n_timesteps -n_refinements $n_refinements \
    -output_freq $output_freq -init_timestep 0 \
    -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)' \
    $LIBMESH_OPTIONS || exit 1
echo " "
echo "***** Finished first" $n_timesteps "steps, now read in" \
    "saved solution and continue *****"
echo " "
$LIBMESH_RUN ./$example_name -read_solution -n_timesteps $n_timesteps \
    -output_freq $output_freq -init_timestep $n_timesteps \
    $LIBMESH_OPTIONS || exit 1
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name \[-read_solution\] -n_timesteps $n_timesteps -init_timestep \[0\|$n_timesteps\] $LIBMESH_OPTIONS
echo "***************************************************************"
