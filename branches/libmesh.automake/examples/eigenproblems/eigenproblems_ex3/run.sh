#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex3

example_dir=examples/eigenproblems/$example_name

message_running "$example_name" 

options="-n_evals 5 -mesh_name drum1 -plotting_index 2"
run_example "$example_name" "$options"

options="-n_evals 5 -mesh_name drum2 -plotting_index 2"
run_example "$example_name" "$options"

message_done_running "$example_name"
