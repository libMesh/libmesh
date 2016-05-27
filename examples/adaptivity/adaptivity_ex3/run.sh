#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

run_example "$example_name" refinement_type=h
run_example "$example_name" refinement_type=p

# Some solvers still give us trouble with too much hp
run_example "$example_name" refinement_type=hp max_r_steps=8
run_example "$example_name" refinement_type=matchedhp max_r_steps=4
