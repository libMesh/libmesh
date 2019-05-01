#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

run_example "$example_name" refinement_type=h
# Some solvers still give us trouble with too much p
run_example "$example_name" refinement_type=p max_r_steps=4
run_example "$example_name" refinement_type=hp max_r_steps=6
run_example "$example_name" refinement_type=matchedhp max_r_steps=4

# Let's get some 1D coverage too; that's cheap.
run_example "$example_name" dimension=1 refinement_type=h
run_example "$example_name" dimension=1 refinement_type=p
run_example "$example_name" dimension=1 refinement_type=hp max_r_steps=6
run_example "$example_name" dimension=1 refinement_type=matchedhp max_r_steps=4

# And try out another element type
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=h
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=p
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=hp max_r_steps=8
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=3 refinement_type=matchedhp max_r_steps=4
