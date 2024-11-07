#!/bin/sh

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=adaptivity_ex4

run_example "$example_name" approx_type=HERMITE
run_example "$example_name" approx_type=HERMITE penalty_dirichlet=false

# HERMITE elements can p refine in 1D
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=FOURTH
run_example "$example_name" dimension=1 approx_type=HERMITE approx_order=FOURTH penalty_dirichlet=false

# CLOUGH elements don't currently support threads
run_example_no_extra_options "$example_name" approx_type=CLOUGH $LIBMESH_OPTIONS --n_threads=1
run_example_no_extra_options "$example_name" approx_type=CLOUGH $LIBMESH_OPTIONS --n_threads=1 penalty_dirichlet=false
#run_example "$example_name" approx_type=SECOND  # Broken?

# Examples to use for benchmarking
benchmark_example 1 "$example_name" dimension=3 approx_type=HERMITE approx_order=THIRD max_r_steps=6
