#!/bin/sh

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=miscellaneous_ex8

options=""

run_example "$example_name" "$options"

# This is useless if we might want to run parallel benchmarks:
# RadialBasisInterpolation::prepare_for_use() has a
# very-embarrassingly-serial Eigen solve
#benchmark_example 1 "$example_name" "n_refinements=1 n_source_points=10000 n_target_points=100"
