#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex3

example_dir=examples/adjoints/$example_name

run_example "$example_name" $options

# This example needs better solvers on larger meshes + processor counts
benchmark_example 1 "$example_name" "max_adaptivesteps=18 -pc_type bjacobi -sub_pc_factor_levels 4"
