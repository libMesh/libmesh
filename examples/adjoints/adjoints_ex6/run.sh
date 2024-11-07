#!/bin/sh

#set -x

. "$LIBMESH_DIR"/examples/run_common.sh

example_name=adjoints_ex6

example_dir=examples/adjoints/$example_name

run_example "$example_name"

benchmark_example 1 "$example_name" max_adaptivesteps=5
