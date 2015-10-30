#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=fem_system_ex4

example_dir=examples/fem_system/$example_name

run_example "$example_name"
