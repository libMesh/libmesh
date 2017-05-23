#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex6

example_dir=examples/adjoints/$example_name

run_example "$example_name"
