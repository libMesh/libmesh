#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex6

run_example "$example_name"

# No benchmarking here since the work is done upstream of libMesh
