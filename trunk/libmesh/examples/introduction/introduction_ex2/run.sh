#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex2

run_example "$example_name"
echo " "
options="eqn_sys.dat"
run_example "$example_name" "$options"
