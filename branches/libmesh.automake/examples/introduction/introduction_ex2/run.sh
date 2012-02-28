#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex2

message_running "$example_name" 
run_example "$example_name"
echo " "
options="eqn_sys.dat"
run_example "$example_name" "$options"
message_done_running "$example_name"
