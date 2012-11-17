#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex5

message_running "$example_name" 

options="-ksp_type cg"

run_example "$example_name" "$options"

message_done_running "$example_name"
