#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex2

message_running "$example_name" 

options="-n 5 -eps_type lapack"
run_example "$example_name" "$options"

message_done_running "$example_name"
