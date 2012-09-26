#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

message_running "$example_name" 
run_example "$example_name"
message_done_running "$example_name"
