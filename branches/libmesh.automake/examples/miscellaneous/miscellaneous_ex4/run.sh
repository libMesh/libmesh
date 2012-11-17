#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex4

message_running "$example_name" 
run_example "$example_name"
message_done_running "$example_name"
