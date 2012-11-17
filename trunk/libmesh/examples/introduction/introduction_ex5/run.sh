#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex5

options="-q 0"

message_running "$example_name" "$options" 
run_example "$example_name" "$options"
message_done_running "$example_name" "$options"
