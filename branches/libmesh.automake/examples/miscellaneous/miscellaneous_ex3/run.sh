#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex3

options="-r 3 -o FIRST"

message_running "$example_name" "$options"
run_example "$example_name" "$options"
message_done_running "$example_name" "$options"
