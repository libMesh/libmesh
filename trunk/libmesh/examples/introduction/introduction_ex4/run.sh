#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex4

message_running "$example_name"

options="-d 1 -n 20"
run_example "$example_name" "$options"

options="-d 2 -n 15"
run_example "$example_name" "$options"

options="-d 3 -n  6"
run_example "$example_name" "$options"

message_done_running "$example_name"
