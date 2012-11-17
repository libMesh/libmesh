#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=eigenproblems_ex1

message_running "$example_name"

options="-n 5"
run_example "$example_name" "$options"

message_done_running "$example_name"
