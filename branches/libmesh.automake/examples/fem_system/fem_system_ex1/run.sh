#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=fem_system_ex1

example_dir=examples/fem_system/$example_name

ln -sf $LIBMESH_DIR/$example_dir/navier.in .
ln -sf $LIBMESH_DIR/$example_dir/fem_system_ex1.in .

message_running "$example_name" 
run_example "$example_name"
message_done_running "$example_name"
