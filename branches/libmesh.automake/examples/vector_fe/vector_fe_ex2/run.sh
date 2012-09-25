#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex2

link_if_needed $LIBMESH_DIR/$example_dir/vector_fe_ex2.in

message_running "$example_name" 

run_example "$example_name" 

message_done_running "$example_name"

discard_link $LIBMESH_DIR/$example_dir/vector_fe_ex2.in
