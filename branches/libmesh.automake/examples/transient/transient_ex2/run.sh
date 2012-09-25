#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=transient_ex2
example_dir=examples/transient/$example_name

link_if_needed $LIBMESH_DIR/$example_dir/pipe-mesh.unv

message_running "$example_name" 

options="pipe-mesh.unv"
run_example "$example_name" "$options"

message_done_running "$example_name"

discard_link $LIBMESH_DIR/$example_dir/pipe-mesh.unv
