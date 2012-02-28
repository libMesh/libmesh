#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

ln -sf $LIBMESH_DIR/$example_dir/lshaped.xda .
ln -sf $LIBMESH_DIR/$example_dir/lshaped3D.xda .
ln -sf $LIBMESH_DIR/$example_dir/$example_name.in

message_running "$example_name" 
run_example "$example_name"
message_done_running "$example_name"
