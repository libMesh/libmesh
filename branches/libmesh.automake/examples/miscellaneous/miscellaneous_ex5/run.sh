#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex5

adj_example_dir=examples/adjoints/adjoints_ex1

ln -sf $LIBMESH_DIR/$adj_example_dir/lshaped.xda .

example_dir=examples/miscellaneous/miscellaneous_ex5

ln -sf $LIBMESH_DIR/$example_dir/lshaped3D.xda .
ln -sf $LIBMESH_DIR/$example_dir/miscellaneous_ex5.in .

message_running "$example_name" 
run_example "$example_name" 
message_done_running "$example_name" 
