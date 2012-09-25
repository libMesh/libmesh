#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex5

example_dir=examples/miscellaneous/miscellaneous_ex5
adj_example_dir=examples/adjoints/adjoints_ex1

link_if_needed $LIBMESH_DIR/$example_dir/miscellaneous_ex5.in
link_if_needed $LIBMESH_DIR/$example_dir/lshaped3D.xda
link_if_needed $LIBMESH_DIR/$adj_example_dir/lshaped.xda

ln -sf $LIBMESH_DIR/$example_dir/lshaped3D.xda .
ln -sf $LIBMESH_DIR/$example_dir/miscellaneous_ex5.in .

message_running "$example_name" 
run_example "$example_name" 
message_done_running "$example_name" 

discard_link $LIBMESH_DIR/$example_dir/miscellaneous_ex5.in
discard_link $LIBMESH_DIR/$example_dir/lshaped3D.xda
discard_link $LIBMESH_DIR/$adj_example_dir/lshaped.xda
