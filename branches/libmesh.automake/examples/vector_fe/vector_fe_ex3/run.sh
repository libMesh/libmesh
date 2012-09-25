#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex3

link_if_needed $LIBMESH_DIR/$example_dir/vector_fe_ex3.in

options="element_type=TRI6 -pc_type jacobi"
message_running "$example_name" "$options"
run_example "$example_name" "$options"
message_done_running "$example_name" "$options"

options="element_type=QUAD8 -pc_type jacobi"
message_running "$example_name" "$options"
run_example "$example_name" "$options"
message_done_running "$example_name" "$options"

options="element_type=QUAD9 -pc_type jacobi"
message_running "$example_name" "$options"
run_example "$example_name" "$options"
message_done_running "$example_name" "$options"

discard_link $LIBMESH_DIR/$example_dir/vector_fe_ex3.in
