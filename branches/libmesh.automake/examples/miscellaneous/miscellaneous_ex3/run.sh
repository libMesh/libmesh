#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex3

link_if_needed $LIBMESH_DIR/examples/adjoints/adjoints_ex1/lshaped.xda

options="-r 3 -o FIRST"

message_running "$example_name" "$options"
run_example "$example_name" "$options"
message_done_running "$example_name" "$options"

discard_link $LIBMESH_DIR/examples/adjoints/adjoints_ex1/lshaped.xda
