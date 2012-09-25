#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=reduced_basis_ex2

example_dir=examples/reduced_basis/$example_name

link_if_needed $LIBMESH_DIR/$example_dir/reduced_basis_ex2.in

message_running "$example_name" 

options="-online_mode 0 -eps_type lapack"
run_example "$example_name" "$options"

options="-online_mode 1"
run_example "$example_name" "$options"

message_done_running "$example_name"

discard_link $LIBMESH_DIR/$example_dir/reduced_basis_ex2.in
