#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=transient_ex1
meshfile=examples/adaptivity/adaptivity_ex2/mesh.xda

link_if_needed $LIBMESH_DIR/$meshfile

message_running "$example_name" 

run_example "$example_name" 

message_done_running "$example_name"

discard_link $LIBMESH_DIR/$meshfile
