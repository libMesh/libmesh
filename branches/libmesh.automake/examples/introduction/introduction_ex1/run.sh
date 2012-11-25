#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex1

options="-d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda"
run_example "$example_name" "$options"

options="-d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda -o output.xda"
run_example "$example_name" "$options"
