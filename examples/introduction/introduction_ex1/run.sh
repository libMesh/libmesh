#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex1

options="-d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda"
run_example "$example_name" "$options"

options="-d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda -o output.xda"
run_example "$example_name" "$options"

# XDR format
options="-d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda -o output.xdr"
run_example "$example_name" "$options"

# Read XDA 0_9_2 format
options="-d 3 output.xda"
run_example "$example_name" "$options"

# Read XDR 0_9_2 format
options="-d 3 output.xdr"
run_example "$example_name" "$options"
