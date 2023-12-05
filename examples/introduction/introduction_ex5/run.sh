#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=introduction_ex5

options="-q QGAUSS"

run_example "$example_name" "$options"
