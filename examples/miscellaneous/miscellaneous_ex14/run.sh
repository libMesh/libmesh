#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=miscellaneous_ex14

example_dir=examples/miscellaneous/miscellaneous_ex14

run_example "$example_name" $srcdir/sample.in

