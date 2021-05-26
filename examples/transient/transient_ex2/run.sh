#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=transient_ex2
example_dir=examples/transient/$example_name

options="pipe-mesh.unv"
run_example "$example_name" "$options"

# No benchmark here; "result_node = 274" makes mesh modification iffy
