#!/bin/sh

#set -x

source "$LIBMESH_DIR"/examples/run_common.sh

example_name=miscellaneous_ex11

run_example "$example_name"

# No benchmarks here - at larger scales this problem spends 95%+ time
# in the linear solver
