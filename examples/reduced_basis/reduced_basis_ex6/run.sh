 #!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=reduced_basis_ex6

example_dir=examples/reduced_basis/$example_name

options="-online_mode 0 -mat_new_nonzero_allocation_err false"
run_example "$example_name" "$options"

options="-online_mode 1 -mat_new_nonzero_allocation_err false"
run_example "$example_name" "$options"
