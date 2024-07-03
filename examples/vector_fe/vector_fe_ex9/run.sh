#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex9

options="cavity=true mms=true element_type=TRI6 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=true mms=true dim=2 element_type=TRI7 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=true element_type=TRI6 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=true dim=2 element_type=TRI7 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=true mms=false element_type=TRI6 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=true mms=false dim=2 element_type=TRI7 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=false element_type=TRI6 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=false dim=2 element_type=TRI7 grid_size=1 -ksp_type preonly -condensed_pc_factor_shift_type NONZERO"
run_example_no_extra_options "$example_name" "$options"
