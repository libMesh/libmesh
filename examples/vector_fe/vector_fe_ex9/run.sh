#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex9

options="cavity=true mms=true element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=true element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"

options="cavity=true mms=false element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=false element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"
