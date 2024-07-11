#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex9

options="cavity=true mms=true element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=true element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"

# options="cavity=false mms=true dim=2 element_type=TRI7 -ksp_type preonly -condensed_pc_type \
# fieldsplit -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic \
# -condensed_pc_fieldsplit_type schur -condensed_pc_fieldsplit_schur_fact_type full \
# -condensed_pc_fieldsplit_schur_precondition selfp -condensed_pc_fieldsplit_0_fields 4,5 \
# -condensed_pc_fieldsplit_1_fields 6 --use_petsc_dm --node-major-dofs"
# run_example_no_extra_options "$example_name" "$options"

options="cavity=true mms=false element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"

options="cavity=false mms=false element_type=TRI6 -ksp_type preonly -condensed_pc_type lu -condensed_pc_factor_shift_type NONZERO -snes_linesearch_type basic"
run_example_no_extra_options "$example_name" "$options"
