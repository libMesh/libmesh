#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex8

options="dim=2 element_type=TRI6"
run_example_no_extra_options "$example_name" "$options"
run_example_no_extra_options "$example_name" "$options family=L2_HIERARCHIC vecfamily=L2_HIERARCHIC_VEC"
run_example_no_extra_options "$example_name" "$options family=L2_HIERARCHIC vecfamily=L2_HIERARCHIC_VEC order=2"
run_example_no_extra_options "$example_name" "$options family=L2_HIERARCHIC vecfamily=L2_HIERARCHIC_VEC order=3"
run_example_no_extra_options "$example_name" "$options family=XYZ"

options="dim=2 element_type=TRI7"
run_example_no_extra_options "$example_name" "$options"
run_example_no_extra_options "$example_name" "$options family=L2_HIERARCHIC vecfamily=L2_HIERARCHIC_VEC"

# Subdividing each hex into 24 tets gets expensive in dbg...
options="dim=3 element_type=TET14 grid_size=4"
run_example_no_extra_options "$example_name" "$options"
run_example_no_extra_options "$example_name" "$options family=L2_HIERARCHIC vecfamily=L2_HIERARCHIC_VEC"
run_example_no_extra_options "$example_name" "$options family=XYZ"
