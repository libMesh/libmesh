#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=vector_fe_ex6

### Neumann boundary conditions

options="dim=2 element_type=TRI6 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=TRI7 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=QUAD8 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=QUAD9 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=TRI6 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=TRI7 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=QUAD8 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=QUAD9 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

# Subdividing each hex into 24 tets gets expensive in dbg...
options="dim=3 element_type=TET14 boundary_condition=neumann grid_size=6"
run_example_no_extra_options "$example_name" "$options"

options="dim=3 element_type=HEX27 boundary_condition=neumann"
run_example_no_extra_options "$example_name" "$options"

### Dirichlet boundary conditions

options="dim=2 element_type=TRI6 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=TRI7 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=QUAD8 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 element_type=QUAD9 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=TRI6 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=TRI7 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=QUAD8 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

options="dim=2 order=2 element_type=QUAD9 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"

# Subdividing each hex into 24 tets gets expensive in dbg...
options="dim=3 element_type=TET14 boundary_condition=dirichlet grid_size=6"
run_example_no_extra_options "$example_name" "$options"

options="dim=3 element_type=HEX27 boundary_condition=dirichlet"
run_example_no_extra_options "$example_name" "$options"
