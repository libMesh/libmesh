#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adjoints_ex1

example_dir=examples/adjoints/$example_name

run_example "$example_name"

# Exercise our different partitioner types, with smaller meshes for
# speed.
run_example "$example_name" mesh_partitioner_type=Metis coarserefinements=2 max_adaptivesteps=4
run_example "$example_name" mesh_partitioner_type=Parmetis coarserefinements=2 max_adaptivesteps=4
run_example "$example_name" mesh_partitioner_type=SFCurves coarserefinements=2 max_adaptivesteps=4
run_example "$example_name" mesh_partitioner_type=Hilbert coarserefinements=2 max_adaptivesteps=4
run_example "$example_name" mesh_partitioner_type=Morton coarserefinements=2 max_adaptivesteps=4
run_example "$example_name" mesh_partitioner_type=Linear coarserefinements=2 max_adaptivesteps=4
run_example "$example_name" mesh_partitioner_type=Centroid coarserefinements=2 max_adaptivesteps=4
