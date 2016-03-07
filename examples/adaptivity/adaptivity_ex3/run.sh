#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

run_example "$example_name" refinement_type=h
run_example "$example_name" refinement_type=p
run_example "$example_name" refinement_type=hp
run_example "$example_name" refinement_type=hpmatched
