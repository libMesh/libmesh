#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=optimization_ex1

example_dir=examples/optimization/$example_name

options="-tao_monitor -tao_view -tao_type nls"
run_example "$example_name" "$options"
