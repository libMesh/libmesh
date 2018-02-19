#!/bin/bash

# This script assumes you have libMesh installed at $LIBMESH_DIR, and
# calls test_installed_headers.sh with specific arguments.
script_full_path=$(dirname "$0")
test_CXXFLAGS="`$LIBMESH_DIR/bin/libmesh-config --cppflags --cxxflags --include`" HEADERS_TO_TEST="`find $LIBMESH_DIR/include/libmesh -name "*.h" -type f -exec basename {} \;`" ${script_full_path}/test_installed_headers.sh
