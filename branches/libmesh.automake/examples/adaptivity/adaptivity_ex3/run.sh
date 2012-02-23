#!/bin/sh

set -x

example_name=adaptivity_ex3
example_dir=examples/adaptivity/$example_name

ln -sf $LIBMESH_DIR/$example_dir/lshaped.xda .
ln -sf $LIBMESH_DIR/$example_dir/lshaped3D.xda .
ln -sf $LIBMESH_DIR/$example_dir/$example_name.in

echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name $LIBMESH_OPTIONS
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
