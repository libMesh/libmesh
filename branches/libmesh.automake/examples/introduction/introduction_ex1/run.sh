#!/bin/sh

#set -x

example_name=introduction_ex1

echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name -d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda $LIBMESH_OPTIONS || exit 1
$LIBMESH_RUN ./$example_name -d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda -o output.xda $LIBMESH_OPTIONS || exit 1
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
