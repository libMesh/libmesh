#!/bin/sh

#set -x

example_name=introduction_ex1
example_path=introduction/$example_name

PATH=.:$PATH

echo "example_name=$example_name" 
cd $example_path || exit 1


echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name -d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda $LIBMESH_OPTIONS
$LIBMESH_RUN ./$example_name -d 3 $LIBMESH_DIR/reference_elements/3D/one_hex27.xda -o output.xda $LIBMESH_OPTIONS
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name
echo "***************************************************************"
