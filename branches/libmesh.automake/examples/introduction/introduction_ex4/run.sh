#!/bin/sh

#set -x

example_name=introduction_ex4

echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name -d 1 -n 20 $LIBMESH_OPTIONS || exit 1
$LIBMESH_RUN ./$example_name -d 2 -n 15 $LIBMESH_OPTIONS || exit 1
$LIBMESH_RUN ./$example_name -d 3 -n  6 $LIBMESH_OPTIONS || exit 1
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
