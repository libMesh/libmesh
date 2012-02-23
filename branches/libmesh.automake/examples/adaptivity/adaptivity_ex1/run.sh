#!/bin/sh

#set -x

example_name=adaptivity_ex1

echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name $LIBMESH_OPTIONS
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
