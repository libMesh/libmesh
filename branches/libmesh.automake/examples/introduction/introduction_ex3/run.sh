#!/bin/sh

#set -x

example_name=introduction_ex3

echo "***************************************************************"
echo "* Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
echo " "
$LIBMESH_RUN ./$example_name $LIBMESH_OPTIONS || exit 1
echo " "
echo "***************************************************************"
echo "* Done Running Example " $LIBMESH_RUN $example_name $LIBMESH_OPTIONS
echo "***************************************************************"
