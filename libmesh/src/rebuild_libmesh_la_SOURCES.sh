#!/bin/bash

sources=`find base error_estimation fe geom mesh numerics parallel partitioning physics quadrature reduced_basis solvers systems utils -name "*.C" -o -name "*.c" -type f | sort`

echo "# Do not edit - automatically generated from $0" > libmesh_la_SOURCES
echo -n "libmesh_la_SOURCES = " >> libmesh_la_SOURCES
for source_with_path in $sources ; do
    
    echo " \\" >> libmesh_la_SOURCES
    echo -n "        "src/$source_with_path >> libmesh_la_SOURCES

done


echo " " >> libmesh_la_SOURCES

