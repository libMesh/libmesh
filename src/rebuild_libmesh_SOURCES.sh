#!/bin/bash

sources=`find base error_estimation fe geom mesh numerics parallel partitioning physics quadrature reduced_basis solution_transfer solvers systems utils -name "*.C" -o -name "*.c" -o -name "*.data" -type f | LC_COLLATE=POSIX sort`

echo "# Do not edit - automatically generated from $0" > libmesh_SOURCES
printf '%s' "libmesh_SOURCES = " >> libmesh_SOURCES
for source_with_path in $sources ; do

    echo " \\" >> libmesh_SOURCES
    printf '%s' "        "src/$source_with_path >> libmesh_SOURCES

done


echo " " >> libmesh_SOURCES

