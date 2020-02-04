#!/bin/bash

# these specific headers are required to build libMesh
# but we do not want to install them!
noinst_blacklist="mesh/nemesis_io_helper.h \
numerics/laspack_matrix.h \
numerics/laspack_vector.h \
solvers/laspack_linear_solver.h \
parallel/parallel_conversion_utils.h \
parallel/parallel_hilbert.h \
partitioning/parmetis_helper.h"

echo "# Do not edit - automatically generated from $0" > include_HEADERS


echo "# These are headers we need to compile the libMesh library but should" >> include_HEADERS
echo "# not be installed.  These header files may have implementation" >> include_HEADERS
echo "# details which are not required for the public interface" >> include_HEADERS

printf '%s' "noinst_HEADERS = " >> include_HEADERS
for header_with_path in $noinst_blacklist ; do
    echo " \\" >> include_HEADERS
    printf '%s' "        "$header_with_path >> include_HEADERS
done
echo >> include_HEADERS
echo >> include_HEADERS



headers=`find base enums error_estimation fe geom mesh numerics parallel partitioning timpi_shims physics quadrature reduced_basis solution_transfer solvers systems utils -name "*.h" -o -name "*specializations" -type f | LC_COLLATE=POSIX sort`
headers="libmesh_config.h $headers"

echo "# These are the headers we actually want to install and make available" >> include_HEADERS
echo "# to a user.  Consider these the public interface of libMesh." >> include_HEADERS
echo "# will get installed into '\$(prefix)/include/libmesh'" >> include_HEADERS
printf '%s' "include_HEADERS = " >> include_HEADERS
for header_with_path in $headers ; do
    add_header=1
    for blacklist in $noinst_blacklist ; do
	if (test "$blacklist" = "$header_with_path"); then
	    add_header=0
	fi
    done
    if (test $add_header -eq 1); then
	echo " \\" >> include_HEADERS
        printf '%s' "        "$header_with_path >> include_HEADERS
    fi
done
echo " " >> include_HEADERS
