AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = -I$(srcdir)/include
AM_CFLAGS    =
AM_CXXFLAGS  =

# generally we do not install contributed headers, but
# for nanoflann we will because of its tight inlining
# in InverseDistanceInterpolation and friends
includedir      = $(prefix)/include/libmesh
include_HEADERS = include/nanoflann.hpp

######################################################################
#
# Building the examples
#
check_PROGRAMS =

check_PROGRAMS                += pointcloud_kdd_radius
pointcloud_kdd_radius_SOURCES  = examples/pointcloud_kdd_radius.cpp examples/utils.h

check_PROGRAMS             += pointcloud_example
pointcloud_example_SOURCES  = examples/pointcloud_example.cpp examples/utils.h

check_PROGRAMS                     += pointcloud_adaptor_example
pointcloud_adaptor_example_SOURCES  = examples/pointcloud_adaptor_example.cpp

check_PROGRAMS                    += vector_of_vectors_example
vector_of_vectors_example_SOURCES  = examples/vector_of_vectors_example.cpp examples/KDTreeVectorOfVectorsAdaptor.h

if LIBMESH_ENABLE_EIGEN
  AM_CPPFLAGS            += $(EIGEN_INCLUDE)
  check_PROGRAMS         += matrix_example
  matrix_example_SOURCES  = examples/matrix_example.cpp
endif

TESTS = $(check_PROGRAMS)
