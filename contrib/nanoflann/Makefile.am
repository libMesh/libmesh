AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = -I$(top_srcdir)/contrib/nanoflann/nanoflann/include
AM_CFLAGS    =

# Nanoflann uses C++ threads now, and sometimes we need special C++
# flags to get those to link properly.  We'll just use opt-mode here.
AM_CXXFLAGS  = $(CXXFLAGS_OPT)

# generally we do not install contributed headers, but
# for nanoflann we will because of its tight inlining
# in InverseDistanceInterpolation and friends
includedir      = $(prefix)/include/libmesh
include_HEADERS = $(srcdir)/nanoflann/include/nanoflann.hpp

######################################################################
#
# Building the examples
#
check_PROGRAMS =

check_PROGRAMS                += pointcloud_kdd_radius
pointcloud_kdd_radius_SOURCES  = $(srcdir)/nanoflann/examples/pointcloud_kdd_radius.cpp nanoflann/examples/utils.h

check_PROGRAMS             += pointcloud_example
pointcloud_example_SOURCES  = $(srcdir)/nanoflann/examples/pointcloud_example.cpp nanoflann/examples/utils.h

check_PROGRAMS                     += pointcloud_adaptor_example
pointcloud_adaptor_example_SOURCES  = $(srcdir)/nanoflann/examples/pointcloud_adaptor_example.cpp

check_PROGRAMS                    += vector_of_vectors_example
vector_of_vectors_example_SOURCES  = $(srcdir)/nanoflann/examples/vector_of_vectors_example.cpp nanoflann/examples/KDTreeVectorOfVectorsAdaptor.h

if LIBMESH_ENABLE_EIGEN
  AM_CPPFLAGS            += $(EIGEN_INCLUDE)
  check_PROGRAMS         += matrix_example
  matrix_example_SOURCES  = $(srcdir)/nanoflann/examples/matrix_example.cpp
endif

TESTS = $(check_PROGRAMS)
