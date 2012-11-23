include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_INCLUDE) -I$(top_builddir)/include \
               $(libmesh_contrib_INCLUDES) $(CPPUNIT_CFLAGS)
LIBS         = $(libmesh_LIBS) $(CPPUNIT_LIBS)

check_PROGRAMS = unit_tests

unit_tests_SOURCES = \
	driver.C \
	base/dof_object_test.h \
	geom/node_test.C \
	geom/point_test.C \
	geom/point_test.h \
	numerics/distributed_vector_test.C \
	numerics/laspack_vector_test.C \
	numerics/numeric_vector_test.h \
	numerics/petsc_vector_test.C \
	numerics/trilinos_epetra_vector_test.C \
	numerics/type_vector_test.h \
	parallel/parallel_test.C \
	quadrature/quadrature_test.C

unit_tests_LDADD = $(top_builddir)/libmesh.la

TESTS = $(check_PROGRAMS)
