AUTOMAKE_OPTIONS = subdir-objects 

AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CPPFLAGS  = $(libmesh_optional_INCLUDES) -I$(top_builddir)/include \
               $(libmesh_contrib_INCLUDES) $(CPPUNIT_CFLAGS) \
	       -DLIBMESH_IS_UNIT_TESTING
AM_LDFLAGS   = $(libmesh_LDFLAGS)
LIBS         = $(libmesh_optional_LIBS) $(CPPUNIT_LIBS)

unit_tests_sources = \
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

check_PROGRAMS = # empty, append below

if LIBMESH_OPT_MODE
  check_PROGRAMS         += unit_tests-opt
  unit_tests_opt_SOURCES  = $(unit_tests_sources)
  unit_tests_opt_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
  unit_tests_opt_CXXFLAGS = $(CXXFLAGS_OPT)
  unit_tests_opt_LDADD    = $(top_builddir)/libmesh_opt.la
endif

# our GLIBC debugging preprocessor flags break this.
if !LIBMESH_ENABLE_GLIBCXX_DEBUGGING
if LIBMESH_DBG_MODE
  check_PROGRAMS         += unit_tests-dbg
  unit_tests_dbg_SOURCES  = $(unit_tests_sources)
  unit_tests_dbg_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
  unit_tests_dbg_CXXFLAGS = $(CXXFLAGS_DBG)
  unit_tests_dbg_LDADD    = $(top_builddir)/libmesh_dbg.la
endif
endif

if LIBMESH_DEVEL_MODE
  check_PROGRAMS           += unit_tests-devel
  unit_tests_devel_SOURCES  = $(unit_tests_sources)
  unit_tests_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
  unit_tests_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
  unit_tests_devel_LDADD    = $(top_builddir)/libmesh_devel.la
endif

if LIBMESH_PROF_MODE
  check_PROGRAMS          += unit_tests-prof
  unit_tests_prof_SOURCES  = $(unit_tests_sources)
  unit_tests_prof_CPPFLAGS = $(CPPFLAGS_PROF) $(AM_CPPFLAGS)
  unit_tests_prof_CXXFLAGS = $(CXXFLAGS_PROF)
  unit_tests_prof_LDADD    = $(top_builddir)/libmesh_prof.la
endif

if LIBMESH_OPROF_MODE
  check_PROGRAMS           += unit_tests-oprof
  unit_tests_oprof_SOURCES  = $(unit_tests_sources)
  unit_tests_oprof_CPPFLAGS = $(CPPFLAGS_OPROF) $(AM_CPPFLAGS)
  unit_tests_oprof_CXXFLAGS = $(CXXFLAGS_OPROF)
  unit_tests_oprof_LDADD    = $(top_builddir)/libmesh_oprof.la
endif

TESTS = $(check_PROGRAMS)

######################################################################
#
# Don't leave code coverage outputs lying around
if CODE_COVERAGE_ENABLED
  CLEANFILES = */*.gcda */*.gcno
endif

