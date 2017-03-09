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
  test_comm.h \
  base/dof_object_test.h \
  base/default_coupling_test.C \
  base/getpot_test.C \
  base/point_neighbor_coupling_test.C \
  base/unique_ptr_test.C \
  fe/fe_bernstein_test.C \
  fe/fe_clough_test.C \
  fe/fe_hermite_test.C \
  fe/fe_hierarchic_test.C \
  fe/fe_l2_hierarchic_test.C \
  fe/fe_l2_lagrange_test.C \
  fe/fe_lagrange_test.C \
  fe/fe_monomial_test.C \
  fe/fe_szabab_test.C \
  fe/fe_test.h \
  fe/fe_xyz_test.C \
  geom/elem_test.C \
  geom/node_test.C \
  geom/point_test.C \
  geom/point_test.h \
  mesh/all_tri.C \
  mesh/boundary_mesh.C \
  mesh/boundary_info.C \
  mesh/checkpoint.C \
  mesh/contains_point.C \
  mesh/mixed_dim_mesh_test.C \
  mesh/nodal_neighbors.C \
  mesh/mesh_extruder.C \
  mesh/slit_mesh_test.C \
  mesh/spatial_dimension_test.C \
  mesh/mapped_subdomain_partitioner_test.C \
  mesh/mesh_function_dfem.C \
  numerics/composite_function_test.C \
  numerics/coupling_matrix_test.C \
  numerics/distributed_vector_test.C \
  numerics/eigen_sparse_vector_test.C \
  numerics/laspack_vector_test.C \
  numerics/numeric_vector_test.h \
  numerics/parsed_fem_function_test.C \
  numerics/parsed_function_test.C \
  numerics/petsc_vector_test.C \
  numerics/trilinos_epetra_vector_test.C \
  numerics/type_vector_test.h \
  numerics/vector_value_test.C \
  numerics/type_tensor_test.C \
  numerics/dense_matrix_test.C \
  parallel/packed_range_test.C \
  parallel/parallel_test.C \
  parallel/parallel_point_test.C \
  quadrature/quadrature_test.C \
  solvers/time_solver_test_common.h \
  solvers/first_order_unsteady_solver_test.C \
  solvers/second_order_unsteady_solver_test.C \
  systems/equation_systems_test.C \
  systems/systems_test.C \
  utils/point_locator_test.C \
  utils/vectormap_test.C

#EXTRA_DIST = base/getpot_test_input.in

if LIBMESH_ENABLE_FPARSER
  unit_tests_sources += \
	fparser/autodiff.C
endif

check_PROGRAMS = # empty, append below
noinst_PROGRAMS = # empty, append below

if LIBMESH_OPT_MODE
  check_PROGRAMS         += unit_tests-opt
  unit_tests_opt_SOURCES  = $(unit_tests_sources)
  unit_tests_opt_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
  unit_tests_opt_CXXFLAGS = $(CXXFLAGS_OPT)
  unit_tests_opt_LDADD    = $(top_builddir)/libmesh_opt.la
endif

# our GLIBC debugging preprocessor flags seem to potentially conflict
# with libcppunit binaries.  Some cppunit versions work fine for us,
# others segfault and/or hang.  By default we will not run
# GLIBCXX-debugging builds with cppunit unless specificallly
# configured to.
if LIBMESH_DBG_MODE
if LIBMESH_ENABLE_GLIBCXX_DEBUGGING
if LIBMESH_ENABLE_GLIBCXX_DEBUGGING_CPPUNIT
  check_PROGRAMS         += unit_tests-dbg
else
  noinst_PROGRAMS        += unit_tests-dbg
endif
else
  check_PROGRAMS         += unit_tests-dbg
endif
  unit_tests_dbg_SOURCES  = $(unit_tests_sources)
  unit_tests_dbg_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
  unit_tests_dbg_CXXFLAGS = $(CXXFLAGS_DBG)
  unit_tests_dbg_LDADD    = $(top_builddir)/libmesh_dbg.la
endif

if LIBMESH_DEVEL_MODE
  check_PROGRAMS           += unit_tests-devel
  noinst_PROGRAMS          += unit_tests-devel
  unit_tests_devel_SOURCES  = $(unit_tests_sources)
  unit_tests_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
  unit_tests_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
  unit_tests_devel_LDADD    = $(top_builddir)/libmesh_devel.la
endif

if LIBMESH_PROF_MODE
  check_PROGRAMS          += unit_tests-prof
  noinst_PROGRAMS         += unit_tests-prof
  unit_tests_prof_SOURCES  = $(unit_tests_sources)
  unit_tests_prof_CPPFLAGS = $(CPPFLAGS_PROF) $(AM_CPPFLAGS)
  unit_tests_prof_CXXFLAGS = $(CXXFLAGS_PROF)
  unit_tests_prof_LDADD    = $(top_builddir)/libmesh_prof.la
endif

if LIBMESH_OPROF_MODE
  check_PROGRAMS           += unit_tests-oprof
  noinst_PROGRAMS          += unit_tests-oprof
  unit_tests_oprof_SOURCES  = $(unit_tests_sources)
  unit_tests_oprof_CPPFLAGS = $(CPPFLAGS_OPROF) $(AM_CPPFLAGS)
  unit_tests_oprof_CXXFLAGS = $(CXXFLAGS_OPROF)
  unit_tests_oprof_LDADD    = $(top_builddir)/libmesh_oprof.la
endif


# Recursive automake builds subdirectories before parent directories.
# But here we need the subdirectory to be able to link to
# already-built parent directory libraries.
# For most of our modes, this means "you can't run make check unless
# make already completed", which is an annoying bug; for the dbg
# noinst case this means "you can't run make unless make already
# completed", which is a devastating bug.

# As a quick fix we'll ensure library dependencies are built manually:
# http://lists.gnu.org/archive/html/automake/2009-03/msg00011.html

FORCE:

.PHONY: FORCE

$(top_builddir)/libmesh_dbg.la: FORCE
	(cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) libmesh_dbg.la)

$(top_builddir)/libmesh_devel.la: FORCE
	(cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) libmesh_devel.la)

$(top_builddir)/libmesh_opt.la: FORCE
	(cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) libmesh_opt.la)

$(top_builddir)/libmesh_prof.la: FORCE
	(cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) libmesh_prof.la)

$(top_builddir)/libmesh_oprof.la: FORCE
	(cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) libmesh_oprof.la)



TESTS = run_unit_tests.sh

CLEANFILES = cube_mesh.xdr checkpoint_splitter.cpr-2-0 checkpoint_splitter.cpr-2-1 slit_mesh.xda slit_solution.xda

######################################################################
#
# Don't leave code coverage outputs lying around
if CODE_COVERAGE_ENABLED
  CLEANFILES += */*.gcda */*.gcno
endif
