AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CPPFLAGS  = $(libmesh_optional_INCLUDES) -I$(top_builddir)/include \
               $(libmesh_contrib_INCLUDES) $(CPPUNIT_CFLAGS) \
	       -DLIBMESH_IS_UNIT_TESTING
AM_LDFLAGS   = $(libmesh_LDFLAGS) $(libmesh_contrib_LDFLAGS)
LIBS         = $(libmesh_optional_LIBS) $(CPPUNIT_LIBS)

# We might have turned on -Werror and/or paranoid warnings
CXXFLAGS_DBG += $(ACSM_ANY_WERROR_FLAG) $(ACSM_ANY_PARANOID_FLAGS)
CXXFLAGS_DEVEL += $(ACSM_ANY_WERROR_FLAG) $(ACSM_ANY_PARANOID_FLAGS)
CXXFLAGS_OPROF += $(ACSM_ANY_WERROR_FLAG) $(ACSM_ANY_PARANOID_FLAGS)
CXXFLAGS_OPT += $(ACSM_ANY_WERROR_FLAG) $(ACSM_ANY_PARANOID_FLAGS)
CXXFLAGS_PROF += $(ACSM_ANY_WERROR_FLAG) $(ACSM_ANY_PARANOID_FLAGS)

unit_tests_sources = \
  driver.C \
  libmesh_cppunit.h \
  stream_redirector.h \
  test_comm.h \
  base/dof_object_test.h \
  base/dof_map_test.C \
  base/default_coupling_test.C \
  base/getpot_test.C \
  base/point_neighbor_coupling_test.C \
  base/overlapping_coupling_test.C \
  base/nonmanifold_coupling_test.C \
  base/multi_evaluable_pred_test.C \
  fe/fe_bernstein_test.C \
  fe/fe_clough_test.C \
  fe/fe_hermite_test.C \
  fe/fe_hierarchic_test.C \
  fe/inf_fe_radial_test.C \
  fe/fe_l2_hierarchic_test.C \
  fe/fe_l2_lagrange_test.C \
  fe/fe_lagrange_test.C \
  fe/fe_monomial_test.C \
  fe/fe_rational_map.C \
  fe/fe_rational_test.C \
  fe/fe_side_test.C \
  fe/fe_szabab_test.C \
  fe/fe_test.h \
  fe/fe_xyz_test.C \
  fe/dual_shape_verification_test.C \
  geom/bbox_test.C \
  geom/edge_test.C \
  geom/elem_test.C \
  geom/elem_test.h \
  geom/node_test.C \
  geom/point_test.C \
  geom/point_test.h \
  geom/side_test.C \
  geom/volume_test.C \
  geom/which_node_am_i_test.C \
  mesh/all_second_order.C \
  mesh/all_tri.C \
  mesh/distort.C \
  mesh/boundary_mesh.C \
  mesh/boundary_info.C \
  mesh/boundary_points.C \
  mesh/checkpoint.C \
  mesh/connected_components.C \
  mesh/contains_point.C \
  mesh/distributed_mesh_test.C \
  mesh/extra_integers.C \
  mesh/exodus_test.C \
  mesh/mesh_assign.C \
  mesh/mesh_base_test.C \
  mesh/mesh_collection.C \
  mesh/mesh_deletions.C \
  mesh/mesh_extruder.C \
  mesh/mesh_function.C \
  mesh/mesh_function_dfem.C \
  mesh/mesh_generation_test.C \
  mesh/mesh_input.C \
  mesh/mesh_smoother_test.C \
  mesh/mesh_stitch.C \
  mesh/mesh_tet_test.C \
  mesh/mesh_triangulation.C \
  mesh/mixed_dim_mesh_test.C \
  mesh/mixed_order_test.C \
  mesh/nodal_neighbors.C \
  mesh/libmesh_poly2tri.C \
  mesh/libmesh_netgen.C \
  mesh/simplex_refinement_test.C \
  mesh/slit_mesh_test.C \
  mesh/spatial_dimension_test.C \
  mesh/mapped_subdomain_partitioner_test.C \
  mesh/write_elemset_data.C \
  mesh/write_sideset_data.C \
  mesh/write_nodeset_data.C \
  mesh/write_edgeset_data.C \
  mesh/write_vec_and_scalar.C \
  numerics/composite_function_test.C \
  numerics/coupling_matrix_test.C \
  numerics/distributed_vector_test.C \
  numerics/eigen_sparse_vector_test.C \
  numerics/laspack_vector_test.C \
  numerics/laspack_matrix_test.C \
  numerics/numeric_vector_test.h \
  numerics/parsed_fem_function_test.C \
  numerics/parsed_function_test.C \
  numerics/petsc_vector_test.C \
  numerics/trilinos_epetra_vector_test.C \
  numerics/type_vector_test.h \
  numerics/vector_value_test.C \
  numerics/type_tensor_test.C \
  numerics/sparse_matrix_test.h \
  numerics/dense_matrix_test.C \
  numerics/petsc_matrix_test.C \
  numerics/diagonal_matrix_test.C \
  numerics/lumped_mass_matrix_test.C \
  numerics/eigen_sparse_matrix_test.C \
  numerics/tensor_traits_test.C \
  parallel/message_tag.C \
  parallel/packed_range_test.C \
  parallel/packing_types_test.C \
  parallel/parallel_sort_test.C \
  parallel/parallel_sync_test.C \
  parallel/parallel_ghost_sync_test.C \
  parallel/parallel_test.C \
  parallel/parallel_point_test.C \
  partitioning/partitioner_test.h \
  partitioning/centroid_partitioner_test.C \
  partitioning/hilbert_sfc_partitioner_test.C \
  partitioning/linear_partitioner_test.C \
  partitioning/metis_partitioner_test.C \
  partitioning/morton_sfc_partitioner_test.C \
  partitioning/parmetis_partitioner_test.C \
  partitioning/sfc_partitioner_test.C \
  quadrature/quadrature_test.C \
  solvers/time_solver_test_common.h \
  solvers/first_order_unsteady_solver_test.C \
  solvers/second_order_unsteady_solver_test.C \
  systems/constraint_operator_test.C \
  systems/equation_systems_test.C \
  systems/periodic_bc_test.C \
  systems/systems_test.C \
  utils/parameters_test.C \
  utils/point_locator_test.C \
  utils/rb_parameters_test.C \
  utils/transparent_comparator.C \
  utils/vectormap_test.C \
  utils/xdr_test.C

# Make sure any data goes in "make dist"
datadir = $(install_dir)
data = matrices/geom_1_extraction_op.m \
       matrices/geom_1_extraction_op.m.gz \
       matrices/geom_2_extraction_op.m \
       matrices/geom_2_extraction_op.m.gz \
       matrices/geom_4_extraction_op.m \
       matrices/geom_4_extraction_op.m.gz \
       matrices/geom_8_extraction_op.m \
       matrices/geom_1_extraction_op.petsc32 \
       matrices/geom_1_extraction_op.petsc64 \
       meshes/1_quad.bxt.gz \
       meshes/25_quad.bxt.gz \
       meshes/bad_64bit_elem_integers.e \
       meshes/block.msh \
       meshes/BlockWithHole_Patch9.bxt.gz \
       meshes/circle.msh \
       meshes/constrain10to4.m \
       meshes/Cube_With_Sidesets.e \
       meshes/exodus_elements/read_exodus_EDGE2.e \
       meshes/exodus_elements/read_exodus_EDGE3.e \
       meshes/exodus_elements/read_exodus_EDGE4.e \
       meshes/exodus_elements/read_exodus_HEX20.e \
       meshes/exodus_elements/read_exodus_HEX27.e \
       meshes/exodus_elements/read_exodus_HEX8.e \
       meshes/exodus_elements/read_exodus_PRISM15.e \
       meshes/exodus_elements/read_exodus_PRISM18.e \
       meshes/exodus_elements/read_exodus_PRISM20.e \
       meshes/exodus_elements/read_exodus_PRISM21.e \
       meshes/exodus_elements/read_exodus_PRISM6.e \
       meshes/exodus_elements/read_exodus_PYRAMID13.e \
       meshes/exodus_elements/read_exodus_PYRAMID14.e \
       meshes/exodus_elements/read_exodus_PYRAMID18.e \
       meshes/exodus_elements/read_exodus_PYRAMID5.e \
       meshes/exodus_elements/read_exodus_QUAD4.e \
       meshes/exodus_elements/read_exodus_QUADSHELL4.e \
       meshes/exodus_elements/read_exodus_QUAD8.e \
       meshes/exodus_elements/read_exodus_QUADSHELL8.e \
       meshes/exodus_elements/read_exodus_QUAD9.e \
       meshes/exodus_elements/read_exodus_QUADSHELL9.e \
       meshes/exodus_elements/read_exodus_TET10.e \
       meshes/exodus_elements/read_exodus_TET14.e \
       meshes/exodus_elements/read_exodus_TET4.e \
       meshes/exodus_elements/read_exodus_TRI3.e \
       meshes/exodus_elements/read_exodus_TRISHELL3.e \
       meshes/exodus_elements/read_exodus_TRI6.e \
       meshes/exodus_elements/read_exodus_TRI7.e \
       meshes/geom_2.exo \
       meshes/geom_4.exo \
       meshes/geom_8.exo \
       meshes/good_32bit_elem_integers.e \
       meshes/mesh_assign_test_mesh.xda \
       meshes/mesh_with_low_order_edge_blocks.e \
       meshes/PlateWithHole_Patch8.bxt.gz \
       meshes/PlateWithHole_Patch8.e \
       meshes/poly_sample_test2.inp.gz \
       meshes/PressurizedCyl3d_Patch1_8Elem.bxt.gz \
       meshes/PressurizedCyl3d_Patch1_8Elem.e \
       meshes/PressurizedCyl_Patch6_256Elem.bxt.gz \
       meshes/shark_tooth_tri6.xda.gz \
       meshes/tensile_sample_test1.inp.gz \
       meshes/tetgen_one_tet10.ele \
       meshes/tetgen_one_tet10.node \
       meshes/two_element_iga_in.e \
       meshes/two_quads_two_blocks.e \
       meshes/non_manifold_junction0.exo \
       meshes/non_manifold_junction1.exo \
       meshes/non_manifold_junction2.exo \
       meshes/non_manifold_junction3.exo \
       meshes/Cluster_34.stl \
       meshes/engraving.stl


unit_tests_data = $(data)

# Why isn't this working automatically?
EXTRA_DIST = $(data)

if LIBMESH_ENABLE_FPARSER
  unit_tests_sources += \
	fparser/autodiff.C
endif

check_PROGRAMS = # empty, append below

# our GLIBC debugging preprocessor flags seem to potentially conflict
# with libcppunit binaries.  Some cppunit versions work fine for us,
# others segfault and/or hang.  By default we will not run
# GLIBCXX-debugging builds with cppunit unless specifically
# configured to.
if LIBMESH_ENABLE_CPPUNIT
if LIBMESH_DBG_MODE
if ACSM_ENABLE_GLIBCXX_DEBUGGING
if ACSM_ENABLE_GLIBCXX_DEBUGGING_CPPUNIT
  check_PROGRAMS         += unit_tests-dbg
endif
else
  check_PROGRAMS         += unit_tests-dbg
endif
  unit_tests_dbg_SOURCES  = $(unit_tests_sources)
  unit_tests_dbg_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
  unit_tests_dbg_CXXFLAGS = $(CXXFLAGS_DBG)
  unit_tests_dbg_LDADD    = $(top_builddir)/libmesh_dbg.la
  unit_tests_dbgdir       = $(datadir)
  unit_tests_dbg_DATA     = $(data)
endif
endif

if LIBMESH_ENABLE_CPPUNIT
if LIBMESH_DEVEL_MODE
  check_PROGRAMS           += unit_tests-devel
  unit_tests_devel_SOURCES  = $(unit_tests_sources)
  unit_tests_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
  unit_tests_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
  unit_tests_devel_LDADD    = $(top_builddir)/libmesh_devel.la
  unit_tests_develdir       = $(datadir)
  unit_tests_devel_DATA     = $(data)
endif
endif

if LIBMESH_ENABLE_CPPUNIT
if LIBMESH_PROF_MODE
  check_PROGRAMS          += unit_tests-prof
  unit_tests_prof_SOURCES  = $(unit_tests_sources)
  unit_tests_prof_CPPFLAGS = $(CPPFLAGS_PROF) $(AM_CPPFLAGS)
  unit_tests_prof_CXXFLAGS = $(CXXFLAGS_PROF)
  unit_tests_prof_LDADD    = $(top_builddir)/libmesh_prof.la
  unit_tests_profdir       = $(datadir)
  unit_tests_prof_DATA     = $(data)
endif
endif

if LIBMESH_ENABLE_CPPUNIT
if LIBMESH_OPROF_MODE
  check_PROGRAMS           += unit_tests-oprof
  unit_tests_oprof_SOURCES  = $(unit_tests_sources)
  unit_tests_oprof_CPPFLAGS = $(CPPFLAGS_OPROF) $(AM_CPPFLAGS)
  unit_tests_oprof_CXXFLAGS = $(CXXFLAGS_OPROF)
  unit_tests_oprof_LDADD    = $(top_builddir)/libmesh_oprof.la
  unit_tests_oprofdir       = $(datadir)
  unit_tests_oprof_DATA     = $(data)
endif
endif

if LIBMESH_ENABLE_CPPUNIT
if LIBMESH_OPT_MODE
  check_PROGRAMS         += unit_tests-opt
  unit_tests_opt_SOURCES  = $(unit_tests_sources)
  unit_tests_opt_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
  unit_tests_opt_CXXFLAGS = $(CXXFLAGS_OPT)
  unit_tests_opt_LDADD    = $(top_builddir)/libmesh_opt.la
  unit_tests_optdir       = $(datadir)
  unit_tests_opt_DATA     = $(data)
endif
endif

# Recursive automake builds subdirectories before parent directories.
# But here we need the subdirectory to be able to link to
# already-built parent directory libraries.
# For most of our modes, this means "you can't run make check unless
# make already completed", which is an annoying bug; for the dbg
# noinst case this means "you can't run make unless make already
# completed", which is a devastating bug.
#
# Putting . in the top-level SUBDIRS works when make is invoked on the
# top level, but can leave `make -C tests` runs linked (or failing to
# link, if there have been ABI changes) against a non-up-to-date
# library file.
#
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



if LIBMESH_ENABLE_CPPUNIT
TESTS = run_unit_tests.sh
endif

CLEANFILES = cube_mesh.xda \
             slit_mesh.xda \
             slit_solution.xda \
             out.e \
             mesh_with_soln.e \
             elemental_from_nodal.e \
             write_elemset_data.e \
             write_sideset_data.e \
             write_nodeset_data.e \
             write_edgeset_data.e \
             read_header_test.e \
             output.dat \
             periodic.e \
             mesh.out.memory.moved.xda \
             mesh.out.file.moved.xda \
             repl_with_nodal_soln.nem.1.0 \
             dist_with_elem_soln.e \
             dist.e \
             dist_with_nodal_soln.e \
             dist_with_nodal_soln.nem.1.0 \
             dist_with_elem_vec.e \
             dist_with_single_elem.nem.1.0 \
             dist_with_elem_soln.nem.1.0 \
             rep.nem.1.0 \
             repl_with_elem_soln.nem.1.0 \
             repl_with_single_elem.nem.1.0 \
             dist_with_elem_vec.nem.1.0 \
             repl_with_elem_vec.nem.1.0 \
             dist.nem.1.0 \
             repl_with_elem_soln.e \
             test_nemesis_read.nem.1.0 \
             rep.e \
             repl_with_nodal_soln.e \
             repl_with_elem_vec.e \
             extra_integers.xda \
             extra_integers.xdr \
             write_elemset_data.xda \
             write_elemset_data.xdr \
             test_helper.xdr \
             SparseMatrixTest_matrix.m \
             geom_1_extraction_op.petsc \
             geom_1_extraction_op.petsc.info \
             side_discontinuous_EDGE3.e \
             side_discontinuous_EDGE3_disc.e \
             side_discontinuous_HEX27.e \
             side_discontinuous_HEX27_disc.e \
             side_discontinuous_QUAD9.e \
             side_discontinuous_QUAD9_disc.e \
             side_discontinuous_TRI6.e \
             side_discontinuous_TRI6_disc.e \
             write_exodus_EDGE2.e \
             write_exodus_EDGE3.e \
             write_exodus_EDGE4.e \
             write_exodus_HEX20.e \
             write_exodus_HEX27.e \
             write_exodus_HEX8.e \
             write_exodus_PRISM15.e \
             write_exodus_PRISM18.e \
             write_exodus_PRISM20.e \
             write_exodus_PRISM21.e \
             write_exodus_PRISM6.e \
             write_exodus_PYRAMID13.e \
             write_exodus_PYRAMID14.e \
             write_exodus_PYRAMID18.e \
             write_exodus_PYRAMID5.e \
             write_exodus_QUAD4.e \
             write_exodus_QUAD8.e \
             write_exodus_QUAD9.e \
             write_exodus_QUADSHELL4.e \
             write_exodus_QUADSHELL8.e \
             write_exodus_QUADSHELL9.e \
             write_exodus_TET10.e \
             write_exodus_TET14.e \
             write_exodus_TET4.e \
             write_exodus_TRI3.e \
             write_exodus_TRI6.e \
             write_exodus_TRI7.e \
             write_exodus_TRISHELL3.e

# need to link any data files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp

.linkstamp:
	-rm -f meshes && $(LN_S) -f $(srcdir)/meshes .
	-rm -f matrices && $(LN_S) -f $(srcdir)/matrices .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += .linkstamp
endif

# CLEANFILES only applies to files, handle directories separately
clean-local:
	rm -rf checkpoint_splitter.cpr \
	       checkpoint_splitter.cpa \
               extra_integers.cpr \
               extra_integers.cpa \
	       .jitcache

######################################################################
#
# Don't leave code coverage outputs lying around
if CODE_COVERAGE_ENABLED
  CLEANFILES += */*.gcda */*.gcno
endif
