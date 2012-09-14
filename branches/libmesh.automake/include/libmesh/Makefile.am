# special handholding for prefix_config.m4 generated files
# so that 'make clean ; make' works as does 'make distcheck'
# libmesh_config.h is made by ./configure, so it should get 
# cleaned by 'make distclean'
DISTCLEANFILES = libmesh_config.h

#
# special things to do when running 'make dist'
dist-hook:
	rm -rf $(distdir)/libmesh_config.h

#
# include the magic script!
EXTRA_DIST = rebuild_makefile.sh

Makefile.am: rebuild_makefile.sh
	$(AM_V_GEN)./rebuild_makefile.sh

BUILT_SOURCES =  auto_ptr.h dirichlet_boundaries.h dof_map.h dof_object.h factory.h getpot.h id_types.h libmesh_base.h libmesh_C_isnan.h libmesh_common.h libmesh_documentation.h libmesh_exceptions.h libmesh.h libmesh_logging.h multi_predicates.h periodic_boundaries.h periodic_boundary_base.h periodic_boundary.h print_trace.h reference_counted_object.h reference_counter.h single_predicates.h sparsity_pattern.h variable.h variant_filter_iterator.h enum_eigen_solver_type.h enum_elem_quality.h enum_elem_type.h enum_fe_family.h enum_inf_map_type.h enum_io_package.h enum_norm_type.h enum_order.h enum_parallel_type.h enum_point_locator_type.h enum_preconditioner_type.h enum_quadrature_type.h enum_solver_package.h enum_solver_type.h enum_subset_solve_mode.h enum_xdr_mode.h adjoint_refinement_estimator.h adjoint_residual_error_estimator.h discontinuity_measure.h error_estimator.h exact_error_estimator.h exact_solution.h fourth_error_estimators.h hp_coarsentest.h hp_selector.h hp_singular.h jump_error_estimator.h kelly_error_estimator.h patch_recovery_error_estimator.h uniform_refinement_estimator.h fe_abstract.h fe_base.h fe_compute_data.h fe.h fe_interface.h fe_macro.h fe_map.h fe_transformation_base.h fe_type.h fe_xyz_map.h h1_fe_transformation.h hcurl_fe_transformation.h inf_fe.h inf_fe_instantiate_1D.h inf_fe_instantiate_2D.h inf_fe_instantiate_3D.h inf_fe_macro.h cell.h cell_hex20.h cell_hex27.h cell_hex8.h cell_hex.h cell_inf.h cell_inf_hex16.h cell_inf_hex18.h cell_inf_hex8.h cell_inf_hex.h cell_inf_prism12.h cell_inf_prism6.h cell_inf_prism.h cell_prism15.h cell_prism18.h cell_prism6.h cell_prism.h cell_pyramid5.h cell_pyramid.h cell_tet10.h cell_tet4.h cell_tet.h edge_edge2.h edge_edge3.h edge_edge4.h edge.h edge_inf_edge2.h elem.h elem_quality.h elem_range.h elem_type.h face.h face_inf_quad4.h face_inf_quad6.h face_inf_quad.h face_quad4.h face_quad8.h face_quad9.h face_quad.h face_tri3.h face_tri6.h face_tri.h node_elem.h node.h node_range.h plane.h point.h remote_elem.h side.h sphere.h stored_range.h surface.h abaqus_io.h boundary_info.h boundary_mesh.h diva_io.h ensight_io.h exodusII_io.h exodusII_io_helper.h fro_io.h gmsh_io.h gmv_io.h gnuplot_io.h inf_elem_builder.h legacy_xdr_io.h matlab_io.h medit_io.h mesh_base.h mesh_communication.h mesh_data.h mesh_function.h mesh_generation.h mesh.h mesh_input.h mesh_inserter_iterator.h mesh_modification.h mesh_output.h mesh_refinement.h mesh_serializer.h mesh_smoother.h mesh_smoother_laplace.h mesh_smoother_vsmoother.h mesh_tetgen_interface.h mesh_tetgen_wrapper.h mesh_tools.h mesh_triangle_holes.h mesh_triangle_interface.h mesh_triangle_wrapper.h nemesis_io.h nemesis_io_helper.h off_io.h parallel_mesh.h patch.h postscript_io.h serial_mesh.h tecplot_io.h tetgen_io.h ucd_io.h unstructured_mesh.h unv_io.h vtk_io.h xdr_head.h xdr_io.h xdr_mesh.h xdr_mgf.h xdr_mhead.h xdr_shead.h xdr_soln.h analytic_function.h const_function.h coupling_matrix.h dense_matrix_base.h dense_matrix.h dense_submatrix.h dense_subvector.h dense_vector_base.h dense_vector.h distributed_vector.h function_base.h laspack_matrix.h laspack_vector.h numeric_vector.h parsed_function.h petsc_macro.h petsc_matrix.h petsc_preconditioner.h petsc_vector.h preconditioner.h raw_accessor.h refinement_selector.h shell_matrix.h sparse_matrix.h sparse_shell_matrix.h sum_shell_matrix.h tensor_shell_matrix.h tensor_tools.h tensor_value.h trilinos_epetra_matrix.h trilinos_epetra_vector.h trilinos_preconditioner.h type_n_tensor.h type_tensor.h type_vector.h vector_value.h wrapped_function.h zero_function.h parallel_algebra.h parallel_bin_sorter.h parallel_conversion_utils.h parallel_ghost_sync.h parallel.h parallel_hilbert.h parallel_histogram.h parallel_sort.h threads_allocators.h threads.h centroid_partitioner.h hilbert_sfc_partitioner.h linear_partitioner.h metis_partitioner.h morton_sfc_partitioner.h parmetis_partitioner.h partitioner.h sfc_partitioner.h diff_physics.h diff_qoi.h quadrature_clough.h quadrature_conical.h quadrature_gauss.h quadrature_gm.h quadrature_grid.h quadrature.h quadrature_jacobi.h quadrature_monomial.h quadrature_rules.h quadrature_simpson.h quadrature_trap.h derived_rb_construction.h derived_rb_evaluation.h rb_assembly_expansion.h rb_construction_base.h rb_construction.h rb_eim_assembly.h rb_eim_construction.h rb_eim_evaluation.h rb_eim_theta.h rb_evaluation.h rb_parameters.h rb_parametrized_function.h rb_parametrized.h rb_scm_construction.h rb_scm_evaluation.h rb_temporal_discretization.h rb_theta_expansion.h rb_theta.h transient_rb_assembly_expansion.h transient_rb_construction.h transient_rb_evaluation.h transient_rb_theta_expansion.h adaptive_time_solver.h diff_solver.h eigen_solver.h eigen_time_solver.h euler2_solver.h euler_solver.h laspack_linear_solver.h linear.h linear_solver.h newton_solver.h nonlinear_solver.h petsc_diff_solver.h petscdmlibmesh.h petsc_dm_nonlinear_solver.h petsc_linear_solver.h petsc_nonlinear_solver.h slepc_eigen_solver.h slepc_macro.h solver.h steady_solver.h time_solver.h trilinos_aztec_linear_solver.h trilinos_nox_nonlinear_solver.h twostep_time_solver.h unsteady_solver.h condensed_eigen_system.h continuation_system.h diff_context.h diff_system.h eigen_system.h elem_assembly.h equation_systems.h explicit_system.h fem_context.h fem_system.h frequency_system.h implicit_system.h linear_implicit_system.h newmark_system.h nonlinear_implicit_system.h parameter_vector.h qoi_set.h sensitivity_data.h steady_system.h system.h system_norm.h system_subset_by_subdomain.h system_subset.h transient_system.h compare_types.h error_vector.h hashword.h ignore_warnings.h location_maps.h mapvector.h null_output_iterator.h number_lookups.h o_f_stream.h ostream_proxy.h o_string_stream.h parameters.h perf_log.h perfmon.h plt_loader.h point_locator_base.h point_locator_list.h point_locator_tree.h pool_allocator.h restore_warnings.h statistics.h string_to_enum.h timestamp.h tree_base.h tree.h tree_node.h utility.h xdr_cxx.h
CLEANFILES = $(BUILT_SOURCES)

auto_ptr.h: ../base/auto_ptr.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/auto_ptr.h auto_ptr.h

dirichlet_boundaries.h: ../base/dirichlet_boundaries.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/dirichlet_boundaries.h dirichlet_boundaries.h

dof_map.h: ../base/dof_map.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/dof_map.h dof_map.h

dof_object.h: ../base/dof_object.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/dof_object.h dof_object.h

factory.h: ../base/factory.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/factory.h factory.h

getpot.h: ../base/getpot.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/getpot.h getpot.h

id_types.h: ../base/id_types.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/id_types.h id_types.h

libmesh_base.h: ../base/libmesh_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh_base.h libmesh_base.h

libmesh_C_isnan.h: ../base/libmesh_C_isnan.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh_C_isnan.h libmesh_C_isnan.h

libmesh_common.h: ../base/libmesh_common.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh_common.h libmesh_common.h

libmesh_documentation.h: ../base/libmesh_documentation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh_documentation.h libmesh_documentation.h

libmesh_exceptions.h: ../base/libmesh_exceptions.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh_exceptions.h libmesh_exceptions.h

libmesh.h: ../base/libmesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh.h libmesh.h

libmesh_logging.h: ../base/libmesh_logging.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/libmesh_logging.h libmesh_logging.h

multi_predicates.h: ../base/multi_predicates.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/multi_predicates.h multi_predicates.h

periodic_boundaries.h: ../base/periodic_boundaries.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/periodic_boundaries.h periodic_boundaries.h

periodic_boundary_base.h: ../base/periodic_boundary_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/periodic_boundary_base.h periodic_boundary_base.h

periodic_boundary.h: ../base/periodic_boundary.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/periodic_boundary.h periodic_boundary.h

print_trace.h: ../base/print_trace.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/print_trace.h print_trace.h

reference_counted_object.h: ../base/reference_counted_object.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/reference_counted_object.h reference_counted_object.h

reference_counter.h: ../base/reference_counter.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/reference_counter.h reference_counter.h

single_predicates.h: ../base/single_predicates.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/single_predicates.h single_predicates.h

sparsity_pattern.h: ../base/sparsity_pattern.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/sparsity_pattern.h sparsity_pattern.h

variable.h: ../base/variable.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/variable.h variable.h

variant_filter_iterator.h: ../base/variant_filter_iterator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/base/variant_filter_iterator.h variant_filter_iterator.h

enum_eigen_solver_type.h: ../enums/enum_eigen_solver_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_eigen_solver_type.h enum_eigen_solver_type.h

enum_elem_quality.h: ../enums/enum_elem_quality.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_elem_quality.h enum_elem_quality.h

enum_elem_type.h: ../enums/enum_elem_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_elem_type.h enum_elem_type.h

enum_fe_family.h: ../enums/enum_fe_family.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_fe_family.h enum_fe_family.h

enum_inf_map_type.h: ../enums/enum_inf_map_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_inf_map_type.h enum_inf_map_type.h

enum_io_package.h: ../enums/enum_io_package.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_io_package.h enum_io_package.h

enum_norm_type.h: ../enums/enum_norm_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_norm_type.h enum_norm_type.h

enum_order.h: ../enums/enum_order.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_order.h enum_order.h

enum_parallel_type.h: ../enums/enum_parallel_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_parallel_type.h enum_parallel_type.h

enum_point_locator_type.h: ../enums/enum_point_locator_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_point_locator_type.h enum_point_locator_type.h

enum_preconditioner_type.h: ../enums/enum_preconditioner_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_preconditioner_type.h enum_preconditioner_type.h

enum_quadrature_type.h: ../enums/enum_quadrature_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_quadrature_type.h enum_quadrature_type.h

enum_solver_package.h: ../enums/enum_solver_package.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_solver_package.h enum_solver_package.h

enum_solver_type.h: ../enums/enum_solver_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_solver_type.h enum_solver_type.h

enum_subset_solve_mode.h: ../enums/enum_subset_solve_mode.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_subset_solve_mode.h enum_subset_solve_mode.h

enum_xdr_mode.h: ../enums/enum_xdr_mode.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/enums/enum_xdr_mode.h enum_xdr_mode.h

adjoint_refinement_estimator.h: ../error_estimation/adjoint_refinement_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/adjoint_refinement_estimator.h adjoint_refinement_estimator.h

adjoint_residual_error_estimator.h: ../error_estimation/adjoint_residual_error_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/adjoint_residual_error_estimator.h adjoint_residual_error_estimator.h

discontinuity_measure.h: ../error_estimation/discontinuity_measure.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/discontinuity_measure.h discontinuity_measure.h

error_estimator.h: ../error_estimation/error_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/error_estimator.h error_estimator.h

exact_error_estimator.h: ../error_estimation/exact_error_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/exact_error_estimator.h exact_error_estimator.h

exact_solution.h: ../error_estimation/exact_solution.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/exact_solution.h exact_solution.h

fourth_error_estimators.h: ../error_estimation/fourth_error_estimators.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/fourth_error_estimators.h fourth_error_estimators.h

hp_coarsentest.h: ../error_estimation/hp_coarsentest.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/hp_coarsentest.h hp_coarsentest.h

hp_selector.h: ../error_estimation/hp_selector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/hp_selector.h hp_selector.h

hp_singular.h: ../error_estimation/hp_singular.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/hp_singular.h hp_singular.h

jump_error_estimator.h: ../error_estimation/jump_error_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/jump_error_estimator.h jump_error_estimator.h

kelly_error_estimator.h: ../error_estimation/kelly_error_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/kelly_error_estimator.h kelly_error_estimator.h

patch_recovery_error_estimator.h: ../error_estimation/patch_recovery_error_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/patch_recovery_error_estimator.h patch_recovery_error_estimator.h

uniform_refinement_estimator.h: ../error_estimation/uniform_refinement_estimator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/error_estimation/uniform_refinement_estimator.h uniform_refinement_estimator.h

fe_abstract.h: ../fe/fe_abstract.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_abstract.h fe_abstract.h

fe_base.h: ../fe/fe_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_base.h fe_base.h

fe_compute_data.h: ../fe/fe_compute_data.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_compute_data.h fe_compute_data.h

fe.h: ../fe/fe.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe.h fe.h

fe_interface.h: ../fe/fe_interface.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_interface.h fe_interface.h

fe_macro.h: ../fe/fe_macro.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_macro.h fe_macro.h

fe_map.h: ../fe/fe_map.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_map.h fe_map.h

fe_transformation_base.h: ../fe/fe_transformation_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_transformation_base.h fe_transformation_base.h

fe_type.h: ../fe/fe_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_type.h fe_type.h

fe_xyz_map.h: ../fe/fe_xyz_map.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/fe_xyz_map.h fe_xyz_map.h

h1_fe_transformation.h: ../fe/h1_fe_transformation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/h1_fe_transformation.h h1_fe_transformation.h

hcurl_fe_transformation.h: ../fe/hcurl_fe_transformation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/hcurl_fe_transformation.h hcurl_fe_transformation.h

inf_fe.h: ../fe/inf_fe.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/inf_fe.h inf_fe.h

inf_fe_instantiate_1D.h: ../fe/inf_fe_instantiate_1D.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/inf_fe_instantiate_1D.h inf_fe_instantiate_1D.h

inf_fe_instantiate_2D.h: ../fe/inf_fe_instantiate_2D.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/inf_fe_instantiate_2D.h inf_fe_instantiate_2D.h

inf_fe_instantiate_3D.h: ../fe/inf_fe_instantiate_3D.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/inf_fe_instantiate_3D.h inf_fe_instantiate_3D.h

inf_fe_macro.h: ../fe/inf_fe_macro.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/fe/inf_fe_macro.h inf_fe_macro.h

cell.h: ../geom/cell.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell.h cell.h

cell_hex20.h: ../geom/cell_hex20.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_hex20.h cell_hex20.h

cell_hex27.h: ../geom/cell_hex27.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_hex27.h cell_hex27.h

cell_hex8.h: ../geom/cell_hex8.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_hex8.h cell_hex8.h

cell_hex.h: ../geom/cell_hex.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_hex.h cell_hex.h

cell_inf.h: ../geom/cell_inf.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf.h cell_inf.h

cell_inf_hex16.h: ../geom/cell_inf_hex16.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_hex16.h cell_inf_hex16.h

cell_inf_hex18.h: ../geom/cell_inf_hex18.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_hex18.h cell_inf_hex18.h

cell_inf_hex8.h: ../geom/cell_inf_hex8.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_hex8.h cell_inf_hex8.h

cell_inf_hex.h: ../geom/cell_inf_hex.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_hex.h cell_inf_hex.h

cell_inf_prism12.h: ../geom/cell_inf_prism12.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_prism12.h cell_inf_prism12.h

cell_inf_prism6.h: ../geom/cell_inf_prism6.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_prism6.h cell_inf_prism6.h

cell_inf_prism.h: ../geom/cell_inf_prism.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_inf_prism.h cell_inf_prism.h

cell_prism15.h: ../geom/cell_prism15.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_prism15.h cell_prism15.h

cell_prism18.h: ../geom/cell_prism18.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_prism18.h cell_prism18.h

cell_prism6.h: ../geom/cell_prism6.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_prism6.h cell_prism6.h

cell_prism.h: ../geom/cell_prism.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_prism.h cell_prism.h

cell_pyramid5.h: ../geom/cell_pyramid5.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_pyramid5.h cell_pyramid5.h

cell_pyramid.h: ../geom/cell_pyramid.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_pyramid.h cell_pyramid.h

cell_tet10.h: ../geom/cell_tet10.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_tet10.h cell_tet10.h

cell_tet4.h: ../geom/cell_tet4.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_tet4.h cell_tet4.h

cell_tet.h: ../geom/cell_tet.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/cell_tet.h cell_tet.h

edge_edge2.h: ../geom/edge_edge2.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/edge_edge2.h edge_edge2.h

edge_edge3.h: ../geom/edge_edge3.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/edge_edge3.h edge_edge3.h

edge_edge4.h: ../geom/edge_edge4.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/edge_edge4.h edge_edge4.h

edge.h: ../geom/edge.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/edge.h edge.h

edge_inf_edge2.h: ../geom/edge_inf_edge2.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/edge_inf_edge2.h edge_inf_edge2.h

elem.h: ../geom/elem.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/elem.h elem.h

elem_quality.h: ../geom/elem_quality.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/elem_quality.h elem_quality.h

elem_range.h: ../geom/elem_range.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/elem_range.h elem_range.h

elem_type.h: ../geom/elem_type.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/elem_type.h elem_type.h

face.h: ../geom/face.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face.h face.h

face_inf_quad4.h: ../geom/face_inf_quad4.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_inf_quad4.h face_inf_quad4.h

face_inf_quad6.h: ../geom/face_inf_quad6.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_inf_quad6.h face_inf_quad6.h

face_inf_quad.h: ../geom/face_inf_quad.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_inf_quad.h face_inf_quad.h

face_quad4.h: ../geom/face_quad4.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_quad4.h face_quad4.h

face_quad8.h: ../geom/face_quad8.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_quad8.h face_quad8.h

face_quad9.h: ../geom/face_quad9.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_quad9.h face_quad9.h

face_quad.h: ../geom/face_quad.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_quad.h face_quad.h

face_tri3.h: ../geom/face_tri3.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_tri3.h face_tri3.h

face_tri6.h: ../geom/face_tri6.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_tri6.h face_tri6.h

face_tri.h: ../geom/face_tri.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/face_tri.h face_tri.h

node_elem.h: ../geom/node_elem.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/node_elem.h node_elem.h

node.h: ../geom/node.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/node.h node.h

node_range.h: ../geom/node_range.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/node_range.h node_range.h

plane.h: ../geom/plane.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/plane.h plane.h

point.h: ../geom/point.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/point.h point.h

remote_elem.h: ../geom/remote_elem.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/remote_elem.h remote_elem.h

side.h: ../geom/side.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/side.h side.h

sphere.h: ../geom/sphere.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/sphere.h sphere.h

stored_range.h: ../geom/stored_range.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/stored_range.h stored_range.h

surface.h: ../geom/surface.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/geom/surface.h surface.h

abaqus_io.h: ../mesh/abaqus_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/abaqus_io.h abaqus_io.h

boundary_info.h: ../mesh/boundary_info.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/boundary_info.h boundary_info.h

boundary_mesh.h: ../mesh/boundary_mesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/boundary_mesh.h boundary_mesh.h

diva_io.h: ../mesh/diva_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/diva_io.h diva_io.h

ensight_io.h: ../mesh/ensight_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/ensight_io.h ensight_io.h

exodusII_io.h: ../mesh/exodusII_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/exodusII_io.h exodusII_io.h

exodusII_io_helper.h: ../mesh/exodusII_io_helper.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/exodusII_io_helper.h exodusII_io_helper.h

fro_io.h: ../mesh/fro_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/fro_io.h fro_io.h

gmsh_io.h: ../mesh/gmsh_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/gmsh_io.h gmsh_io.h

gmv_io.h: ../mesh/gmv_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/gmv_io.h gmv_io.h

gnuplot_io.h: ../mesh/gnuplot_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/gnuplot_io.h gnuplot_io.h

inf_elem_builder.h: ../mesh/inf_elem_builder.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/inf_elem_builder.h inf_elem_builder.h

legacy_xdr_io.h: ../mesh/legacy_xdr_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/legacy_xdr_io.h legacy_xdr_io.h

matlab_io.h: ../mesh/matlab_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/matlab_io.h matlab_io.h

medit_io.h: ../mesh/medit_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/medit_io.h medit_io.h

mesh_base.h: ../mesh/mesh_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_base.h mesh_base.h

mesh_communication.h: ../mesh/mesh_communication.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_communication.h mesh_communication.h

mesh_data.h: ../mesh/mesh_data.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_data.h mesh_data.h

mesh_function.h: ../mesh/mesh_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_function.h mesh_function.h

mesh_generation.h: ../mesh/mesh_generation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_generation.h mesh_generation.h

mesh.h: ../mesh/mesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh.h mesh.h

mesh_input.h: ../mesh/mesh_input.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_input.h mesh_input.h

mesh_inserter_iterator.h: ../mesh/mesh_inserter_iterator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_inserter_iterator.h mesh_inserter_iterator.h

mesh_modification.h: ../mesh/mesh_modification.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_modification.h mesh_modification.h

mesh_output.h: ../mesh/mesh_output.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_output.h mesh_output.h

mesh_refinement.h: ../mesh/mesh_refinement.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_refinement.h mesh_refinement.h

mesh_serializer.h: ../mesh/mesh_serializer.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_serializer.h mesh_serializer.h

mesh_smoother.h: ../mesh/mesh_smoother.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_smoother.h mesh_smoother.h

mesh_smoother_laplace.h: ../mesh/mesh_smoother_laplace.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_smoother_laplace.h mesh_smoother_laplace.h

mesh_smoother_vsmoother.h: ../mesh/mesh_smoother_vsmoother.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_smoother_vsmoother.h mesh_smoother_vsmoother.h

mesh_tetgen_interface.h: ../mesh/mesh_tetgen_interface.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_tetgen_interface.h mesh_tetgen_interface.h

mesh_tetgen_wrapper.h: ../mesh/mesh_tetgen_wrapper.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_tetgen_wrapper.h mesh_tetgen_wrapper.h

mesh_tools.h: ../mesh/mesh_tools.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_tools.h mesh_tools.h

mesh_triangle_holes.h: ../mesh/mesh_triangle_holes.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_triangle_holes.h mesh_triangle_holes.h

mesh_triangle_interface.h: ../mesh/mesh_triangle_interface.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_triangle_interface.h mesh_triangle_interface.h

mesh_triangle_wrapper.h: ../mesh/mesh_triangle_wrapper.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/mesh_triangle_wrapper.h mesh_triangle_wrapper.h

nemesis_io.h: ../mesh/nemesis_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/nemesis_io.h nemesis_io.h

nemesis_io_helper.h: ../mesh/nemesis_io_helper.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/nemesis_io_helper.h nemesis_io_helper.h

off_io.h: ../mesh/off_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/off_io.h off_io.h

parallel_mesh.h: ../mesh/parallel_mesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/parallel_mesh.h parallel_mesh.h

patch.h: ../mesh/patch.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/patch.h patch.h

postscript_io.h: ../mesh/postscript_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/postscript_io.h postscript_io.h

serial_mesh.h: ../mesh/serial_mesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/serial_mesh.h serial_mesh.h

tecplot_io.h: ../mesh/tecplot_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/tecplot_io.h tecplot_io.h

tetgen_io.h: ../mesh/tetgen_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/tetgen_io.h tetgen_io.h

ucd_io.h: ../mesh/ucd_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/ucd_io.h ucd_io.h

unstructured_mesh.h: ../mesh/unstructured_mesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/unstructured_mesh.h unstructured_mesh.h

unv_io.h: ../mesh/unv_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/unv_io.h unv_io.h

vtk_io.h: ../mesh/vtk_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/vtk_io.h vtk_io.h

xdr_head.h: ../mesh/xdr_head.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_head.h xdr_head.h

xdr_io.h: ../mesh/xdr_io.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_io.h xdr_io.h

xdr_mesh.h: ../mesh/xdr_mesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_mesh.h xdr_mesh.h

xdr_mgf.h: ../mesh/xdr_mgf.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_mgf.h xdr_mgf.h

xdr_mhead.h: ../mesh/xdr_mhead.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_mhead.h xdr_mhead.h

xdr_shead.h: ../mesh/xdr_shead.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_shead.h xdr_shead.h

xdr_soln.h: ../mesh/xdr_soln.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/mesh/xdr_soln.h xdr_soln.h

analytic_function.h: ../numerics/analytic_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/analytic_function.h analytic_function.h

const_function.h: ../numerics/const_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/const_function.h const_function.h

coupling_matrix.h: ../numerics/coupling_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/coupling_matrix.h coupling_matrix.h

dense_matrix_base.h: ../numerics/dense_matrix_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/dense_matrix_base.h dense_matrix_base.h

dense_matrix.h: ../numerics/dense_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/dense_matrix.h dense_matrix.h

dense_submatrix.h: ../numerics/dense_submatrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/dense_submatrix.h dense_submatrix.h

dense_subvector.h: ../numerics/dense_subvector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/dense_subvector.h dense_subvector.h

dense_vector_base.h: ../numerics/dense_vector_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/dense_vector_base.h dense_vector_base.h

dense_vector.h: ../numerics/dense_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/dense_vector.h dense_vector.h

distributed_vector.h: ../numerics/distributed_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/distributed_vector.h distributed_vector.h

function_base.h: ../numerics/function_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/function_base.h function_base.h

laspack_matrix.h: ../numerics/laspack_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/laspack_matrix.h laspack_matrix.h

laspack_vector.h: ../numerics/laspack_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/laspack_vector.h laspack_vector.h

numeric_vector.h: ../numerics/numeric_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/numeric_vector.h numeric_vector.h

parsed_function.h: ../numerics/parsed_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/parsed_function.h parsed_function.h

petsc_macro.h: ../numerics/petsc_macro.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/petsc_macro.h petsc_macro.h

petsc_matrix.h: ../numerics/petsc_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/petsc_matrix.h petsc_matrix.h

petsc_preconditioner.h: ../numerics/petsc_preconditioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/petsc_preconditioner.h petsc_preconditioner.h

petsc_vector.h: ../numerics/petsc_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/petsc_vector.h petsc_vector.h

preconditioner.h: ../numerics/preconditioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/preconditioner.h preconditioner.h

raw_accessor.h: ../numerics/raw_accessor.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/raw_accessor.h raw_accessor.h

refinement_selector.h: ../numerics/refinement_selector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/refinement_selector.h refinement_selector.h

shell_matrix.h: ../numerics/shell_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/shell_matrix.h shell_matrix.h

sparse_matrix.h: ../numerics/sparse_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/sparse_matrix.h sparse_matrix.h

sparse_shell_matrix.h: ../numerics/sparse_shell_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/sparse_shell_matrix.h sparse_shell_matrix.h

sum_shell_matrix.h: ../numerics/sum_shell_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/sum_shell_matrix.h sum_shell_matrix.h

tensor_shell_matrix.h: ../numerics/tensor_shell_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/tensor_shell_matrix.h tensor_shell_matrix.h

tensor_tools.h: ../numerics/tensor_tools.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/tensor_tools.h tensor_tools.h

tensor_value.h: ../numerics/tensor_value.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/tensor_value.h tensor_value.h

trilinos_epetra_matrix.h: ../numerics/trilinos_epetra_matrix.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/trilinos_epetra_matrix.h trilinos_epetra_matrix.h

trilinos_epetra_vector.h: ../numerics/trilinos_epetra_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/trilinos_epetra_vector.h trilinos_epetra_vector.h

trilinos_preconditioner.h: ../numerics/trilinos_preconditioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/trilinos_preconditioner.h trilinos_preconditioner.h

type_n_tensor.h: ../numerics/type_n_tensor.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/type_n_tensor.h type_n_tensor.h

type_tensor.h: ../numerics/type_tensor.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/type_tensor.h type_tensor.h

type_vector.h: ../numerics/type_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/type_vector.h type_vector.h

vector_value.h: ../numerics/vector_value.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/vector_value.h vector_value.h

wrapped_function.h: ../numerics/wrapped_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/wrapped_function.h wrapped_function.h

zero_function.h: ../numerics/zero_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/numerics/zero_function.h zero_function.h

parallel_algebra.h: ../parallel/parallel_algebra.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_algebra.h parallel_algebra.h

parallel_bin_sorter.h: ../parallel/parallel_bin_sorter.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_bin_sorter.h parallel_bin_sorter.h

parallel_conversion_utils.h: ../parallel/parallel_conversion_utils.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_conversion_utils.h parallel_conversion_utils.h

parallel_ghost_sync.h: ../parallel/parallel_ghost_sync.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_ghost_sync.h parallel_ghost_sync.h

parallel.h: ../parallel/parallel.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel.h parallel.h

parallel_hilbert.h: ../parallel/parallel_hilbert.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_hilbert.h parallel_hilbert.h

parallel_histogram.h: ../parallel/parallel_histogram.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_histogram.h parallel_histogram.h

parallel_sort.h: ../parallel/parallel_sort.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/parallel_sort.h parallel_sort.h

threads_allocators.h: ../parallel/threads_allocators.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/threads_allocators.h threads_allocators.h

threads.h: ../parallel/threads.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/parallel/threads.h threads.h

centroid_partitioner.h: ../partitioning/centroid_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/centroid_partitioner.h centroid_partitioner.h

hilbert_sfc_partitioner.h: ../partitioning/hilbert_sfc_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/hilbert_sfc_partitioner.h hilbert_sfc_partitioner.h

linear_partitioner.h: ../partitioning/linear_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/linear_partitioner.h linear_partitioner.h

metis_partitioner.h: ../partitioning/metis_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/metis_partitioner.h metis_partitioner.h

morton_sfc_partitioner.h: ../partitioning/morton_sfc_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/morton_sfc_partitioner.h morton_sfc_partitioner.h

parmetis_partitioner.h: ../partitioning/parmetis_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/parmetis_partitioner.h parmetis_partitioner.h

partitioner.h: ../partitioning/partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/partitioner.h partitioner.h

sfc_partitioner.h: ../partitioning/sfc_partitioner.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/partitioning/sfc_partitioner.h sfc_partitioner.h

diff_physics.h: ../physics/diff_physics.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/physics/diff_physics.h diff_physics.h

diff_qoi.h: ../physics/diff_qoi.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/physics/diff_qoi.h diff_qoi.h

quadrature_clough.h: ../quadrature/quadrature_clough.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_clough.h quadrature_clough.h

quadrature_conical.h: ../quadrature/quadrature_conical.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_conical.h quadrature_conical.h

quadrature_gauss.h: ../quadrature/quadrature_gauss.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_gauss.h quadrature_gauss.h

quadrature_gm.h: ../quadrature/quadrature_gm.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_gm.h quadrature_gm.h

quadrature_grid.h: ../quadrature/quadrature_grid.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_grid.h quadrature_grid.h

quadrature.h: ../quadrature/quadrature.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature.h quadrature.h

quadrature_jacobi.h: ../quadrature/quadrature_jacobi.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_jacobi.h quadrature_jacobi.h

quadrature_monomial.h: ../quadrature/quadrature_monomial.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_monomial.h quadrature_monomial.h

quadrature_rules.h: ../quadrature/quadrature_rules.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_rules.h quadrature_rules.h

quadrature_simpson.h: ../quadrature/quadrature_simpson.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_simpson.h quadrature_simpson.h

quadrature_trap.h: ../quadrature/quadrature_trap.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/quadrature/quadrature_trap.h quadrature_trap.h

derived_rb_construction.h: ../reduced_basis/derived_rb_construction.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/derived_rb_construction.h derived_rb_construction.h

derived_rb_evaluation.h: ../reduced_basis/derived_rb_evaluation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/derived_rb_evaluation.h derived_rb_evaluation.h

rb_assembly_expansion.h: ../reduced_basis/rb_assembly_expansion.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_assembly_expansion.h rb_assembly_expansion.h

rb_construction_base.h: ../reduced_basis/rb_construction_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_construction_base.h rb_construction_base.h

rb_construction.h: ../reduced_basis/rb_construction.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_construction.h rb_construction.h

rb_eim_assembly.h: ../reduced_basis/rb_eim_assembly.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_eim_assembly.h rb_eim_assembly.h

rb_eim_construction.h: ../reduced_basis/rb_eim_construction.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_eim_construction.h rb_eim_construction.h

rb_eim_evaluation.h: ../reduced_basis/rb_eim_evaluation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_eim_evaluation.h rb_eim_evaluation.h

rb_eim_theta.h: ../reduced_basis/rb_eim_theta.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_eim_theta.h rb_eim_theta.h

rb_evaluation.h: ../reduced_basis/rb_evaluation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_evaluation.h rb_evaluation.h

rb_parameters.h: ../reduced_basis/rb_parameters.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_parameters.h rb_parameters.h

rb_parametrized_function.h: ../reduced_basis/rb_parametrized_function.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_parametrized_function.h rb_parametrized_function.h

rb_parametrized.h: ../reduced_basis/rb_parametrized.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_parametrized.h rb_parametrized.h

rb_scm_construction.h: ../reduced_basis/rb_scm_construction.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_scm_construction.h rb_scm_construction.h

rb_scm_evaluation.h: ../reduced_basis/rb_scm_evaluation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_scm_evaluation.h rb_scm_evaluation.h

rb_temporal_discretization.h: ../reduced_basis/rb_temporal_discretization.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_temporal_discretization.h rb_temporal_discretization.h

rb_theta_expansion.h: ../reduced_basis/rb_theta_expansion.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_theta_expansion.h rb_theta_expansion.h

rb_theta.h: ../reduced_basis/rb_theta.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/rb_theta.h rb_theta.h

transient_rb_assembly_expansion.h: ../reduced_basis/transient_rb_assembly_expansion.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/transient_rb_assembly_expansion.h transient_rb_assembly_expansion.h

transient_rb_construction.h: ../reduced_basis/transient_rb_construction.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/transient_rb_construction.h transient_rb_construction.h

transient_rb_evaluation.h: ../reduced_basis/transient_rb_evaluation.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/transient_rb_evaluation.h transient_rb_evaluation.h

transient_rb_theta_expansion.h: ../reduced_basis/transient_rb_theta_expansion.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/reduced_basis/transient_rb_theta_expansion.h transient_rb_theta_expansion.h

adaptive_time_solver.h: ../solvers/adaptive_time_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/adaptive_time_solver.h adaptive_time_solver.h

diff_solver.h: ../solvers/diff_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/diff_solver.h diff_solver.h

eigen_solver.h: ../solvers/eigen_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/eigen_solver.h eigen_solver.h

eigen_time_solver.h: ../solvers/eigen_time_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/eigen_time_solver.h eigen_time_solver.h

euler2_solver.h: ../solvers/euler2_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/euler2_solver.h euler2_solver.h

euler_solver.h: ../solvers/euler_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/euler_solver.h euler_solver.h

laspack_linear_solver.h: ../solvers/laspack_linear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/laspack_linear_solver.h laspack_linear_solver.h

linear.h: ../solvers/linear.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/linear.h linear.h

linear_solver.h: ../solvers/linear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/linear_solver.h linear_solver.h

newton_solver.h: ../solvers/newton_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/newton_solver.h newton_solver.h

nonlinear_solver.h: ../solvers/nonlinear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/nonlinear_solver.h nonlinear_solver.h

petsc_diff_solver.h: ../solvers/petsc_diff_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/petsc_diff_solver.h petsc_diff_solver.h

petscdmlibmesh.h: ../solvers/petscdmlibmesh.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/petscdmlibmesh.h petscdmlibmesh.h

petsc_dm_nonlinear_solver.h: ../solvers/petsc_dm_nonlinear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/petsc_dm_nonlinear_solver.h petsc_dm_nonlinear_solver.h

petsc_linear_solver.h: ../solvers/petsc_linear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/petsc_linear_solver.h petsc_linear_solver.h

petsc_nonlinear_solver.h: ../solvers/petsc_nonlinear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/petsc_nonlinear_solver.h petsc_nonlinear_solver.h

slepc_eigen_solver.h: ../solvers/slepc_eigen_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/slepc_eigen_solver.h slepc_eigen_solver.h

slepc_macro.h: ../solvers/slepc_macro.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/slepc_macro.h slepc_macro.h

solver.h: ../solvers/solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/solver.h solver.h

steady_solver.h: ../solvers/steady_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/steady_solver.h steady_solver.h

time_solver.h: ../solvers/time_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/time_solver.h time_solver.h

trilinos_aztec_linear_solver.h: ../solvers/trilinos_aztec_linear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/trilinos_aztec_linear_solver.h trilinos_aztec_linear_solver.h

trilinos_nox_nonlinear_solver.h: ../solvers/trilinos_nox_nonlinear_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/trilinos_nox_nonlinear_solver.h trilinos_nox_nonlinear_solver.h

twostep_time_solver.h: ../solvers/twostep_time_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/twostep_time_solver.h twostep_time_solver.h

unsteady_solver.h: ../solvers/unsteady_solver.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/solvers/unsteady_solver.h unsteady_solver.h

condensed_eigen_system.h: ../systems/condensed_eigen_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/condensed_eigen_system.h condensed_eigen_system.h

continuation_system.h: ../systems/continuation_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/continuation_system.h continuation_system.h

diff_context.h: ../systems/diff_context.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/diff_context.h diff_context.h

diff_system.h: ../systems/diff_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/diff_system.h diff_system.h

eigen_system.h: ../systems/eigen_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/eigen_system.h eigen_system.h

elem_assembly.h: ../systems/elem_assembly.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/elem_assembly.h elem_assembly.h

equation_systems.h: ../systems/equation_systems.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/equation_systems.h equation_systems.h

explicit_system.h: ../systems/explicit_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/explicit_system.h explicit_system.h

fem_context.h: ../systems/fem_context.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/fem_context.h fem_context.h

fem_system.h: ../systems/fem_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/fem_system.h fem_system.h

frequency_system.h: ../systems/frequency_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/frequency_system.h frequency_system.h

implicit_system.h: ../systems/implicit_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/implicit_system.h implicit_system.h

linear_implicit_system.h: ../systems/linear_implicit_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/linear_implicit_system.h linear_implicit_system.h

newmark_system.h: ../systems/newmark_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/newmark_system.h newmark_system.h

nonlinear_implicit_system.h: ../systems/nonlinear_implicit_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/nonlinear_implicit_system.h nonlinear_implicit_system.h

parameter_vector.h: ../systems/parameter_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/parameter_vector.h parameter_vector.h

qoi_set.h: ../systems/qoi_set.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/qoi_set.h qoi_set.h

sensitivity_data.h: ../systems/sensitivity_data.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/sensitivity_data.h sensitivity_data.h

steady_system.h: ../systems/steady_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/steady_system.h steady_system.h

system.h: ../systems/system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/system.h system.h

system_norm.h: ../systems/system_norm.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/system_norm.h system_norm.h

system_subset_by_subdomain.h: ../systems/system_subset_by_subdomain.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/system_subset_by_subdomain.h system_subset_by_subdomain.h

system_subset.h: ../systems/system_subset.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/system_subset.h system_subset.h

transient_system.h: ../systems/transient_system.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/systems/transient_system.h transient_system.h

compare_types.h: ../utils/compare_types.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/compare_types.h compare_types.h

error_vector.h: ../utils/error_vector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/error_vector.h error_vector.h

hashword.h: ../utils/hashword.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/hashword.h hashword.h

ignore_warnings.h: ../utils/ignore_warnings.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/ignore_warnings.h ignore_warnings.h

location_maps.h: ../utils/location_maps.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/location_maps.h location_maps.h

mapvector.h: ../utils/mapvector.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/mapvector.h mapvector.h

null_output_iterator.h: ../utils/null_output_iterator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/null_output_iterator.h null_output_iterator.h

number_lookups.h: ../utils/number_lookups.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/number_lookups.h number_lookups.h

o_f_stream.h: ../utils/o_f_stream.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/o_f_stream.h o_f_stream.h

ostream_proxy.h: ../utils/ostream_proxy.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/ostream_proxy.h ostream_proxy.h

o_string_stream.h: ../utils/o_string_stream.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/o_string_stream.h o_string_stream.h

parameters.h: ../utils/parameters.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/parameters.h parameters.h

perf_log.h: ../utils/perf_log.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/perf_log.h perf_log.h

perfmon.h: ../utils/perfmon.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/perfmon.h perfmon.h

plt_loader.h: ../utils/plt_loader.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/plt_loader.h plt_loader.h

point_locator_base.h: ../utils/point_locator_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/point_locator_base.h point_locator_base.h

point_locator_list.h: ../utils/point_locator_list.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/point_locator_list.h point_locator_list.h

point_locator_tree.h: ../utils/point_locator_tree.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/point_locator_tree.h point_locator_tree.h

pool_allocator.h: ../utils/pool_allocator.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/pool_allocator.h pool_allocator.h

restore_warnings.h: ../utils/restore_warnings.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/restore_warnings.h restore_warnings.h

statistics.h: ../utils/statistics.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/statistics.h statistics.h

string_to_enum.h: ../utils/string_to_enum.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/string_to_enum.h string_to_enum.h

timestamp.h: ../utils/timestamp.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/timestamp.h timestamp.h

tree_base.h: ../utils/tree_base.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/tree_base.h tree_base.h

tree.h: ../utils/tree.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/tree.h tree.h

tree_node.h: ../utils/tree_node.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/tree_node.h tree_node.h

utility.h: ../utils/utility.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/utility.h utility.h

xdr_cxx.h: ../utils/xdr_cxx.h
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/include/utils/xdr_cxx.h xdr_cxx.h

