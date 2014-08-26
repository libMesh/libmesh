example_name  = reduced_basis_ex6
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/reduced_basis/ex6
data          = reduced_basis_ex6.C assembly.h eim_classes.h rb_classes.h reduced_basis_ex6.in eim.in rb.in run.sh
sources       = $(data) run.sh

CLEANFILES = RB_sol.e rb_data/RB_F_001.xdr rb_data/RB_F_000.xdr rb_data/RB_A_005.xdr rb_data/RB_A_004.xdr rb_data/RB_A_003.xdr rb_data/RB_A_002.xdr rb_data/RB_A_001.xdr rb_data/RB_A_000.xdr rb_data/parameter_ranges.xdr rb_data/n_bfs.xdr rb_data/greedy_params.xdr rb_data/Fq_innerprods.xdr rb_data/Fq_Aq_innerprods.xdr rb_data/bf_header.xdr rb_data/bf4.xdr rb_data/bf3.xdr rb_data/bf2.xdr rb_data/bf1.xdr rb_data/bf0.xdr rb_data/Aq_Aq_innerprods.xdr eim_data/RB_inner_product_matrix.xdr eim_data/parameter_ranges.xdr eim_data/n_bfs.xdr eim_data/interpolation_points_var.xdr eim_data/interpolation_points.xdr eim_data/interpolation_matrix.xdr eim_data/greedy_params.xdr eim_data/Fq_innerprods.xdr eim_data/Fq_Aq_innerprods.xdr eim_data/extra_interpolation_point_var.xdr eim_data/extra_interpolation_point.xdr eim_data/extra_interpolation_matrix_row.xdr eim_data/bf_header.xdr eim_data/bf2.xdr eim_data/bf1.xdr eim_data/bf0.xdr eim_data/Aq_Aq_innerprods.xdr

clean-local:
	-rm -rf rb_data eim_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex6.in && $(LN_S) -f $(srcdir)/reduced_basis_ex6.in .
	-rm -f eim.in               && $(LN_S) -f $(srcdir)/eim.in .
	-rm -f rb.in                && $(LN_S) -f $(srcdir)/rb.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex6.in eim.in rb.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
