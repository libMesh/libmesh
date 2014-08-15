example_name  = reduced_basis_ex1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/reduced_basis/ex1
data          = reduced_basis_ex1.C assembly.h  rb_classes.h reduced_basis_ex1.in run.sh
sources       = $(data) run.sh

CLEANFILES = RB_sol.e offline_data/RB_F_000.xdr offline_data/RB_A_002.xdr offline_data/RB_A_001.xdr offline_data/RB_A_000.xdr offline_data/parameter_ranges.xdr offline_data/output_003_dual_innerprods.xdr offline_data/output_003_000.xdr offline_data/output_002_dual_innerprods.xdr offline_data/output_002_000.xdr offline_data/output_001_dual_innerprods.xdr offline_data/output_001_000.xdr offline_data/output_000_dual_innerprods.xdr offline_data/output_000_000.xdr offline_data/n_bfs.xdr offline_data/greedy_params.xdr offline_data/Fq_innerprods.xdr offline_data/Fq_Aq_innerprods.xdr offline_data/bf_header.xdr offline_data/bf9.xdr offline_data/bf8.xdr offline_data/bf7.xdr offline_data/bf6.xdr offline_data/bf5.xdr offline_data/bf4.xdr offline_data/bf3.xdr offline_data/bf2.xdr offline_data/bf19.xdr offline_data/bf18.xdr offline_data/bf17.xdr offline_data/bf16.xdr offline_data/bf15.xdr offline_data/bf14.xdr offline_data/bf13.xdr offline_data/bf12.xdr offline_data/bf11.xdr offline_data/bf10.xdr offline_data/bf1.xdr offline_data/bf0.xdr offline_data/Aq_Aq_innerprods.xdr bf0.e

clean-local:
	-rm -rf offline_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex1.in && $(LN_S) -f $(srcdir)/reduced_basis_ex1.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex1.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
