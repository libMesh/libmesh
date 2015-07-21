example_name = reduced_basis_ex3
install_dir  = $(examples_install_path)/reduced_basis/ex3
data         = reduced_basis_ex3.C assembly.h rb_classes.h reduced_basis_ex3.in run.sh
sources      = $(data) run.sh

CLEANFILES = bf0.e offline_data/* RB_sol.e trans_rb_eval.bin

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

clean-local:
	-rm -rf offline_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex3.in && $(LN_S) -f $(srcdir)/reduced_basis_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
