example_name = reduced_basis_ex2
install_dir  = $(examples_install_path)/reduced_basis/ex2
data         = reduced_basis_ex2.C assembly.h rb_classes.h reduced_basis_ex2.in run.sh
sources      = $(data) run.sh


CLEANFILES = bf0.e RB_sol.e rb_data/* scm_data/* rb_eval.bin rb_scm_eval.bin

clean-local:
	-rm -rf rb_data scm_data

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
if LIBMESH_ENABLE_GLPK
  check_SCRIPTS = run.sh
endif
endif

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex2.in && $(LN_S) -f $(srcdir)/reduced_basis_ex2.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex2.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
