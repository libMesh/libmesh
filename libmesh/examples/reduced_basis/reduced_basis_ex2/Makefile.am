check_PROGRAMS            = reduced_basis_ex2
reduced_basis_ex2dir      = $(examples_install_path)/reduced_basis/ex2
reduced_basis_ex2_DATA    = reduced_basis_ex2.C assembly.h rb_classes.h reduced_basis_ex2.in
reduced_basis_ex2_SOURCES = $(reduced_basis_ex2_DATA) run.sh
reduced_basis_ex2_LDADD   = $(top_builddir)/libmesh.la


CLEANFILES = bf0.e RB_sol.e rb_data/* scm_data/*

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
	-rm -f reduced_basis_ex2.in && $(LN_S) $(srcdir)/reduced_basis_ex2.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex2.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
