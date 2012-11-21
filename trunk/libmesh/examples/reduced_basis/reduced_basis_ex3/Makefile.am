check_PROGRAMS            = reduced_basis_ex3
reduced_basis_ex3dir      = $(examples_install_path)/reduced_basis/ex3
reduced_basis_ex3_DATA    = reduced_basis_ex3.C assembly.h rb_classes.h reduced_basis_ex3.in run.sh $(top_srcdir)/contrib/utils/Makefile
reduced_basis_ex3_SOURCES = $(reduced_basis_ex3_DATA) run.sh
reduced_basis_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = bf0.e offline_data/* RB_sol.e

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
	-rm -f reduced_basis_ex3.in && $(LN_S) $(srcdir)/reduced_basis_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
