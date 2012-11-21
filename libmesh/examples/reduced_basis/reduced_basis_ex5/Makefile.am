check_PROGRAMS            = reduced_basis_ex5
check_SCRIPTS             = run.sh
reduced_basis_ex5dir      = $(examples_install_path)/reduced_basis/ex5
reduced_basis_ex5_DATA    = reduced_basis_ex5.C assembly.h assembly.C rb_classes.h reduced_basis_ex5.in run.sh $(top_srcdir)/contrib/utils/Makefile
reduced_basis_ex5_SOURCES = $(reduced_basis_ex5_DATA) run.sh
reduced_basis_ex5_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = RB_sol.e offline_data/*

clean-local:
	-rm -rf offline_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex5.in && $(LN_S) $(srcdir)/reduced_basis_ex5.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex5.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
