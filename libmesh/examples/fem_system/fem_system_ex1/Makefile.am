check_PROGRAMS         = fem_system_ex1
check_SCRIPTS          = run.sh
fem_system_ex1dir      = $(examples_install_path)/fem_system/ex1
fem_system_ex1_DATA    = fem_system_ex1.C fem_system_ex1.in navier.in naviersystem.C naviersystem.h run.sh $(top_srcdir)/contrib/utils/Makefile
fem_system_ex1_SOURCES = $(fem_system_ex1_DATA) run.sh
fem_system_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f navier.in && $(LN_S) $(srcdir)/navier.in .
	-rm -f fem_system_ex1.in && $(LN_S) $(srcdir)/fem_system_ex1.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += navier.in fem_system_ex1.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
