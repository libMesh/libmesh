example_name  = fem_system_ex1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/fem_system/ex1
data          = fem_system_ex1.C fem_system_ex1.in navier.in naviersystem.C naviersystem.h run.sh
sources       = $(data) run.sh

CLEANFILES = out_*.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f navier.in && $(LN_S) -f $(srcdir)/navier.in .
	-rm -f fem_system_ex1.in && $(LN_S) -f $(srcdir)/fem_system_ex1.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += navier.in fem_system_ex1.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
