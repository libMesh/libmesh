example_name  = fem_system_ex4
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/fem_system/ex4
data          = fem_system_ex4.C fem_system_ex4.in heatsystem.C heatsystem.h run.sh
sources       = $(data) run.sh

CLEANFILES = out_*.e out.gmv

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f fem_system_ex4.in && $(LN_S) -f $(srcdir)/fem_system_ex4.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += fem_system_ex4.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
