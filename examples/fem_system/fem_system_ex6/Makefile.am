example_name  = fem_system_ex6
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/fem_system/ex6
data          = fem_system_ex6.C fem_system_ex6.in run.sh
sources       = $(data)

CLEANFILES =

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f fem_system_ex6.in && $(LN_S) -f $(srcdir)/fem_system_ex6.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += fem_system_ex6.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
