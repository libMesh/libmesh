example_name  = fem_system_ex3
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/fem_system/ex3
data          = fem_system_ex3.C fem_system_ex3.in elasticity_system.C elasticity_system.h run.sh
sources       = $(data) run.sh

CLEANFILES = out.e* out.newmark.* out.steady.*

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f fem_system_ex3.in && $(LN_S) -f $(srcdir)/fem_system_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += fem_system_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
