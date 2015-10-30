example_name  = fem_system_ex2
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/fem_system/ex2
data          = fem_system_ex2.C solid.in solid_system.C solid_system.h nonlinear_neohooke_cc.C nonlinear_neohooke_cc.h run.sh
sources       = $(data) run.sh

CLEANFILES = *vtu

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f solid.in && $(LN_S) -f $(srcdir)/solid.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += solid.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
