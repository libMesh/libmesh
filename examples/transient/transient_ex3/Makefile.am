example_name  = transient_ex3
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/transient/ex3
data          = advection_2D.in advection_system.C advection_system.h claw_system.C claw_system.h transient_ex3.C run.sh
sources       = $(data) run.sh

CLEANFILES = claw_solution.*.e

if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f advection_2D.in && $(LN_S) -f $(srcdir)/advection_2D.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += advection_2D.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
