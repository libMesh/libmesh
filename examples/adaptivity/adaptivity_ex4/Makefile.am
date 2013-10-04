example_name  = adaptivity_ex4
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adaptivity/ex4
data          = adaptivity_ex4.C adaptivity_ex4.in run.sh
sources       = $(data) run.sh

CLEANFILES = 2D_hermite_adaptive.m 2D_hermite_adaptive.e

# need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f adaptivity_ex4.in && $(LN_S) $(srcdir)/adaptivity_ex4.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += adaptivity_ex4.in .linkstamp
endif


##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
