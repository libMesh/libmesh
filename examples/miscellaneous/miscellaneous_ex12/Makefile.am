example_name  = miscellaneous_ex12
install_dir   = $(examples_install_path)/miscellaneous/ex12
data          = miscellaneous_ex12.C run.sh cylinder.exo
sources       = $(data) run.sh
check_SCRIPTS = run.sh

CLEANFILES = out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f cylinder.exo && $(LN_S) -f $(srcdir)/cylinder.exo .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += cylinder.exo .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
