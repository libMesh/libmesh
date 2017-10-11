example_name  = miscellaneous_ex13
install_dir   = $(examples_install_path)/miscellaneous/ex13
data          = miscellaneous_ex13.C run.sh cylinder.exo
sources       = $(data) run.sh
check_SCRIPTS = run.sh

CLEANFILES = out.e cylinder.xdr

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
