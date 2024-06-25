example_name  = miscellaneous_ex16
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/miscellaneous/ex16
data          = miscellaneous_ex16.C exact_solution.C miscellaneous_ex16.in run.sh
sources       = $(data)

CLEANFILES = out.pvtu

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f miscellaneous_ex16.in && $(LN_S) -f $(srcdir)/miscellaneous_ex16.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += miscellaneous_ex16.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
