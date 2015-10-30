example_name  = optimization_ex1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/optimization/ex1
data          = optimization_ex1.C optimization_ex1.in run.sh
sources       = $(data) run.sh

CLEANFILES = optimization_soln.exo

clean-local:
	-rm -rf offline_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f optimization_ex1.in && $(LN_S) -f $(srcdir)/optimization_ex1.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += optimization_ex1.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
