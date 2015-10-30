example_name  = miscellaneous_ex9
install_dir   = $(examples_install_path)/miscellaneous/ex9
data          = augment_sparsity_on_interface.h miscellaneous_ex9.exo augment_sparsity_on_interface.C miscellaneous_ex9.C run.sh
sources       = $(data) run.sh

check_SCRIPTS = run.sh
CLEANFILES = solution.exo

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f miscellaneous_ex9.exo && $(LN_S) -f $(srcdir)/miscellaneous_ex9.exo
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += miscellaneous_ex9.exo .linkstamp
endif


##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
