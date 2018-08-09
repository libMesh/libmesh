example_name  = adaptivity_ex3
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adaptivity/ex3
data          = adaptivity_ex3.C adaptivity_ex3.in lshaped.xda lshaped3D.xda run.sh
sources       = $(data) run.sh

CLEANFILES = lshaped.e *.m error_*.e error_*.n*

# need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f lshaped.xda       && $(LN_S) -f $(srcdir)/lshaped.xda .
	-rm -f lshaped3D.xda     && $(LN_S) -f $(srcdir)/lshaped3D.xda .
	-rm -f adaptivity_ex3.in && $(LN_S) -f $(srcdir)/adaptivity_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += lshaped.xda lshaped3D.xda adaptivity_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
