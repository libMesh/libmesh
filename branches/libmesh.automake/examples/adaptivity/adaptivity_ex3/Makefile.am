check_PROGRAMS         = adaptivity_ex3
check_SCRIPTS          = run.sh
adaptivity_ex3dir      = $(examples_install_path)/adaptivity/ex3
adaptivity_ex3_DATA    = adaptivity_ex3.C adaptivity_ex3.in lshaped.xda lshaped3D.xda run.sh $(top_srcdir)/contrib/utils/Makefile
adaptivity_ex3_SOURCES = $(adaptivity_ex3_DATA) run.sh
adaptivity_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = lshaped.e biquadratic_h_adaptive.m error_*.e

# need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f lshaped.xda       && $(LN_S) $(srcdir)/lshaped.xda .
	-rm -f lshaped3D.xda     && $(LN_S) $(srcdir)/lshaped3D.xda .
	-rm -f adaptivity_ex3.in && $(LN_S) $(srcdir)/adaptivity_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += lshaped.xda lshaped3D.xda adaptivity_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
