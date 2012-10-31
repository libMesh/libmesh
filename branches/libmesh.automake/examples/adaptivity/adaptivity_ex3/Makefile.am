include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_INCLUDE) -I$(top_builddir)/include
LIBS         = $(libmesh_LIBS)

if LIBMESH_ENABLE_FPARSER
  AM_CPPFLAGS += -I$(top_srcdir)/contrib/fparser
endif

check_PROGRAMS         = adaptivity_ex3
check_SCRIPTS          = run.sh
adaptivity_ex3dir      = $(examples_install_path)/adaptivity/ex3
adaptivity_ex3_DATA    = adaptivity_ex3.C adaptivity_ex3.in lshaped.xda lshaped3D.xda
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

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


