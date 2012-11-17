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

check_PROGRAMS         = fem_system_ex1
check_SCRIPTS          = run.sh
fem_system_ex1dir      = $(examples_install_path)/fem_system/ex1
fem_system_ex1_DATA    = fem_system_ex1.C fem_system_ex1.in navier.in naviersystem.C naviersystem.h
fem_system_ex1_SOURCES = $(fem_system_ex1_DATA) run.sh
fem_system_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f navier.in && $(LN_S) $(srcdir)/navier.in .
	-rm -f fem_system_ex1.in && $(LN_S) $(srcdir)/fem_system_ex1.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += navier.in fem_system_ex1.in .linkstamp
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


