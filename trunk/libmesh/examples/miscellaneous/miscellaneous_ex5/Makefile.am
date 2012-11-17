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

# need to link lshaped.xda
BUILT_SOURCES = lshaped.xda
lshaped.xda:
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/examples/adjoints/adjoints_ex1/lshaped.xda .

check_PROGRAMS            = miscellaneous_ex5
miscellaneous_ex5dir      = $(examples_install_path)/miscellaneous/ex5
miscellaneous_ex5_DATA    = miscellaneous_ex5.C ../../adjoints/adjoints_ex1/lshaped.xda lshaped3D.xda miscellaneous_ex5.in
miscellaneous_ex5_SOURCES = $(miscellaneous_ex5_DATA) run.sh
miscellaneous_ex5_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = lshaped.xda lshaped_dg.e 

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f miscellaneous_ex5.in && $(LN_S) $(srcdir)/miscellaneous_ex5.in .
	-rm -f lshaped3D.xda        && $(LN_S) $(srcdir)/lshaped3D.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += miscellaneous_ex5.in lshaped3D.xda .linkstamp
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


