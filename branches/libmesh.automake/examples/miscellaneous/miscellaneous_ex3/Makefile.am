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

check_PROGRAMS            = miscellaneous_ex3
miscellaneous_ex3dir      = $(examples_install_path)/miscellaneous/ex3
miscellaneous_ex3_DATA    = miscellaneous_ex3.C ../../adjoints/adjoints_ex1/lshaped.xda
miscellaneous_ex3_SOURCES = $(miscellaneous_ex3_DATA) run.sh
miscellaneous_ex3_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = lshaped.xda out.e 

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


