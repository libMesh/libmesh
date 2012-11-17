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

check_PROGRAMS        = transient_ex2
check_SCRIPTS         = run.sh
transient_ex2dir      = $(examples_install_path)/transient/ex2
transient_ex2_DATA    = transient_ex2.C pipe-mesh.unv
transient_ex2_SOURCES = $(transient_ex2_DATA) run.sh
transient_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = pressure_node.res out.*.gmv

if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f pipe-mesh.unv && $(LN_S) $(srcdir)/pipe-mesh.unv .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += pipe-mesh.unv .linkstamp
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


