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

check_PROGRAMS            = reduced_basis_ex7
check_SCRIPTS             = run.sh
reduced_basis_ex7dir      = $(examples_install_path)/reduced_basis/ex7
reduced_basis_ex7_DATA    = assembly.h horn.geo horn.msh rb_classes.h reduced_basis_ex7.C reduced_basis_ex7.in
reduced_basis_ex7_SOURCES = $(reduced_basis_ex7_DATA) run.sh
reduced_basis_ex7_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES =

clean-local:
	-rm -rf rb_data eim_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex7.in && $(LN_S) $(srcdir)/reduced_basis_ex7.in .
	-rm -f horn.msh             && $(LN_S) $(srcdir)/horn.msh .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex7.in horn.in .linkstamp
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


