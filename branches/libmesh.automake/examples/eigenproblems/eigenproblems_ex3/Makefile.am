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

check_PROGRAMS            = eigenproblems_ex3
eigenproblems_ex3dir      = $(examples_install_path)/eigenproblems/ex3
eigenproblems_ex3_DATA    = eigenproblems_ex3.C drum1_mesh.e drum2_mesh.e
eigenproblems_ex3_SOURCES = $(eigenproblems_ex3_DATA) run.sh
eigenproblems_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES                = drum1_evals.txt drum2_evals.txt drum1_evec.e drum2_evec.e

# Only *run* the tests if SLEPc is actually enabled
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f drum1_mesh.e && $(LN_S) $(srcdir)/drum1_mesh.e .
	-rm -f drum2_mesh.e && $(LN_S) $(srcdir)/drum2_mesh.e .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += drum*_mesh.xda .linkstamp
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


