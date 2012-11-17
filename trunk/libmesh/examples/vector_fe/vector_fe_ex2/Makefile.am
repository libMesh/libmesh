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

check_PROGRAMS        = vector_fe_ex2
check_SCRIPTS         = run.sh
vector_fe_ex2dir      = $(examples_install_path)/vector_fe/ex2
vector_fe_ex2_DATA    = laplace.in laplace_system.C solution_function.h vector_fe_ex2.in laplace_exact_solution.h laplace_system.h vector_fe_ex2.C

vector_fe_ex2_SOURCES = $(vector_fe_ex2_DATA) run.sh
vector_fe_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.e out.inp

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex2.in && $(LN_S) $(srcdir)/vector_fe_ex2.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex2.in .linkstamp
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


