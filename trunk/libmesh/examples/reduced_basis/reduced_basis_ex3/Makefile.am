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

check_PROGRAMS            = reduced_basis_ex3
reduced_basis_ex3dir      = $(examples_install_path)/reduced_basis/ex3
reduced_basis_ex3_DATA    = reduced_basis_ex3.C assembly.h rb_classes.h reduced_basis_ex3.in
reduced_basis_ex3_SOURCES = $(reduced_basis_ex3_DATA) run.sh
reduced_basis_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = bf0.e offline_data/* RB_sol.e

clean-local:
	-rm -rf offline_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex3.in && $(LN_S) $(srcdir)/reduced_basis_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex3.in .linkstamp
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


# Only *run* the tests if SLEPc is actually enabled
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh

# how to make documentation
CLEANFILES += $(check_PROGRAMS).php stdout.log
doc: $(check_PROGRAMS).php
html: $(check_PROGRAMS).php
$(check_PROGRAMS).php: $(check_PROGRAMS) \
                        $(top_srcdir)/contrib/bin/ex2html.sh \
                       Makefile
	$(MAKE) -s run > stdout.log
	$(AM_V_GEN) PATH=$(abs_top_srcdir)/contrib/bin:$(PATH) \
	  ex2html.sh $(check_PROGRAMS) $(abs_srcdir)

endif
