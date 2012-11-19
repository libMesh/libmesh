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

check_PROGRAMS            = eigenproblems_ex1
eigenproblems_ex1dir      = $(examples_install_path)/eigenproblems/ex1
eigenproblems_ex1_DATA    = eigenproblems_ex1.C
eigenproblems_ex1_SOURCES = $(eigenproblems_ex1_DATA) run.sh
eigenproblems_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.e

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


# Only *run* the test if SLEPc is actually enabled
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
