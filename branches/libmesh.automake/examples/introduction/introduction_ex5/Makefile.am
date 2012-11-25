BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS           = introduction_ex5
check_SCRIPTS            = run.sh
introduction_ex5dir      = $(examples_install_path)/introduction/ex5
introduction_ex5_DATA    = introduction_ex5.C exact_solution.C run.sh $(top_srcdir)/contrib/utils/Makefile
introduction_ex5_SOURCES = $(introduction_ex5_DATA) run.sh
introduction_ex5_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES = out_0.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
