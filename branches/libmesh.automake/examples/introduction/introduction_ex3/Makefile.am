check_PROGRAMS           = introduction_ex3
check_SCRIPTS            = run.sh
introduction_ex3dir      = $(examples_install_path)/introduction/ex3
introduction_ex3_DATA    = introduction_ex3.C exact_solution.C run.sh $(top_srcdir)/contrib/utils/Makefile
introduction_ex3_SOURCES = $(introduction_ex3_DATA) run.sh
introduction_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_0.vtu out.pvtu

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
