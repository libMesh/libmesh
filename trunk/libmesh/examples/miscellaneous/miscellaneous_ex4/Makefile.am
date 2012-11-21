check_PROGRAMS            = miscellaneous_ex4
miscellaneous_ex4dir      = $(examples_install_path)/miscellaneous/ex4
miscellaneous_ex4_DATA    = miscellaneous_ex4.C run.sh $(top_srcdir)/contrib/utils/Makefile
miscellaneous_ex4_SOURCES = $(miscellaneous_ex4_DATA) run.sh
miscellaneous_ex4_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = out_0.vtu out.pvtu

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
