check_PROGRAMS            = miscellaneous_ex7
miscellaneous_ex7dir      = $(examples_install_path)/miscellaneous/ex7
miscellaneous_ex7_DATA    = biharmonic.C biharmonic.h biharmonic_jr.C biharmonic_jr.h miscellaneous_ex7.C
miscellaneous_ex7_SOURCES = $(miscellaneous_ex7_DATA) run.sh
miscellaneous_ex7_LDADD   = $(top_builddir)/libmesh.la
check_SCRIPTS             = run.sh

CLEANFILES = bih.1.e 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
