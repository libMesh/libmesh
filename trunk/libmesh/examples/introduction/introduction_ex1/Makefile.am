check_PROGRAMS           = introduction_ex1
check_SCRIPTS            = run.sh
introduction_ex1dir      = $(examples_install_path)/introduction/ex1
introduction_ex1_DATA    = introduction_ex1.C run.sh $(top_srcdir)/contrib/utils/Makefile
introduction_ex1_SOURCES = $(introduction_ex1_DATA) run.sh
introduction_ex1_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES               = output.xda

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
