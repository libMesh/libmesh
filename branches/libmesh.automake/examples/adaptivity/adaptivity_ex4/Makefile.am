check_PROGRAMS         = adaptivity_ex4
check_SCRIPTS          = run.sh
adaptivity_ex4dir      = $(examples_install_path)/adaptivity/ex4
adaptivity_ex4_DATA    = adaptivity_ex4.C adaptivity_ex4.in run.sh $(top_srcdir)/contrib/utils/Makefile
adaptivity_ex4_SOURCES = $(adaptivity_ex4_DATA) run.sh
adaptivity_ex4_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = 2D_hermite_adaptive.m 2D_hermite_adaptive.e 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
