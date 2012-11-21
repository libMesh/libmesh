check_PROGRAMS           = introduction_ex2
check_SCRIPTS            = run.sh
introduction_ex2dir      = $(examples_install_path)/introduction/ex2
introduction_ex2_DATA    = introduction_ex2.C
introduction_ex2_SOURCES = $(introduction_ex2_DATA) run.sh
introduction_ex2_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES               = eqn_sys.dat

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
