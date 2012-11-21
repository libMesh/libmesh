check_PROGRAMS         = adaptivity_ex1
check_SCRIPTS          = run.sh
adaptivity_ex1dir      = $(examples_install_path)/adaptivity/ex1
adaptivity_ex1_DATA    = adaptivity_ex1.C
adaptivity_ex1_SOURCES = $(adaptivity_ex1_DATA) run.sh
adaptivity_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = gnuplot_script_data gnuplot_script 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
