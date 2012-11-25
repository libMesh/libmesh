check_PROGRAMS        = vector_fe_ex1
check_SCRIPTS         = run.sh
vector_fe_ex1dir      = $(examples_install_path)/vector_fe/ex1
vector_fe_ex1_DATA    = vector_fe_ex1.C exact_solution.C run.sh $(top_srcdir)/contrib/utils/Makefile
vector_fe_ex1_SOURCES = $(vector_fe_ex1_DATA) run.sh
vector_fe_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.gmv out.e 

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
