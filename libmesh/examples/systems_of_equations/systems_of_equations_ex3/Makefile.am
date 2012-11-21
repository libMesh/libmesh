check_PROGRAMS                   = systems_of_equations_ex3
check_SCRIPTS                    = run.sh
systems_of_equations_ex3dir      = $(examples_install_path)/systems_of_equations/ex3
systems_of_equations_ex3_DATA    = systems_of_equations_ex3.C run.sh $(top_srcdir)/contrib/utils/Makefile
systems_of_equations_ex3_SOURCES = $(systems_of_equations_ex3_DATA) run.sh
systems_of_equations_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
