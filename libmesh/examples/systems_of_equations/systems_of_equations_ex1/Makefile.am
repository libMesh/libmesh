check_PROGRAMS                   = systems_of_equations_ex1
check_SCRIPTS                    = run.sh
systems_of_equations_ex1dir      = $(examples_install_path)/systems_of_equations/ex1
systems_of_equations_ex1_DATA    = systems_of_equations_ex1.C
systems_of_equations_ex1_SOURCES = $(systems_of_equations_ex1_DATA) run.sh
systems_of_equations_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
