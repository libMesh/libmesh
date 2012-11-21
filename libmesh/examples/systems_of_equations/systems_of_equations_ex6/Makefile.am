check_PROGRAMS                   = systems_of_equations_ex6
check_SCRIPTS                    = run.sh
systems_of_equations_ex6dir      = $(examples_install_path)/systems_of_equations/ex6
systems_of_equations_ex6_DATA    = systems_of_equations_ex6.C
systems_of_equations_ex6_SOURCES = $(systems_of_equations_ex6_DATA) run.sh
systems_of_equations_ex6_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = displacement.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
