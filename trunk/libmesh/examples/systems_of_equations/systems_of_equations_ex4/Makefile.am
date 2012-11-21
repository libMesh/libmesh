check_PROGRAMS                   = systems_of_equations_ex4
check_SCRIPTS                    = run.sh
systems_of_equations_ex4dir      = $(examples_install_path)/systems_of_equations/ex4
systems_of_equations_ex4_DATA    = systems_of_equations_ex4.C
systems_of_equations_ex4_SOURCES = $(systems_of_equations_ex4_DATA) run.sh
systems_of_equations_ex4_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = displacement.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
