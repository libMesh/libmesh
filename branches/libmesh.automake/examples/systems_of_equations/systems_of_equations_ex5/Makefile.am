check_PROGRAMS                   = systems_of_equations_ex5
check_SCRIPTS                    = run.sh
systems_of_equations_ex5dir      = $(examples_install_path)/systems_of_equations/ex5
systems_of_equations_ex5_DATA    = systems_of_equations_ex5.C run.sh $(top_srcdir)/contrib/utils/Makefile
systems_of_equations_ex5_SOURCES = $(systems_of_equations_ex5_DATA) run.sh
systems_of_equations_ex5_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = displacement.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
