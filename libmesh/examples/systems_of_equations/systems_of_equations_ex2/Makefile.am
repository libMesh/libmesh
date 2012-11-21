check_PROGRAMS                   = systems_of_equations_ex2
check_SCRIPTS                    = run.sh
systems_of_equations_ex2dir      = $(examples_install_path)/systems_of_equations/ex2
systems_of_equations_ex2_DATA    = systems_of_equations_ex2.C gmv_attributes.attr movie.sh
systems_of_equations_ex2_SOURCES = $(systems_of_equations_ex2_DATA) run.sh
systems_of_equations_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
