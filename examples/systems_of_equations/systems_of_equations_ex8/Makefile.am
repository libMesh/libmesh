example_name  = systems_of_equations_ex8
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/systems_of_equations/ex8
data          = systems_of_equations_ex8.C run.sh systems_of_equations_ex8.in systems_of_equations_ex8.exo linear_elasticity_with_contact.h linear_elasticity_with_contact.C
sources       = $(data)

CLEANFILES = solution.exo
# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f systems_of_equations_ex8.in && $(LN_S) -f $(srcdir)/systems_of_equations_ex8.in .
	-rm -f systems_of_equations_ex8.exo && $(LN_S) -f $(srcdir)/systems_of_equations_ex8.exo .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += systems_of_equations_ex8.in systems_of_equations_ex8.exo .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
