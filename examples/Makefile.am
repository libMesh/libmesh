SUBDIRS = \
	introduction/introduction_ex1 \
	introduction/introduction_ex2 \
	introduction/introduction_ex3 \
	introduction/introduction_ex4 \
	introduction/introduction_ex5 \
	adaptivity/adaptivity_ex1 \
	adaptivity/adaptivity_ex2 \
	adaptivity/adaptivity_ex3 \
	adaptivity/adaptivity_ex4 \
	adaptivity/adaptivity_ex5 \
	adjoints/adjoints_ex1 \
	adjoints/adjoints_ex2 \
	adjoints/adjoints_ex3 \
	adjoints/adjoints_ex4 \
	adjoints/adjoints_ex5 \
	eigenproblems/eigenproblems_ex1 \
	eigenproblems/eigenproblems_ex2 \
	eigenproblems/eigenproblems_ex3 \
	fem_system/fem_system_ex1 \
	fem_system/fem_system_ex2 \
	solution_transfer/solution_transfer_ex1 \
	miscellaneous/miscellaneous_ex1 \
	miscellaneous/miscellaneous_ex2 \
	miscellaneous/miscellaneous_ex3 \
	miscellaneous/miscellaneous_ex4 \
	miscellaneous/miscellaneous_ex5 \
	miscellaneous/miscellaneous_ex6 \
	miscellaneous/miscellaneous_ex7 \
	miscellaneous/miscellaneous_ex8 \
	miscellaneous/miscellaneous_ex9 \
	subdomains/subdomains_ex1 \
	subdomains/subdomains_ex2 \
	systems_of_equations/systems_of_equations_ex1 \
	systems_of_equations/systems_of_equations_ex2 \
	systems_of_equations/systems_of_equations_ex3 \
	systems_of_equations/systems_of_equations_ex4 \
	systems_of_equations/systems_of_equations_ex5 \
	systems_of_equations/systems_of_equations_ex6 \
	reduced_basis/reduced_basis_ex1 \
	reduced_basis/reduced_basis_ex2 \
	reduced_basis/reduced_basis_ex3 \
	reduced_basis/reduced_basis_ex4 \
	reduced_basis/reduced_basis_ex5 \
	reduced_basis/reduced_basis_ex6 \
	reduced_basis/reduced_basis_ex7 \
	transient/transient_ex1 \
	transient/transient_ex2 \
	vector_fe/vector_fe_ex1 \
	vector_fe/vector_fe_ex2 \
	vector_fe/vector_fe_ex3 \
        vector_fe/vector_fe_ex4

AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_optional_INCLUDES) -I$(top_builddir)/include
AM_LDFLAGS   = $(libmesh_LDFLAGS)
LIBS         = $(libmesh_optional_LIBS)

examples_source_path  = $(top_srcdir)/examples
examples_install_path = $(prefix)/examples

#
# anything defined in 'data_DATA' will be installed in 'datadir'
datadir   = $(examples_install_path)
data_DATA = run_common.sh

EXTRA_DIST = $(data_DATA)

# include our header checking script when doing 'make dist'
EXTRA_DIST += ../contrib/bin/test_installed_examples.sh Make.common


######################################################################
#
# What to do for 'make installcheck'
#bindir=@bindir@/@host@_@METHOD@
#libdir=@libdir@/@host@_@METHOD@
pkgconfigdir=$(libdir)/pkgconfig

installcheck-local:
	@echo " "
	@echo "Checking for self-sufficient examples..."
	@echo " "
	@CXX="$(CXX)" \
	  examples_install_path=$(examples_install_path) \
	  PKG_CONFIG=$(PKG_CONFIG) \
	  PKG_CONFIG_PATH=$(pkgconfigdir) \
	  LIBMESH_CONFIG_PATH=$(bindir) \
	  $(top_srcdir)/contrib/bin/test_installed_examples.sh

######################################################################
#
# What to do for 'make doc'
doc:
	for dir in $(SUBDIRS) ; do \
	  cd $$dir ; \
	  $(MAKE) doc ; \
	  cp *.php $(abs_top_srcdir)/doc/html ; \
	  cd - ; \
	done

######################################################################
#
# Don't leave code coverage outputs lying around
if CODE_COVERAGE_ENABLED
  CLEANFILES = */*/*.gcda */*/*.gcno
endif
