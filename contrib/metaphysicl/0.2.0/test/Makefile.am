check_PROGRAMS  =
check_PROGRAMS += compare_types_unit
check_PROGRAMS += derivs_unit
check_PROGRAMS += divgrad_unit
check_PROGRAMS += dynamic_sparse_vector_navier_unit
check_PROGRAMS += dynamic_sparse_vector_pde_unit
check_PROGRAMS += identities_unit
check_PROGRAMS += instantiations_unit
check_PROGRAMS += main_unit
check_PROGRAMS += shadow_dynamic_sparse_vector_navier_unit
check_PROGRAMS += shadow_dynamic_sparse_vector_pde_unit
check_PROGRAMS += shadow_sparse_struct_navier_unit
check_PROGRAMS += shadow_sparse_struct_pde_unit
check_PROGRAMS += shadow_sparse_vector_navier_unit
check_PROGRAMS += shadow_sparse_vector_pde_unit
check_PROGRAMS += shadow_vector_navier_unit
check_PROGRAMS += shadow_vector_pde_unit
check_PROGRAMS += sparse_derivs_unit
check_PROGRAMS += sparse_identities_unit
check_PROGRAMS += sparse_struct_navier_unit
check_PROGRAMS += sparse_struct_pde_unit
check_PROGRAMS += sparse_vector_navier_unit
check_PROGRAMS += sparse_vector_pde_unit
check_PROGRAMS += testheaders_unit
check_PROGRAMS += testopt_unit
check_PROGRAMS += vector_navier_unit
check_PROGRAMS += vector_pde_unit

if CXX11_ENABLED
  check_PROGRAMS += dualnamedarray_unit
  check_PROGRAMS += namedindexarray_unit
  check_PROGRAMS += physics_unit
endif

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/graphs/include
AM_CPPFLAGS += -I$(top_srcdir)/src/numerics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #metaphysicl_version.h

AM_LDFLAGS = $(top_builddir)/src/libmetaphysicl.la

# Sources for these tests
compare_types_unit_SOURCES = compare_types_unit.C
derivs_unit_SOURCES = derivs_unit.C
divgrad_unit_SOURCES = divgrad_unit.C
dualnamedarray_unit_SOURCES = dualnamedarray_unit.C
dynamic_sparse_vector_navier_unit_SOURCES =  dynamic_sparse_vector_navier_unit.C
dynamic_sparse_vector_navier_unit_SOURCES += navier_unit.h
dynamic_sparse_vector_navier_unit_SOURCES += testing.h
dynamic_sparse_vector_pde_unit_SOURCES =  dynamic_sparse_vector_pde_unit.C
dynamic_sparse_vector_pde_unit_SOURCES += pde_unit.h
dynamic_sparse_vector_pde_unit_SOURCES += testing.h
identities_unit_SOURCES = identities_unit.C
instantiations_unit_SOURCES = instantiations_unit.C
main_unit_SOURCES = main_unit.C
namedindexarray_unit_SOURCES =  namedindexarray_unit.C
shadow_dynamic_sparse_vector_navier_unit_SOURCES =  shadow_dynamic_sparse_vector_navier_unit.C
shadow_dynamic_sparse_vector_navier_unit_SOURCES += navier_unit.h
shadow_dynamic_sparse_vector_navier_unit_SOURCES += testing.h
shadow_dynamic_sparse_vector_pde_unit_SOURCES =  shadow_dynamic_sparse_vector_pde_unit.C
shadow_dynamic_sparse_vector_pde_unit_SOURCES += pde_unit.h
shadow_dynamic_sparse_vector_pde_unit_SOURCES += testing.h
shadow_sparse_struct_navier_unit_SOURCES =  shadow_sparse_struct_navier_unit.C
shadow_sparse_struct_navier_unit_SOURCES += navier_unit.h
shadow_sparse_struct_navier_unit_SOURCES += testing.h
shadow_sparse_struct_pde_unit_SOURCES =  shadow_sparse_struct_pde_unit.C
shadow_sparse_struct_pde_unit_SOURCES += pde_unit.h
shadow_sparse_struct_pde_unit_SOURCES += testing.h
shadow_sparse_vector_navier_unit_SOURCES =  shadow_sparse_vector_navier_unit.C
shadow_sparse_vector_navier_unit_SOURCES += navier_unit.h
shadow_sparse_vector_navier_unit_SOURCES += testing.h
shadow_sparse_vector_pde_unit_SOURCES =  shadow_sparse_vector_pde_unit.C
shadow_sparse_vector_pde_unit_SOURCES += pde_unit.h
shadow_sparse_vector_pde_unit_SOURCES += testing.h
shadow_vector_navier_unit_SOURCES =  shadow_vector_navier_unit.C
shadow_vector_navier_unit_SOURCES += navier_unit.h
shadow_vector_navier_unit_SOURCES += testing.h
shadow_vector_pde_unit_SOURCES =  shadow_vector_pde_unit.C
shadow_vector_pde_unit_SOURCES += pde_unit.h
shadow_vector_pde_unit_SOURCES += testing.h
sparse_derivs_unit_SOURCES = sparse_derivs_unit.C
sparse_identities_unit_SOURCES = sparse_identities_unit.C
sparse_struct_navier_unit_SOURCES =  sparse_struct_navier_unit.C
sparse_struct_navier_unit_SOURCES += navier_unit.h
sparse_struct_navier_unit_SOURCES += testing.h
sparse_struct_pde_unit_SOURCES =  sparse_struct_pde_unit.C
sparse_struct_pde_unit_SOURCES += pde_unit.h
sparse_struct_pde_unit_SOURCES += testing.h
sparse_vector_navier_unit_SOURCES =  sparse_vector_navier_unit.C
sparse_vector_navier_unit_SOURCES += navier_unit.h
sparse_vector_navier_unit_SOURCES += testing.h
sparse_vector_pde_unit_SOURCES =  sparse_vector_pde_unit.C
sparse_vector_pde_unit_SOURCES += pde_unit.h
sparse_vector_pde_unit_SOURCES += testing.h
physics_unit_SOURCES = physics_unit.C
testheaders_unit_SOURCES = testheaders_unit.C
testopt_unit_SOURCES = testopt_unit.C
vector_navier_unit_SOURCES =  vector_navier_unit.C
vector_navier_unit_SOURCES += navier_unit.h
vector_navier_unit_SOURCES += testing.h
vector_pde_unit_SOURCES =  vector_pde_unit.C
vector_pde_unit_SOURCES += pde_unit.h
vector_pde_unit_SOURCES += testing.h

#Define tests to actually be run
TESTS  =
TESTS += compare_types_unit
TESTS += derivs_unit
TESTS += divgrad_unit
#TESTS += dynamic_sparse_vector_navier_unit
TESTS += dynamic_sparse_vector_pde_unit
TESTS += identities_unit
TESTS += instantiations_unit
TESTS += main_unit
#TESTS += shadow_dynamic_sparse_vector_navier_unit
TESTS += shadow_dynamic_sparse_vector_pde_unit
TESTS += shadow_sparse_struct_navier_unit
TESTS += shadow_sparse_struct_pde_unit
TESTS += shadow_sparse_vector_navier_unit
TESTS += shadow_sparse_vector_pde_unit
TESTS += shadow_vector_navier_unit
TESTS += shadow_vector_pde_unit
TESTS += sparse_derivs_unit
TESTS += sparse_identities_unit
TESTS += sparse_struct_navier_unit
TESTS += sparse_struct_pde_unit
TESTS += sparse_vector_navier_unit
TESTS += sparse_vector_pde_unit
TESTS += testheaders_unit
TESTS += testopt_unit
TESTS += vector_navier_unit
TESTS += vector_pde_unit

if CXX11_ENABLED
  TESTS += dualnamedarray_unit
  TESTS += namedindexarray_unit
  TESTS += physics_unit
endif

#-------------
# MASA support
#-------------
if MASA_ENABLED
  AM_CPPFLAGS += $(MASA_CXXFLAGS)
  AM_LDFLAGS  += $(MASA_LIBS)
endif

if VEXCL_ENABLED
  AM_CPPFLAGS += $(VEXCL_CPPFLAGS)
  AM_LDFLAGS  += $(VEXCL_LDFLAGS)
  LIBS        += $(VEXCL_LIBS)
endif

CLEANFILES =

if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

if METAPHYSICL_LICENSESTAMPEXISTS
BUILT_SOURCES   = .license.stamp

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

CLEANFILES += .license.stamp
endif

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
