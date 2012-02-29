#SUBDIRS = introduction

include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)
INCLUDES    = $(libmesh_INCLUDE)
LIBS        = $(libmesh_LIBS)

EXTRA_DIST = run_common.sh

check_PROGRAMS =
check_SCRIPTS = 
CLEANFILES = biquadratic_h_adaptive.m *.gmv* *.xda *.e *.e.* *.vtu *.pvtu adaptivity_ex3.in  general.in  l-shaped.in navier.in fem_system_ex1.in #FIXME - tighten up!


examples_source_path=$(top_srcdir)/examples
examples_install_path=$(prefix)/examples


#
# Introductory examples
#
check_PROGRAMS += introduction_ex1
check_SCRIPTS  += introduction/introduction_ex1/run.sh
introduction_ex1dir      = $(examples_install_path)/introduction/ex1
introduction_ex1_DATA    = introduction/introduction_ex1/introduction_ex1.C
introduction_ex1_SOURCES = $(introduction_ex1_DATA) introduction/introduction_ex1/run.sh
introduction_ex1_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES += output.xda

check_PROGRAMS += introduction_ex2
check_SCRIPTS  += introduction/introduction_ex2/run.sh
introduction_ex2dir      = $(examples_install_path)/introduction/ex2
introduction_ex2_DATA    = introduction/introduction_ex2/introduction_ex2.C
introduction_ex2_SOURCES = $(introduction_ex2_DATA) introduction/introduction_ex2/run.sh
introduction_ex2_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES += eqn_sys.dat

check_PROGRAMS += introduction_ex3
check_SCRIPTS  += introduction/introduction_ex3/run.sh
introduction_ex3dir      = $(examples_install_path)/introduction/ex3
introduction_ex3_DATA    = introduction/introduction_ex3/introduction_ex3.C introduction/introduction_ex3/exact_solution.C 
introduction_ex3_SOURCES = $(introduction_ex3_DATA) introduction/introduction_ex3/run.sh
introduction_ex3_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += introduction_ex4
check_SCRIPTS  += introduction/introduction_ex4/run.sh
introduction_ex4dir      = $(examples_install_path)/introduction/ex4
introduction_ex4_DATA    = introduction/introduction_ex4/introduction_ex4.C introduction/introduction_ex3/exact_solution.C 
introduction_ex4_SOURCES = $(introduction_ex4_DATA) introduction/introduction_ex4/run.sh
introduction_ex4_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += introduction_ex5
check_SCRIPTS  += introduction/introduction_ex5/run.sh
introduction_ex5dir      = $(examples_install_path)/introduction/ex5
introduction_ex5_DATA    = introduction/introduction_ex5/introduction_ex5.C introduction/introduction_ex3/exact_solution.C
introduction_ex5_SOURCES = $(introduction_ex5_DATA) introduction/introduction_ex5/run.sh
introduction_ex5_LDADD   = $(top_builddir)/libmesh.la
#CLEANFILES += 

CLEANFILES += out_1 out_1_data gnupl* *.xdr 2D_hermite_adaptive.m # FIXME - make this stricter

#
# Adaptive Mesh Refinement
# 
check_PROGRAMS += adaptivity_ex1
check_SCRIPTS  += adaptivity/adaptivity_ex1/run.sh
adaptivity_ex1dir      = $(examples_install_path)/adaptivity/ex1
adaptivity_ex1_DATA    = adaptivity/adaptivity_ex1/adaptivity_ex1.C
adaptivity_ex1_SOURCES = $(adaptivity_ex1_DATA) adaptivity/adaptivity_ex1/run.sh
adaptivity_ex1_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += adaptivity_ex2
check_SCRIPTS  += adaptivity/adaptivity_ex2/run.sh
adaptivity_ex2dir      = $(examples_install_path)/adaptivity/ex2
adaptivity_ex2_DATA    = adaptivity/adaptivity_ex2/adaptivity_ex2.C transient/transient_ex1/exact_solution.C adaptivity/adaptivity_ex2/mesh.xda
adaptivity_ex2_SOURCES = $(adaptivity_ex2_DATA) adaptivity/adaptivity_ex2/run.sh
adaptivity_ex2_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += adaptivity_ex3
check_SCRIPTS  += adaptivity/adaptivity_ex3/run.sh
adaptivity_ex3dir      = $(examples_install_path)/adaptivity/ex3
adaptivity_ex3_DATA    = adaptivity/adaptivity_ex3/adaptivity_ex3.C adaptivity/adaptivity_ex3/adaptivity_ex3.in adaptivity/adaptivity_ex3/lshaped.xda adaptivity/adaptivity_ex3/lshaped3D.xda
adaptivity_ex3_SOURCES = $(adaptivity_ex3_DATA) adaptivity/adaptivity_ex3/run.sh
adaptivity_ex3_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += adaptivity_ex4
check_SCRIPTS  += adaptivity/adaptivity_ex4/run.sh
adaptivity_ex4dir      = $(examples_install_path)/adaptivity/ex4
adaptivity_ex4_DATA    = adaptivity/adaptivity_ex4/adaptivity_ex4.C adaptivity/adaptivity_ex4/adaptivity_ex4.in
adaptivity_ex4_SOURCES = $(adaptivity_ex4_DATA) adaptivity/adaptivity_ex4/run.sh
adaptivity_ex4_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES += biquadratic_h_adaptive.m

check_PROGRAMS += adaptivity_ex5
check_SCRIPTS  += adaptivity/adaptivity_ex5/run.sh
adaptivity_ex5dir      = $(examples_install_path)/adaptivity/ex5
adaptivity_ex5_DATA    = adaptivity/adaptivity_ex5/adaptivity_ex5.C transient/transient_ex1/exact_solution.C
adaptivity_ex5_SOURCES = $(adaptivity_ex5_DATA) adaptivity/adaptivity_ex5/run.sh
adaptivity_ex5_LDADD   = $(top_builddir)/libmesh.la


#
# Adjoints
# 
check_PROGRAMS += adjoints_ex1
check_SCRIPTS  += adjoints/adjoints_ex1/run.sh
adjoints_ex1dir      = $(examples_install_path)/adjoints/ex1
adjoints_ex1_DATA    = adjoints/adjoints_ex1/adjoints_ex1.C adjoints/adjoints_ex1/L-shaped.C adjoints/adjoints_ex1/L-shaped.h adjoints/adjoints_ex1/element_postprocess.C adjoints/adjoints_ex1/element_qoi_derivative.C adjoints/adjoints_ex1/femparameters.C adjoints/adjoints_ex1/femparameters.h adjoints/adjoints_ex1/general.in adjoints/adjoints_ex1/l-shaped.in adjoints/adjoints_ex1/lshaped.xda adjoints/adjoints_ex1/side_postprocess.C adjoints/adjoints_ex1/side_qoi_derivative.C
adjoints_ex1_SOURCES = $(adjoints_ex1_DATA) adjoints/adjoints_ex1/run.sh
adjoints_ex1_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += adjoints_ex2
check_SCRIPTS  += adjoints/adjoints_ex2/run.sh
adjoints_ex2dir      = $(examples_install_path)/adjoints/ex2
adjoints_ex2_DATA    = adjoints/adjoints_ex2/adjoints_ex2.C adjoints/adjoints_ex1/element_qoi_derivative.C adjoints/adjoints_ex1/lshaped.xda adjoints/adjoints_ex2/L-shaped.C adjoints/adjoints_ex2/L-shaped.h adjoints/adjoints_ex2/femparameters.C adjoints/adjoints_ex2/femparameters.h adjoints/adjoints_ex2/general.in adjoints/adjoints_ex2/l-shaped.in adjoints/adjoints_ex2/element_qoi.C
adjoints_ex2_SOURCES = $(adjoints_ex2_DATA) adjoints/adjoints_ex2/run.sh
adjoints_ex2_LDADD   = $(top_builddir)/libmesh.la


#
# Eigenproblems
# 
check_PROGRAMS += eigenproblems_ex1
eigenproblems_ex1dir      = $(examples_install_path)/eigenproblems/ex1
eigenproblems_ex1_DATA    = eigenproblems/eigenproblems_ex1/eigenproblems_ex1.C
eigenproblems_ex1_SOURCES = $(eigenproblems_ex1_DATA) eigenproblems/eigenproblems_ex1/run.sh
eigenproblems_ex1_LDADD   = $(top_builddir)/libmesh.la

check_PROGRAMS += eigenproblems_ex2
eigenproblems_ex2dir      = $(examples_install_path)/eigenproblems/ex2
eigenproblems_ex2_DATA    = eigenproblems/eigenproblems_ex2/eigenproblems_ex2.C
eigenproblems_ex2_SOURCES = $(eigenproblems_ex2_DATA) eigenproblems/eigenproblems_ex2/run.sh
eigenproblems_ex2_LDADD   = $(top_builddir)/libmesh.la

if LIBMESH_ENABLE_SLEPC #FIXME: Remove this when example is working again.
  check_PROGRAMS += eigenproblems_ex3
  eigenproblems_ex3dir      = $(examples_install_path)/eigenproblems/ex3
  eigenproblems_ex3_DATA    = eigenproblems/eigenproblems_ex3/eigenproblems_ex3.C eigenproblems/eigenproblems_ex3/drum1_mesh.e  eigenproblems/eigenproblems_ex3/drum2_mesh.e
  eigenproblems_ex3_SOURCES = $(eigenproblems_ex3_DATA) eigenproblems/eigenproblems_ex3/run.sh
  eigenproblems_ex3_LDADD   = $(top_builddir)/libmesh.la
endif

# Only *run* the tests if SLEPc is actually enabled
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS  += eigenproblems/eigenproblems_ex1/run.sh
  check_SCRIPTS  += eigenproblems/eigenproblems_ex2/run.sh
  check_SCRIPTS  += eigenproblems/eigenproblems_ex3/run.sh
endif

#
# FEMSystem
#
check_PROGRAMS += fem_system_ex1
check_SCRIPTS  += fem_system/fem_system_ex1/run.sh
fem_system_ex1dir      = $(examples_install_path)/fem_system/ex1
fem_system_ex1_DATA    = fem_system/fem_system_ex1/fem_system_ex1.C fem_system/fem_system_ex1/fem_system_ex1.in  fem_system/fem_system_ex1/navier.in  fem_system/fem_system_ex1/naviersystem.C  fem_system/fem_system_ex1/naviersystem.h
fem_system_ex1_SOURCES = $(fem_system_ex1_DATA) fem_system/fem_system_ex1/run.sh
fem_system_ex1_LDADD   = $(top_builddir)/libmesh.la


#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
