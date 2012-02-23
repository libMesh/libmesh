#SUBDIRS = introduction

include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)
INCLUDES    = $(libmesh_INCLUDE)
LIBS        = $(libmesh_LIBS) 

check_PROGRAMS =
check_SCRIPTS = 
CLEANFILES = biquadratic_h_adaptive.m *.gmv* *.xda


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
CLEANFILES += out_1 out_1_data gnupl*

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
adaptivity_ex4_DATA    = adaptivity/adaptivity_ex4/adaptivity_ex4.C
adaptivity_ex4_SOURCES = $(adaptivity_ex4_DATA) adaptivity/adaptivity_ex4/run.sh
adaptivity_ex4_LDADD   = $(top_builddir)/libmesh.la
CLEANFILES += biquadratic_h_adaptive.m

check_PROGRAMS += adaptivity_ex5
check_SCRIPTS  += adaptivity/adaptivity_ex5/run.sh
adaptivity_ex5dir      = $(examples_install_path)/adaptivity/ex5
adaptivity_ex5_DATA    = adaptivity/adaptivity_ex5/adaptivity_ex5.C transient/transient_ex1/exact_solution.C
adaptivity_ex5_SOURCES = $(adaptivity_ex5_DATA) adaptivity/adaptivity_ex5/run.sh
adaptivity_ex5_LDADD   = $(top_builddir)/libmesh.la


# adaptivity/adaptivity_ex3/adaptivity_ex3.in:
# 	$(LN_S) $(examples_source_path)/$< $@

# adaptivity/adaptivity_ex3/lshaped.xda:
# 	$(LN_S) $(examples_source_path)/$< $@

# adaptivity/adaptivity_ex3/lshaped3D.xda:
# 	$(LN_S) $(examples_source_path)/$< $@


#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)
