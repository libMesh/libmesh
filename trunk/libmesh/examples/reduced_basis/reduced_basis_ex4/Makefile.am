include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_INCLUDE) -I$(top_builddir)/include
LIBS         = $(libmesh_LIBS)

if LIBMESH_ENABLE_FPARSER
  AM_CPPFLAGS += -I$(top_srcdir)/contrib/fparser
endif

check_PROGRAMS            = reduced_basis_ex4
check_SCRIPTS             = run.sh
reduced_basis_ex4dir      = $(examples_install_path)/reduced_basis/ex4
reduced_basis_ex4_DATA    = reduced_basis_ex4.C assembly.h eim_classes.h rb_classes.h reduced_basis_ex4.in eim.in rb.in
reduced_basis_ex4_SOURCES = $(reduced_basis_ex4_DATA) run.sh
reduced_basis_ex4_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = RB_sol.e rb_data/RB_F_019.xdr rb_data/RB_F_018.xdr rb_data/RB_F_017.xdr rb_data/RB_F_016.xdr rb_data/RB_F_015.xdr rb_data/RB_F_014.xdr rb_data/RB_F_013.xdr rb_data/RB_F_012.xdr rb_data/RB_F_011.xdr rb_data/RB_F_010.xdr rb_data/RB_F_009.xdr rb_data/RB_F_008.xdr rb_data/RB_F_007.xdr rb_data/RB_F_006.xdr rb_data/RB_F_005.xdr rb_data/RB_F_004.xdr rb_data/RB_F_003.xdr rb_data/RB_F_002.xdr rb_data/RB_F_001.xdr rb_data/RB_F_000.xdr rb_data/RB_A_000.xdr rb_data/parameter_ranges.xdr rb_data/n_bfs.xdr rb_data/greedy_params.xdr rb_data/Fq_innerprods.xdr rb_data/Fq_Aq_innerprods.xdr rb_data/bf_header.xdr rb_data/bf9.xdr rb_data/bf8.xdr rb_data/bf7.xdr rb_data/bf6.xdr rb_data/bf5.xdr rb_data/bf4.xdr rb_data/bf3.xdr rb_data/bf2.xdr rb_data/bf14.xdr rb_data/bf13.xdr rb_data/bf12.xdr rb_data/bf11.xdr rb_data/bf10.xdr rb_data/bf1.xdr rb_data/bf0.xdr rb_data/Aq_Aq_innerprods.xdr eim_data/RB_inner_product_matrix.xdr eim_data/parameter_ranges.xdr eim_data/n_bfs.xdr eim_data/interpolation_points_var.xdr eim_data/interpolation_points.xdr eim_data/interpolation_matrix.xdr eim_data/greedy_params.xdr eim_data/Fq_innerprods.xdr eim_data/Fq_Aq_innerprods.xdr eim_data/extra_interpolation_point_var.xdr eim_data/extra_interpolation_point.xdr eim_data/extra_interpolation_matrix_row.xdr eim_data/bf_header.xdr eim_data/bf9.xdr eim_data/bf8.xdr eim_data/bf7.xdr eim_data/bf6.xdr eim_data/bf5.xdr eim_data/bf4.xdr eim_data/bf3.xdr eim_data/bf2.xdr eim_data/bf19.xdr eim_data/bf18.xdr eim_data/bf17.xdr eim_data/bf16.xdr eim_data/bf15.xdr eim_data/bf14.xdr eim_data/bf13.xdr eim_data/bf12.xdr eim_data/bf11.xdr eim_data/bf10.xdr eim_data/bf1.xdr eim_data/bf0.xdr eim_data/Aq_Aq_innerprods.xdr 

clean-local:
	-rm -rf rb_data eim_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex4.in && $(LN_S) $(srcdir)/reduced_basis_ex4.in .
	-rm -f eim.in               && $(LN_S) $(srcdir)/eim.in .
	-rm -f rb.in                && $(LN_S) $(srcdir)/rb.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex4.in eim.in rb.in .linkstamp
endif

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


# how to make documentation
CLEANFILES += $(check_PROGRAMS).php stdout.log
doc: $(check_PROGRAMS).php
html: $(check_PROGRAMS).php
$(check_PROGRAMS).php: $(check_PROGRAMS) \
                        $(top_srcdir)/contrib/bin/ex2html.sh \
                       Makefile
	$(MAKE) -s run > stdout.log
	$(AM_V_GEN) PATH=$(abs_top_srcdir)/contrib/bin:$(PATH) \
	  ex2html.sh $(check_PROGRAMS) $(abs_srcdir)

