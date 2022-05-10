pkg_sources = $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.cc
pkg_libs    =

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common

poly2tri_headers = poly2tri/poly2tri.h \
                   poly2tri/common/shapes.h \
                   poly2tri/common/utils.h \
                   poly2tri/sweep/advancing_front.h \
                   poly2tri/sweep/cdt.h \
                   poly2tri/sweep/sweep_context.h \
                   poly2tri/sweep/sweep.h

nobase_include_HEADERS = $(poly2tri_headers)

poly2tri/poly2tri.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/poly2tri.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/
	$(AM_V_GEN) $(SED) -e 's#include "#include "poly2tri/#' < $< > $@

poly2tri/common/shapes.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/common/
	$(AM_V_GEN) rm -f $@ && cp -f $< $@

poly2tri/common/utils.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/utils.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/common/
	$(AM_V_GEN) $(SED) -e 's#include "#include "poly2tri/#' < $< > $@

poly2tri/sweep/advancing_front.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) $(SED) -e 's#include "../#include "poly2tri/#' < $< > $@

poly2tri/sweep/cdt.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) $(SED) -e 's#include "#include "poly2tri/sweep/#' < $< > $@

poly2tri/sweep/sweep_context.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) rm -f $@ && cp -f $< $@

poly2tri/sweep/sweep.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) rm -f $@ && cp -f $< $@
