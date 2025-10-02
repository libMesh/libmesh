pkg_sources = $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/poly2tri.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/utils.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.h
pkg_libs    =

EXTRA_DIST = $(top_srcdir)/contrib/poly2tri/poly2tri/AUTHORS \
             $(top_srcdir)/contrib/poly2tri/poly2tri/LICENSE \
             $(top_srcdir)/contrib/poly2tri/poly2tri/README.md

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common

poly2tri_headers = poly2tri/poly2tri.h \
                   poly2tri/common/shapes.h \
                   poly2tri/common/utils.h \
                   poly2tri/sweep/advancing_front.h \
                   poly2tri/sweep/cdt.h \
                   poly2tri/sweep/sweep_context.h \
                   poly2tri/sweep/sweep.h

# includedir is an automake variable, and we are just using its default value here
includedir             = $(prefix)/include

# nobase_include_HEADERS are headers which should be _installed_ while maintaining relative paths.
# nodist_include_HEADERS are headers which should _not_ be included as part of the source distribution.
#
# In our case, we _do_ want to install the modified poly2tri headers
# that we create during the build process, but we don't want to
# distribute them.  If you separately specify the same set of headers
# as both nobase_ and nodist_, this apparently creates a conflict for
# automake, and it ends up including the files in the distribution
# anyway. But, according to Google AI: "If you need both
# behaviors--keeping the subdirectory structure AND not
# distributing--the correct variable name is
# nobase_nodist_include_HEADERS. The nobase_ prefix must come first."
nobase_nodist_include_HEADERS = $(poly2tri_headers)

BUILT_SOURCES = $(poly2tri_headers)

clean-local:
	rm -f $(BUILT_SOURCES)

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
