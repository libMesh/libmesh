BUILT_SOURCES =
DISTCLEANFILES =

#
# doxygen documentation - be sure to generate symlinks in include/libmesh!
doc html:
	@cd $(top_builddir)/include/libmesh && $(MAKE)
	$(DOXYGEN) $(top_builddir)/doc/Doxyfile

clean-local:
	-rm -rf doxygen

#  need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f developers.html       && $(LN_S) $(srcdir)/developers.html .
	-rm -f examples              && $(LN_S) $(srcdir)/examples .
	-rm -f examples.html         && $(LN_S) $(srcdir)/examples.html .
	-rm -f externalsoftware.html && $(LN_S) $(srcdir)/externalsoftware.html .
	-rm -f images                && $(LN_S) $(srcdir)/images .
	-rm -f index.html            && $(LN_S) $(srcdir)/index.html .
	-rm -f installation.html     && $(LN_S) $(srcdir)/installation.html .
	-rm -f publications.html     && $(LN_S) $(srcdir)/publications.html .
	-rm -f styles.css            && $(LN_S) $(srcdir)/styles.css  .
	-rm -f support.html          && $(LN_S) $(srcdir)/support.html .
	$(AM_V_GEN)touch .linkstamp



  BUILT_SOURCES += \
	developers.html examples examples.html externalsoftware.html \
	images index.html installation.html publications.html \
	styles.css support.html

endif

DISTCLEANFILES += $(BUILT_SOURCES)
