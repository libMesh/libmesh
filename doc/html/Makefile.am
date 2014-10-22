BUILT_SOURCES =
CLEANFILES =
DISTCLEANFILES =

EXTRA_DIST = \
	developers.html examples examples.html externalsoftware.html \
	images index.html installation.html publications.html \
	styles.css doxygen_stylesheet.css support.html src

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
	-rm -f developers.html        && cp -r $(srcdir)/developers.html .
	-rm -f examples               && cp -r $(srcdir)/examples .
	-rm -f examples.html          && cp -r $(srcdir)/examples.html .
	-rm -f externalsoftware.html  && cp -r $(srcdir)/externalsoftware.html .
	-rm -f images                 && cp -r $(srcdir)/images .
	-rm -f index.html             && cp -r $(srcdir)/index.html .
	-rm -f installation.html      && cp -r $(srcdir)/installation.html .
	-rm -f publications.html      && cp -r $(srcdir)/publications.html .
	-rm -f styles.css             && cp -r $(srcdir)/styles.css  .
	-rm -f doxygen_stylesheet.css && cp -r $(srcdir)/doxygen_stylesheet.css .
	-rm -f support.html           && cp -r $(srcdir)/support.html .
	$(AM_V_GEN)touch .linkstamp

  BUILT_SOURCES += \
	developers.html examples examples.html externalsoftware.html \
	images index.html installation.html publications.html \
	styles.css doxygen_stylesheet.css support.html

endif

CLEANFILES += $(BUILT_SOURCES)
