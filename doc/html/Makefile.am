BUILT_SOURCES = developers.html \
                examples.html \
                externalsoftware.html \
                index.html \
                installation.html \
                publications.html \
                support.html

CLEANFILES = $(BUILT_SOURCES)

SOURCES = styles.css doxygen_stylesheet.css \
	  src/developers.html src/doxygen_header.html src/simple_footer.html \
	  src/examples.html src/index.html src/simple_header.html \
	  src/externalsoftware.html src/installation.html src/support.html \
	  src/doxygen_footer.html src/publications.html


EXTRA_DIST = examples images src styles.css doxygen_stylesheet.css



developers.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

examples.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

externalsoftware.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

index.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

installation.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

publications.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

support.html: src/$@ src/simple_header.html src/simple_footer.html
	-cat $(srcdir)/src/simple_header.html $(srcdir)/src/$@ $(srcdir)/src/simple_footer.html > $@.html.tmp
	$(AM_V_GEN)mv $@.html.tmp $@

#
# doxygen documentation - be sure to generate symlinks in include/libmesh!
doc html: doxygen
doxygen:
	@$(MAKE) -C $(top_builddir)/include/libmesh
	$(DOXYGEN) $(top_builddir)/doc/Doxyfile

#  need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -rf examples              && cp -r $(srcdir)/examples . && chmod -R u+rw examples
	-rm -rf images                && cp -r $(srcdir)/images .   && chmod -R u+rw images
	-rm -f styles.css             && cp -r $(srcdir)/styles.css  .
	-rm -f doxygen_stylesheet.css && cp -r $(srcdir)/doxygen_stylesheet.css .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += *.css .linkstamp

clean-local:
	-rm -rf doxygen images examples

else

clean-local:
	-rm -rf doxygen

endif

# Mark doxygen as a phony target to build the documentation even when the
# doxygen directory already exists.
.PHONY: doxygen

