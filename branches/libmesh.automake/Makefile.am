AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS  = -I m4
SUBDIRS          = include src

# # additional files which must be included in 'make dist'
# EXTRA_DIST = \
# 	   LICENSE \
# 	   data \
# 	   dist_version

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/missing

# Additional files to be deleted by 'make distclean'
DISTCLEANFILES =  \
	   _configs.sed

# # Additional files to be deleted by 'make maintainer-clean'
# MAINTAINERCLEANFILES =  \
# 	   $(AUX_DIST)                              \
# 	   Makefile.in                              \
# 	   aclocal.m4                               \
# 	   aminclude.am                             \
# 	   autom4te.cache/*                         \
# 	   config.h	                            \
# 	   config.h.in	                            \
# 	   config.log                               \
# 	   config.status                            \
# 	   config.sub                               \
# 	   configure                                \
# 	   data/Makefile.in                         \
# 	   src/Makefile.in                          \
# 	   src/physics/properties/Makefile.in       \
# 	   test/Makefile.in                         \
# 	   test/common/Makefile.in 		    \
# 	   test/comp_ns/Makefile.in		    \
# 	   test/unit/Makefile.in                    


#
.PHONY: doc svnexpand

#
# doxygen documentation
doc:
	mkdir -p $(top_srcdir)/doc/html
	$(DOXYGEN) $(top_builddir)/doc/Doxyfile

#
# Mark source files for SVN keyword expansion
svnexpand:
	svn -R propset svn:keywords "Date Author Revision HeadURL Id" $(top_srcdir)/include $(top_srcdir)/src


#
# special things to do when running 'make dist'
dist-hook:
	rm -rf `find $(distdir) -name .svn`
	rm -rf `find $(distdir)/test -name logs`
	rm -rf `find $(distdir)/test -name output`
	rm -rf `find $(distdir)/contrib/ablation -type d`
	rm -rf `find $(distdir)/contrib/chemay -type d`
	rm -rf `find $(distdir)/contrib/grvy -type d`
	rm -rf `find $(distdir)/contrib/masa -type d`
	rm -rf $(distdir)/contrib/install


#
# Run unit tests only
check-unit:
	@$(MAKE) && cd $(top_builddir)/test/unit && $(MAKE) check

#
# Run regression tests only.  Don't forget to make compare first though!
check-regression:
	@$(MAKE) 
	@cd $(top_builddir)/test/common  && $(MAKE) check 
	@cd $(top_builddir)/test/comp_ns && $(MAKE) check
