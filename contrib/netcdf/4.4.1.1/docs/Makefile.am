## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005-2011, see the COPYRIGHT file for more information.
# This file builds the netcdf documentation.


# These files will be included with the dist.
EXTRA_DIST = netcdf.m4 DoxygenLayout.xml Doxyfile.in footer.html \
    mainpage.dox tutorial.dox dispatch.dox \
    guide.dox types.dox cdl.dox \
    architecture.dox internal.dox windows-binaries.md \
    building-with-cmake.md CMakeLists.txt \
    groups.dox install.md notes.md install-fortran.md \
    all-error-codes.md cmake_faq.md credits.md software.md \
    auth.html obsolete/fan_utils.html bestpractices.md

# Turn off parallel builds in this directory.
.NOTPARALLEL:

# To build documentation you must have M4 in their path.
M4=m4

SUBDIRS = images

if BUILD_DOCS

# Copy man pages.
#directory = $(top_srcdir)/docs/man/man3/
#dist_man_MANS = $(docsdirectory)/man_page_1.3 $(directory)/man_page_2.3
# $(directory)/man_page_1.3: doxyfile.stamp
# $(directory)/man_page_2.3: doxyfile.stamp

# Note: in order to work with distcheck,
# the Doxyfile needs to be preprocessed
# to insert actual location of $(srcdir)
#BUILT_SOURCES = Doxyfile.tmp

#Doxyfile.tmp: Doxyfile
#	sed -e 's|$$[({]rootdir[})]|$(abs_top_srcdir)|g' <Doxyfile > ./Doxyfile.tmp

all-local: doxyfile.stamp
clean-local:
	rm -rf $(top_builddir)/docs/man
	rm -rf $(top_builddir)/docs/html
	rm -rf $(top_builddir)/docs/latex
endif

# Timestamp to prevent rebuilds.
# We must do this twice. The first time
# builds the tag files. The second time
# includes them in the documentation.

doxyfile.stamp:
	$(DOXYGEN) Doxyfile \
	cp auth.html obsolete/fan_utils.html html

CLEANFILES = doxyfile.stamp html latex man

# This builds the docs from source, if necessary, and tars up
# everything needed for the website. Run this and copy the resulting
# tarball to the /contents/netcdf/docs directory to update the on-line
# docs.
web-tarball: doxyfile.stamp
	cd html; tar cf ../netcdf_docs.tar *
	gzip -f netcdf_docs.tar

# When oc2/auth.html.in is changed, you should generate auth.html
# using the following process.
oc2-auth:
	cat ${top_srcdir}/oc2/auth.html.in \
	| sed -e '/<OC>/d' \
	| sed -e 's|^<NC>||' \
	| sed -e 's|zz|netcdf|g' -e 's|ZZ|netCDF|g' \
	| sed -e '/stylesheet/r${top_srcdir}/oc2/oc.css' -e '/stylesheet/d' \
	> auth.html
