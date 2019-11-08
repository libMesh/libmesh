
# AUTOMAKE_OPTIONS = foreign # Conflicts with libMesh options
AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS  = -I m4

SUBDIRS          = src test
EXTRA_DIST       = COPYING

# Tools in the auxiliary directory

AUX_DIST         = build-aux/install-sh
AUX_DIST        += build-aux/missing
AUX_DIST        += build-aux/config.guess
AUX_DIST        += build-aux/config.sub
AUX_DIST        += build-aux/depcomp
AUX_DIST        += build-aux/ltmain.sh

# Support for pkgconfig
pkgconfigdir   = $(libdir)/pkgconfig
pkgconfig_DATA = timpi.pc

# Support 'make install prefix=/other/path' with pkgconfig
install-data-hook:
	@if (test "x$(prefix)" != "x@prefix@"); then \
	  oldprefix="@prefix@" ; \
	  newprefix="$(prefix)" ; \
	  cd $(DESTDIR)$(libdir)/pkgconfig ; \
	  for genfile in $(pkgconfig_DATA); do \
	  echo " " ; \
	    echo " *** replacing $$oldprefix" ; \
	    echo " ***      with $$newprefix" ; \
	    echo " *** in generated file $$genfile" ; \
	    echo " " ; \
	    cd $(DESTDIR)$(prefix) && cat $$genfile | $(SED) "s,$$oldprefix,$$newprefix,g" > $$genfile.replaced ; \
	    cd $(DESTDIR)$(prefix) && mv $$genfile.replaced $$genfile ; \
          done ; \
	fi

# Additional files to be deleted by 'make distclean'
DISTCLEANFILES  = _configs.sed

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES  = aclocal.m4
MAINTAINERCLEANFILES += aminclude.am
MAINTAINERCLEANFILES += autom4te.cache/*
MAINTAINERCLEANFILES += $(AUX_DIST)
MAINTAINERCLEANFILES += config.log
MAINTAINERCLEANFILES += config.status
MAINTAINERCLEANFILES += config.sub
MAINTAINERCLEANFILES += configure
MAINTAINERCLEANFILES += Makefile.in
MAINTAINERCLEANFILES += src/Makefile.in
MAINTAINERCLEANFILES += test/Makefile.in
MAINTAINERCLEANFILES += m4/libtool.m4
MAINTAINERCLEANFILES += m4/ltoptions.m4
MAINTAINERCLEANFILES += m4/ltsugar.m4
MAINTAINERCLEANFILES += m4/ltversion.m4
MAINTAINERCLEANFILES += m4/lt~obsolete.m4
