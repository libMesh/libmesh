AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = metaphysicl.dox metaphysicl.page txt_common fig_common

DISTCLEANFILES  = metaphysicl.dox
DISTCLEANFILES += Makefile

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

