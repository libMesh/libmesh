# $Id: Makefile,v 1.7 2003-01-31 16:47:26 benkirk Exp $
#
# This is the Makefile for the libMesh library and helper
# applications.  This file is specific to the project.
# See the file Make.common for path configurations.
# Make.common is what should be included by applications
# using the library.


# include the library options determined by configure
include Make.common


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# header files
headerfiles 	:= $(wildcard include/*.h)

#
# source files
srcfiles 	:= $(wildcard src/*/*.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

#
# logged files -- all files you might want to log information for
loggedfiles     := $(srcfiles) $(filter-out include/mesh_config.h, $(headerfiles))
###############################################################################



###############################################################################
# Target:
#
targ_dir   := lib/$(hosttype)_$(METHOD)

targ 	   := $(pwd)/$(targ_dir)
target 	   := $(targ)/libmesh.a

ifeq ($(enable-shared),yes) 
  target   := $(targ)/libmesh.so
endif

ifeq ($(syn-mode),on)
  all:: $(objects)
else
  all:: $(target)
endif
###############################################################################

#
# 
#
.PHONY: clean clobber distclean doc log cvsweb TODO

#
# static library
#
$(targ)/libmesh.a: $(objects)
	@$(MAKE) -C contrib
	@$(shell mkdir -p $(targ_dir))
	@echo ""
	@echo "Linking "$@
	@$(AR) rv $(target) $(objects)

#
# shared library
#
$(targ)/libmesh.so: $(objects)
	@$(MAKE) -C contrib
	@$(shell mkdir -p $(targ_dir))
	@echo ""
	@echo "Linking "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $(target) $(objects)

#
# Build the examples
#
examples: $(target)
	@$(MAKE) -C examples


#
# Useful for checking make rules
#
echo:
	@echo -e "Source Files:\n$(srcfiles)\n"
	@echo -e "Object Files:\n$(objects)\n"
	@echo -e "Target:\n$(target)\n"
	@echo -e "CFLAGS:\n$(CXXFLAGS)\n"
	@echo -e "CXXFLAGS:\n$(CXXFLAGS)\n"
	@echo -e "INCLUDE:\n$(INCLUDE)\n"
	@echo -e "LIBS:\n$(LIBS)\n"

#
# Remove project object files for the current mode
#
clean:
	@rm -f *~ include/*~ src/*/*~ src/*/*.$(obj-suffix)

#
# Make clean, remove contributed objects, and remove binaries
#
clobber:
	@$(MAKE) clean
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@rm -rf $(targ_dir) bin/grid2grid bin/meshtool bin/testexodus

#
# Make clobber, remove documentation, removes all libraries
# Should restore to a pristine state, except for files you
# have added
distclean:
	$(MAKE) clobber
	@$(MAKE) -C contrib $(MAKECMDGOALS)
	@$(MAKE) -C examples $(MAKECMDGOALS)
	@rm -rf doc/html doc/latex doc/kdoc/*.html \
               doc/man/man3 doc/cvshtml/*.html doc/cvshtml/diff \
	       src/*/*.o src/*/*.go src/*/*.pgo \
	       lib/*_opt lib/*_dbg lib/*_pro



#
# kdoc and doxygen documentation
#
doc:
	$(doxygen) ./doc/Doxyfile



log: $(loggedfiles)
	cvs log $(loggedfiles) > cvs_log

#
# CVS documentation
#
cvsweb:
	./contrib/bin/cvs2html -f -p -o doc/cvshtml/index.html -v -a -b -n 3 -C crono.html

#
# Meshtool utility program
#
meshtool: $(target) src/apps/meshtool.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/apps/meshtool.cc $< -o bin/$@ $(LIBS) $(LDFLAGS)

#
# Read_Dat utility program
#
read_dat: $(target) src/apps/read_dat.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/apps/read_dat.cc $< -o bin/$@ $(LIBS) $(LDFLAGS)

#
# foo
#
foo: $(target) src/apps/foo.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/apps/foo.cc $< -o bin/$@ $(LIBS) $(LDFLAGS)


#
# grid2grid
#
grid2grid: $(target) src/apps/grid2grid.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/apps/grid2grid.cc $< -o bin/$@ $(LIBS) $(LDFLAGS)

#
# Testexodus FE program
#
testexodus: $(target) src/apps/testexodus.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/apps/testexodus.cc $< -o bin/$@ $(LIBS) $(LDFLAGS)

#
# Make a TODO list
#
TODO:
	@egrep -i '// *todo' $(srcfiles) $(includes) \
	| perl -pi -e 's#(.*)(TODO:?)(\[.+\])#\3 \1\2\3#i;' \
	| perl -pi -e 's#\s*//\s*TODO:\s*(\[.+\])\s*#\n\1     #i;'



#
# Dependencies
#
.depend:
	@$(perl) ./contrib/bin/make_dependencies.pl -I./include "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"



#
# Make.common target.  Tell the user to run configue.
#
Make.common:
	@echo -e "You must run ./configure first!"
	exit 1

###############################################################################
include .depend


# Local Variables:
# mode: makefile
# End:
