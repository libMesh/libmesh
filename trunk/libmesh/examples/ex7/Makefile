# $Id: Makefile,v 1.8 2003-04-08 22:54:09 benkirk Exp $


# The location of the mesh library
meshdir := ../..

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(meshdir)/Make.common

# include the dependency list
include .depend


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ex7


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects) $(meshdir)/include/mesh_config.h
	@echo "Linking "$@"..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $@ $(LIBS) $(LDFLAGS)


echo:
	@echo "obj= $(objects)"


# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) out0*.gmv eqn_sys.dat

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

ifeq ($(enable-complex),yes)
run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(target)
	@echo "***************************************************************"
	@echo " "
	@$(target) -f .7
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(target)
	@echo "***************************************************************"
else
run: $(target)
	@echo "***************************************************************"
	@echo "*** Skipping Example " $(target) ", only good with --enable-complex"
	@echo "***************************************************************"
endif
#
# Dependencies
#
.depend:
	@$(perl) $(meshdir)/contrib/bin/make_dependencies.pl -I$(meshdir)/include "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
