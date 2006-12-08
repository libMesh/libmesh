# $Id: Makefile,v 1.16 2006-12-08 03:19:07 roystgnr Exp $


# The location of the mesh library
meshdir := ../..

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(meshdir)/Make.common


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
target 	   := ./ex2-$(METHOD)


all:: $(target)

# Production rules:  how to make the target - depends on library configuration 
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $@ $(LIBS) $(LDFLAGS) $(EXTERNAL_FLAGS)


# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) eqn_sys.dat

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(target)
	@echo "***************************************************************"
	@echo " "
	@$(target)
	@echo " "
	@$(target) eqn_sys.dat
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(target)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(meshdir)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard ../../include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
