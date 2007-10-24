# $Id: Makefile,v 1.6 2005-08-23 19:31:17 knezed01 Exp $


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
srcfiles 	:= $(wildcard *.C) ../ex9/exact_solution.C

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	    := ./ex10

# Specify the number of timesteps to do when we call "make run"
n_timesteps := 25

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
	@rm -f $(target) out_*.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o


run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(target)
	@echo "***************************************************************"
	@echo " "
	@$(target) -n_timesteps $(n_timesteps) -init_timestep 0
	@echo " "
	@echo "***** Finished first" $(n_timesteps) "steps, now read in" \
	       "saved solution and continue *****"
	@echo " "
	@$(target) -read_solution -n_timesteps $(n_timesteps) \
                   -init_timestep $(n_timesteps)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(target)
	@echo "***************************************************************"
#
# Dependencies
#
.depend:
	@$(perl) $(meshdir)/contrib/bin/make_dependencies.pl $(foreach i, $(wildcard ../../include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
