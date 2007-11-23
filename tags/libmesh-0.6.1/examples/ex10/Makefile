# $Id: Makefile,v 1.13 2007-10-22 22:53:18 roystgnr Exp $


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
srcfiles 	:= $(wildcard *.C) ../ex9/exact_solution.C

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	    := ./ex10-$(METHOD)

# Specify the number of timesteps to do when we call "make run"
n_timesteps := 25

all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

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
	@echo "* Running Example " $(LIBMESHRUN) $(target) \[-read_solution\] -n_timesteps $(n_timesteps) -init_timestep \[0\|$(n_timesteps)\] $(LIBMESHOPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESHRUN) $(target) -n_timesteps $(n_timesteps) \
                   -init_timestep 0 $(LIBMESHOPTIONS)
	@echo " "
	@echo "***** Finished first" $(n_timesteps) "steps, now read in" \
	       "saved solution and continue *****"
	@echo " "
	@$(LIBMESHRUN) $(target) -read_solution -n_timesteps $(n_timesteps) \
                   -init_timestep $(n_timesteps) $(LIBMESHOPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESHRUN) $(target) \[-read_solution\] -n_timesteps $(n_timesteps) -init_timestep \[0\|$(n_timesteps)\] $(LIBMESHOPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(meshdir)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(meshdir)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
