# $Id$

# The location of the mesh library
LIBMESH_DIR ?= ..

# Include the library options determined by configure.  This gives
# us access to the $(METHOD) variable for building the examples.
include $(LIBMESH_DIR)/Make.common

dirlist := $(wildcard */*ex?/ */*ex??/)


# Build the example programs.  This target has been deprecated
# in favor of one that make can parallelize.
# examples:
# 	@echo "----------------------------------------------"
# 	@echo "--------- Building Example Programs ----------"
# 	@echo "----------------------------------------------"
# 	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done
# 	@echo "----------------------------------------------"
# 	@echo "------- Done Building Example Programs -------"
# 	@echo "----------------------------------------------"


.PHONY: all clean clobber distclean examples link run .depend

# This target allows the examples to be compiled in parallel
examples: all

all:
	@for dir in $(dirlist); do if [ -f "$${dir}"/Makefile ]; then $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; fi; done

#
# Clean up the example programs
clean:
	@for dir in $(dirlist); do if [ -f "$${dir}"/Makefile ]; then $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; fi; done

#
# Really clean up the example programs
clobber:
	@for dir in $(dirlist); do if [ -f "$${dir}"/Makefile ]; then $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; fi; done

# Make the example program directories look
# like a clean distribution
distclean:
	@for dir in $(dirlist); do if [ -f "$${dir}"/Makefile ]; then $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; fi; done

# Link the example programs (deprecated, use "make examples" instead)
link: all

# Run the example programs to see if shared lib's get properly loaded.
# This target should actually be "serialized" in the sense of 'make', we don't
# want all the examples to try and run simultaneously!
run: $(wildcard examples/ex*/*.C)
	@for dir in $(dirlist); do if [ -f "$${dir}"/Makefile ]; then $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; fi; done

#
# Rebuild the dependencies for the examples
.depend: $(wildcard examples/ex*/*.C)
	@for dir in $(dirlist); do if [ -f "$${dir}"/Makefile ]; then $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; fi; done

echo:
	@echo "dirlist = $(dirlist)"
