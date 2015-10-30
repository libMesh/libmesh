# The location of the mesh library
LIBMESH_DIR ?= ..

# Include the library options determined by configure.  This gives
# us access to the $(METHOD) variable for building the examples.
include $(LIBMESH_DIR)/Make.common

makefilelist := $(wildcard */*ex?/Makefile */*ex??/Makefile)
dirlist := $(patsubst %/Makefile, %/, $(makefilelist))
exec_names := $(patsubst %/, %-$(METHOD), $(dirlist))

.PHONY: all clean clobber distclean examples link run .depend

# This "examples" target allows the examples to be compiled in parallel...
examples: $(exec_names)

# This "examples" target is not parallelizable by make...
# examples: all

all:
	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done

#
# Clean up the example programs
clean:
	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done

#
# Really clean up the example programs
clobber:
	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done

# Make the example program directories look
# like a clean distribution
distclean:
	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done

# Run the example programs to see if shared lib's get properly loaded.
# This target should actually be "serialized" in the sense of 'make', we don't
# want all the examples to try and run simultaneously!
run: $(wildcard examples/ex*/*.C)
	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done

#
# Rebuild the dependencies for the examples
.depend: $(wildcard examples/ex*/*.C)
	@for dir in $(dirlist); do $(MAKE) -C "$${dir}" $(MAKECMDGOALS)||exit 1; done

echo:
# @echo "dirlist = $(dirlist)"
	@echo "exec_names = $(exec_names)"

# Generic target for making an example executable. The "%" sign here is e.g.
# reduced_basis/reduced_basis_ex2.  The patsubst call elides the $(METHOD) flag.
%-$(METHOD):
	$(MAKE) -C $(patsubst %-$(METHOD), %, $@) $(MAKECMDGOALS)||exit 1
