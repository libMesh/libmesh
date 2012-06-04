# The location of the mesh library
LIBMESH_DIR ?= ../../..

include ../../../Make.common

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))

libmesh_INCLUDE += -I../GKlib -I. 


.PHONY: clean clobber distclean

# Metis does not cast pointers to (void **) when using GKFree.
# This causes a warning of "Incompatible pointer type" when a
# prototype exists, and an "Implicit declaration" when a prototype
# doesn't exist (see proto.h) which is the way it stands now.
# We don't want to see the warnings so we will filter out the
# '-Wimplicit' flag when compiling this package.
libmesh_CFLAGS := $(filter-out '-Wimplicit',$(libmesh_CFLAGS))

target := $(EXTERNAL_LIBDIR)/libmetis$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/libmetis$(shared_libext)
endif

ifeq ($(enable-metis),no)
  all::
	@echo "<<< metis is disabled, not building >>>"
else
  all:: $(target)
endif

ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/libmetis$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/libmetis$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

$(EXTERNAL_LIBDIR)/libmetis$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(libmesh_CC) $(libmesh_CSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)

clean:
	@rm -f $(OBJS) *~ .depend

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/libmetis.*

distclean:
	@$(MAKE) clobber
	@rm -f ../../lib/*/libmetis.*

# include the dependency list
include .depend


#
# Dependencies
#
.depend: $(SRC) *.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I.  "-S\$$(obj-suffix)" $(SRC) > .depend

###############################################################################
