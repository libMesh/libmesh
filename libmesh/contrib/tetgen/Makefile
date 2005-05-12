# makefile for TetGen

include ../../Make.common

INCLUDE = -I.

# Package name and additional flags required for library compilation
# overwrite the CXXFLAGS included from libmesh/Make.common to suppress
# a bunch of warnings
pkgname  := tetgen
CXXFLAGS := -DNDEBUG -DTRILIBRARY

# (Relative) Path to location of shared libs
libpath := ../lib

SRC     := tetgen.C predicates.C
OBJS   	:= $(patsubst %.C, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

target := $(EXTERNAL_LIBDIR)/lib$(pkgname).a

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(SHARED_LIBEXT)
endif

ifeq ($(enable-triangle),no)
  all::
	@echo "<<< $(pkgname) is disabled, not building >>>"
else
  all:: $(target)
endif


# Build static library
ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/lib$(pkgname).a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/lib$(pkgname).a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

# Build dynamic library
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(SHARED_LIBEXT): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)


clean:
	@rm -f $(OBJS) tetgen *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/lib$(pkgname).*

distclean:
	@$(MAKE) clobber
	@rm -f $(libpath)/*/lib$(pkgname).*

