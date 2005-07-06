# makefile for TetGen

include ../../Make.common

INCLUDE = -I.

# Package name and additional flags required for library compilation
# overwrite the CXXFLAGS included from libmesh/Make.common to suppress
# a bunch of warnings
pkgname  := tetgen
CXXFLAGS := -DNDEBUG -DTRILIBRARY -fPIC

# (Relative) Path to location of shared libs
libpath := ../lib

SRC     := tetgen.C predicates.C
OBJS   	:= $(patsubst %.C, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext)
endif

ifeq ($(enable-tetgen),no)
  all::
	@echo "<<< $(pkgname) is disabled, not building >>>"
else
  all:: $(target)
endif


# Build static library
ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

# Build dynamic library
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)


clean:
	@rm -f $(OBJS) tetgen *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/lib$(pkgname).*

distclean:
	@$(MAKE) clobber
	@rm -f $(libpath)/*/lib$(pkgname).*

