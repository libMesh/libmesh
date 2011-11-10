# makefile for TetGen

include ../../Make.common

INCLUDE = -I.

# Package name and additional flags required for library compilation
pkgname  := tetgen

# Set compilation flags
#libmesh_CXXFLAGS := -DNDEBUG -DTETLIBRARY -fPIC

# [JWP] Should we not inherit whatever libmesh_CXXFLAGS were set by Make.common and add our own?
libmesh_CXXFLAGS += -DTETLIBRARY -fPIC

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
	@$(libmesh_CXX) $(libmesh_CXXSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)


clean:
	@rm -f $(OBJS) tetgen *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/lib$(pkgname).*

distclean:
	@$(MAKE) clobber
	@rm -f $(libpath)/*/lib$(pkgname).*


# Dependencies for tetgen object file
tetgen.$(obj-suffix):\
    tetgen.C\
    tetgen.h

predicates.$(obj-suffix):\
    predicates.C\
    tetgen.h
