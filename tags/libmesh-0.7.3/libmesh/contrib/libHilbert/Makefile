include ../../Make.common

libmesh_INCLUDE = -I./include 

# Package name and additional flags required for library compilation
pkgname  := Hilbert

# (Relative) Path to location of shared libs
libpath := ../lib

# List of source files
SRC     := $(shell find ./src -name "*.C")
OBJS   	:= $(patsubst %.C, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext)
endif


ifeq ($(enable-libhilbert),no)
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
	@$(libmesh_CC) $(libmesh_CSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)


clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f src/*.o
	@rm -f $(EXTERNAL_LIBDIR)/lib$(pkgname).*

distclean:
	@$(MAKE) clobber
	@rm -f $(libpath)/*/lib$(pkgname).*

# Dependencies for libHilbert files
src/%.$(obj-suffix): ../../include/base/libmesh_config.h

src/BigBitVec.$(obj-suffix):\
    src/BigBitVec.C\
    include/Hilbert/BigBitVec.hpp\
    include/Hilbert/Common.hpp\
    include/Hilbert/FixBitVec.hpp
src/FixBitVec.$(obj-suffix):\
    src/FixBitVec.C\
    include/Hilbert/Common.hpp\
    include/Hilbert/FixBitVec.hpp
src/Hilbert.$(obj-suffix):\
    src/Hilbert.C\
    include/Hilbert.hpp\
    include/Hilbert/Algorithm.hpp\
    include/Hilbert/BigBitVec.hpp\
    include/Hilbert/Common.hpp\
    include/Hilbert/FixBitVec.hpp\
    include/Hilbert/GetBits.hpp\
    include/Hilbert/GetLocation.hpp\
    include/Hilbert/GrayCodeRank.hpp\
    include/Hilbert/SetBits.hpp\
    include/Hilbert/SetLocation.hpp
src/hilbert_mpi_ops.$(obj-suffix):\
    src/hilbert_mpi_ops.C\
    include/Hilbert.hpp\
    include/Hilbert/BigBitVec.hpp\
    include/Hilbert/Common.hpp\
    include/Hilbert/FixBitVec.hpp\
    include/hilbert.h
