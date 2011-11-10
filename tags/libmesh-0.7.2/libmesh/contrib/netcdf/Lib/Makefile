include ../../../Make.common

libmesh_INCLUDE = -I.

# Package name and additional flags required for library compilation
pkgname  := netcdf

# (Relative) Path to location of shared libs
libpath := ../lib

# The source files for netcdf
SRC     := attr.c ncx.c putget.c dim.c error.c libvers.c nc.c string.c v1hpg.c var.c posixio.c v2i.c
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext)
endif


ifeq ($(enable-netcdf),no)
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
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/lib$(pkgname).*

distclean:
	@$(MAKE) clobber
	@rm -f $(libpath)/*/lib$(pkgname).*

# NetCDF file dependences.  Generated with:
# perl ../../bin/make_dependencies.pl -I. -S\$\(obj-suffix\) *.c

attr.$(obj-suffix):\
    attr.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    rnd.h
dim.$(obj-suffix):\
    dim.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    rnd.h
error.$(obj-suffix):\
    error.c\
    netcdf.h
ffio.$(obj-suffix):\
    ffio.c\
    fbits.h\
    ncio.h\
    netcdf.h\
    rnd.h
libvers.$(obj-suffix):\
    libvers.c\
    netcdf.h
nc.$(obj-suffix):\
    nc.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    rnd.h
ncx.$(obj-suffix):\
    ncx.c\
    ncx.h\
    rnd.h
posixio.$(obj-suffix):\
    posixio.c\
    fbits.h\
    ncio.h\
    netcdf.h\
    rnd.h
putget.$(obj-suffix):\
    putget.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    onstack.h\
    rnd.h
string.$(obj-suffix):\
    string.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    rnd.h
t_nc.$(obj-suffix):\
    t_nc.c\
    netcdf.h
v1hpg.$(obj-suffix):\
    v1hpg.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    rnd.h
v2i.$(obj-suffix):\
    v2i.c\
    onstack.h
var.$(obj-suffix):\
    var.c\
    fbits.h\
    nc.h\
    ncio.h\
    ncx.h\
    netcdf.h\
    rnd.h
