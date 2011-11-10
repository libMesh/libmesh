include ../../../Make.common

libmesh_INCLUDE = -I../../exodusii/Lib/include -I. -I../../netcdf/Lib

# Package name and additional flags required for library compilation
pkgname  := nemesis

# (Relative) Path to location of shared libs
libpath := ../../lib

# The names of all the source files, except ne_test, which contains a main()
SRC     := $(shell find ./ -name "*.c" -a ! -name "ne_test.c")
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext)
endif


ifeq ($(enable-nemesis),no)
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

# Almost all the files depend on the following headers
nemesis_header_deps := ne_nemesisI.h ne_nemesisI_int.h ../../exodusii/Lib/include/exodusII.h ../../exodusii/Lib/include/exodusII_ext.h ../../exodusii/Lib/include/exodusII_int.h ../../netcdf/Lib/netcdf.h 

# Script to generate dependencies for all source files
# for i in *.c; do echo "`basename $i c`\$(obj-suffix): $i \$(nemesis_header_deps)"; done

ne_gcmp.$(obj-suffix): ne_gcmp.c $(nemesis_header_deps)
ne_gebig.$(obj-suffix): ne_gebig.c $(nemesis_header_deps)
ne_gelcm.$(obj-suffix): ne_gelcm.c $(nemesis_header_deps)
ne_gelm.$(obj-suffix): ne_gelm.c $(nemesis_header_deps)
ne_gelt.$(obj-suffix): ne_gelt.c $(nemesis_header_deps)
ne_gennm.$(obj-suffix): ne_gennm.c $(nemesis_header_deps)
ne_genv.$(obj-suffix): ne_genv.c $(nemesis_header_deps)
ne_giinf.$(obj-suffix): ne_giinf.c $(nemesis_header_deps)
ne_ginig.$(obj-suffix): ne_ginig.c $(nemesis_header_deps)
ne_glbp.$(obj-suffix): ne_glbp.c $(nemesis_header_deps)
ne_gncm.$(obj-suffix): ne_gncm.c $(nemesis_header_deps)
ne_gncor.$(obj-suffix): ne_gncor.c $(nemesis_header_deps)
ne_gneat.$(obj-suffix): ne_gneat.c $(nemesis_header_deps)
ne_gnelc.$(obj-suffix): ne_gnelc.c $(nemesis_header_deps)
ne_gnm.$(obj-suffix): ne_gnm.c $(nemesis_header_deps)
ne_gnnnm.$(obj-suffix): ne_gnnnm.c $(nemesis_header_deps)
ne_gnns.$(obj-suffix): ne_gnns.c $(nemesis_header_deps)
ne_gnnsd.$(obj-suffix): ne_gnnsd.c $(nemesis_header_deps)
ne_gnnv.$(obj-suffix): ne_gnnv.c $(nemesis_header_deps)
ne_gnspg.$(obj-suffix): ne_gnspg.c $(nemesis_header_deps)
ne_gnss.$(obj-suffix): ne_gnss.c $(nemesis_header_deps)
ne_gnssd.$(obj-suffix): ne_gnssd.c $(nemesis_header_deps)
ne_gsspg.$(obj-suffix): ne_gsspg.c $(nemesis_header_deps)
ne_pcmp.$(obj-suffix): ne_pcmp.c $(nemesis_header_deps)
ne_pcmpc.$(obj-suffix): ne_pcmpc.c $(nemesis_header_deps)
ne_pebig.$(obj-suffix): ne_pebig.c $(nemesis_header_deps)
ne_pelcm.$(obj-suffix): ne_pelcm.c $(nemesis_header_deps)
ne_pelm.$(obj-suffix): ne_pelm.c $(nemesis_header_deps)
ne_pennm.$(obj-suffix): ne_pennm.c $(nemesis_header_deps)
ne_pevs.$(obj-suffix): ne_pevs.c $(nemesis_header_deps)
ne_piinf.$(obj-suffix): ne_piinf.c $(nemesis_header_deps)
ne_pinig.$(obj-suffix): ne_pinig.c $(nemesis_header_deps)
ne_plbp.$(obj-suffix): ne_plbp.c $(nemesis_header_deps)
ne_plbpc.$(obj-suffix): ne_plbpc.c $(nemesis_header_deps)
ne_pncm.$(obj-suffix): ne_pncm.c $(nemesis_header_deps)
ne_pncor.$(obj-suffix): ne_pncor.c $(nemesis_header_deps)
ne_pneat.$(obj-suffix): ne_pneat.c $(nemesis_header_deps)
ne_pnelc.$(obj-suffix): ne_pnelc.c $(nemesis_header_deps)
ne_pnm.$(obj-suffix): ne_pnm.c $(nemesis_header_deps)
ne_pnnnm.$(obj-suffix): ne_pnnnm.c $(nemesis_header_deps)
ne_pnns.$(obj-suffix): ne_pnns.c $(nemesis_header_deps)
ne_pnnsd.$(obj-suffix): ne_pnnsd.c $(nemesis_header_deps)
ne_pnspg.$(obj-suffix): ne_pnspg.c $(nemesis_header_deps)
ne_pnss.$(obj-suffix): ne_pnss.c $(nemesis_header_deps)
ne_pnssd.$(obj-suffix): ne_pnssd.c $(nemesis_header_deps)
ne_pnvs.$(obj-suffix): ne_pnvs.c $(nemesis_header_deps)
ne_psspg.$(obj-suffix): ne_psspg.c $(nemesis_header_deps)
ne_test.$(obj-suffix): ne_test.c $(nemesis_header_deps) # Not built in general, contains a main program
ne_util.$(obj-suffix): ne_util.c $(nemesis_header_deps)

