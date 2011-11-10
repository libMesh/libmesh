include ../../../../Make.common

libmesh_INCLUDE = -I../include -I../../../netcdf/Lib

# Package name and additional flags required for library compilation
pkgname  := exodusii

# (Relative) Path to location of shared libs
libpath := ../../lib

# Only one source file for triangle
SRC     := $(shell find ./ -name "*.c")
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext)
endif


ifeq ($(enable-exodus),no)
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

# Dependencies for ExodusII files.  Every c file depends on the following headers:
ex_header_deps := ../../../netcdf/Lib/netcdf.h ../include/exodusII.h ../include/exodusII_ext.h ../include/exodusII_int.h

# Script to generate the following lines automatically, unfortunately
# I don't think there is a nice Makfile way to specify all of these in 
# one line?
# for i in *.c; do echo "`basename $i c`\$(obj-suffix): $i \$(ex_header_deps)"; done

ex_conv.$(obj-suffix): ex_conv.c $(ex_header_deps)
ex_utils.$(obj-suffix): ex_utils.c $(ex_header_deps)
exclos.$(obj-suffix): exclos.c $(ex_header_deps)
excn2s.$(obj-suffix): excn2s.c $(ex_header_deps)
excopy.$(obj-suffix): excopy.c $(ex_header_deps)
excre.$(obj-suffix): excre.c $(ex_header_deps)
exerr.$(obj-suffix): exerr.c $(ex_header_deps)
exgatm.$(obj-suffix): exgatm.c $(ex_header_deps)
exgatn.$(obj-suffix): exgatn.c $(ex_header_deps)
exgatt.$(obj-suffix): exgatt.c $(ex_header_deps)
exgattp.$(obj-suffix): exgattp.c $(ex_header_deps)
exgblk.$(obj-suffix): exgblk.c $(ex_header_deps)
exgcns.$(obj-suffix): exgcns.c $(ex_header_deps)
exgcon.$(obj-suffix): exgcon.c $(ex_header_deps)
exgconn.$(obj-suffix): exgconn.c $(ex_header_deps)
exgcor.$(obj-suffix): exgcor.c $(ex_header_deps)
exgcset.$(obj-suffix): exgcset.c $(ex_header_deps)
exgcss.$(obj-suffix): exgcss.c $(ex_header_deps)
exgcssc.$(obj-suffix): exgcssc.c $(ex_header_deps)
exgean.$(obj-suffix): exgean.c $(ex_header_deps)
exgeat.$(obj-suffix): exgeat.c $(ex_header_deps)
exgebi.$(obj-suffix): exgebi.c $(ex_header_deps)
exgecpp.$(obj-suffix): exgecpp.c $(ex_header_deps)
exgelb.$(obj-suffix): exgelb.c $(ex_header_deps)
exgelc.$(obj-suffix): exgelc.c $(ex_header_deps)
exgem.$(obj-suffix): exgem.c $(ex_header_deps)
exgenm.$(obj-suffix): exgenm.c $(ex_header_deps)
exgev.$(obj-suffix): exgev.c $(ex_header_deps)
exgevid.$(obj-suffix): exgevid.c $(ex_header_deps)
exgevt.$(obj-suffix): exgevt.c $(ex_header_deps)
exgfrm.$(obj-suffix): exgfrm.c $(ex_header_deps)
exggv.$(obj-suffix): exggv.c $(ex_header_deps)
exggvt.$(obj-suffix): exggvt.c $(ex_header_deps)
exgidm.$(obj-suffix): exgidm.c $(ex_header_deps)
exgids.$(obj-suffix): exgids.c $(ex_header_deps)
exginf.$(obj-suffix): exginf.c $(ex_header_deps)
exgini.$(obj-suffix): exgini.c $(ex_header_deps)
exginix.$(obj-suffix): exginix.c $(ex_header_deps)
exgmap.$(obj-suffix): exgmap.c $(ex_header_deps)
exgmp.$(obj-suffix): exgmp.c $(ex_header_deps)
exgnam.$(obj-suffix): exgnam.c $(ex_header_deps)
exgnams.$(obj-suffix): exgnams.c $(ex_header_deps)
exgnconn.$(obj-suffix): exgnconn.c $(ex_header_deps)
exgncor.$(obj-suffix): exgncor.c $(ex_header_deps)
exgnm.$(obj-suffix): exgnm.c $(ex_header_deps)
exgnmap.$(obj-suffix): exgnmap.c $(ex_header_deps)
exgnnm.$(obj-suffix): exgnnm.c $(ex_header_deps)
exgnnv.$(obj-suffix): exgnnv.c $(ex_header_deps)
exgnoatt.$(obj-suffix): exgnoatt.c $(ex_header_deps)
exgnp.$(obj-suffix): exgnp.c $(ex_header_deps)
exgns.$(obj-suffix): exgns.c $(ex_header_deps)
exgnsd.$(obj-suffix): exgnsd.c $(ex_header_deps)
exgnsi.$(obj-suffix): exgnsi.c $(ex_header_deps)
exgnstt.$(obj-suffix): exgnstt.c $(ex_header_deps)
exgnsv.$(obj-suffix): exgnsv.c $(ex_header_deps)
exgnsvid.$(obj-suffix): exgnsvid.c $(ex_header_deps)
exgnv.$(obj-suffix): exgnv.c $(ex_header_deps)
exgnvar.$(obj-suffix): exgnvar.c $(ex_header_deps)
exgnvid.$(obj-suffix): exgnvid.c $(ex_header_deps)
exgnvt.$(obj-suffix): exgnvt.c $(ex_header_deps)
exgnvv.$(obj-suffix): exgnvv.c $(ex_header_deps)
exgoatt.$(obj-suffix): exgoatt.c $(ex_header_deps)
exgoea.$(obj-suffix): exgoea.c $(ex_header_deps)
exgotv.$(obj-suffix): exgotv.c $(ex_header_deps)
exgp.$(obj-suffix): exgp.c $(ex_header_deps)
exgpa.$(obj-suffix): exgpa.c $(ex_header_deps)
exgpem.$(obj-suffix): exgpem.c $(ex_header_deps)
exgpn.$(obj-suffix): exgpn.c $(ex_header_deps)
exgqa.$(obj-suffix): exgqa.c $(ex_header_deps)
exgset.$(obj-suffix): exgset.c $(ex_header_deps)
exgsetd.$(obj-suffix): exgsetd.c $(ex_header_deps)
exgsetp.$(obj-suffix): exgsetp.c $(ex_header_deps)
exgsnl.$(obj-suffix): exgsnl.c $(ex_header_deps)
exgsp.$(obj-suffix): exgsp.c $(ex_header_deps)
exgss.$(obj-suffix): exgss.c $(ex_header_deps)
exgssc.$(obj-suffix): exgssc.c $(ex_header_deps)
exgssd.$(obj-suffix): exgssd.c $(ex_header_deps)
exgssi.$(obj-suffix): exgssi.c $(ex_header_deps)
exgssn.$(obj-suffix): exgssn.c $(ex_header_deps)
exgsstt.$(obj-suffix): exgsstt.c $(ex_header_deps)
exgssv.$(obj-suffix): exgssv.c $(ex_header_deps)
exgssvid.$(obj-suffix): exgssvid.c $(ex_header_deps)
exgtim.$(obj-suffix): exgtim.c $(ex_header_deps)
exgtt.$(obj-suffix): exgtt.c $(ex_header_deps)
exgvan.$(obj-suffix): exgvan.c $(ex_header_deps)
exgvar.$(obj-suffix): exgvar.c $(ex_header_deps)
exgvarnam.$(obj-suffix): exgvarnam.c $(ex_header_deps)
exgvarnams.$(obj-suffix): exgvarnams.c $(ex_header_deps)
exgvarparam.$(obj-suffix): exgvarparam.c $(ex_header_deps)
exgvart.$(obj-suffix): exgvart.c $(ex_header_deps)
exgvartab.$(obj-suffix): exgvartab.c $(ex_header_deps)
exgvid.$(obj-suffix): exgvid.c $(ex_header_deps)
exgvnm.$(obj-suffix): exgvnm.c $(ex_header_deps)
exgvp.$(obj-suffix): exgvp.c $(ex_header_deps)
exgvtt.$(obj-suffix): exgvtt.c $(ex_header_deps)
exgvv.$(obj-suffix): exgvv.c $(ex_header_deps)
exinq.$(obj-suffix): exinq.c $(ex_header_deps)
exopen.$(obj-suffix): exopen.c $(ex_header_deps)
exopts.$(obj-suffix): exopts.c $(ex_header_deps)
expatn.$(obj-suffix): expatn.c $(ex_header_deps)
expatt.$(obj-suffix): expatt.c $(ex_header_deps)
expattp.$(obj-suffix): expattp.c $(ex_header_deps)
expblk.$(obj-suffix): expblk.c $(ex_header_deps)
expcab.$(obj-suffix): expcab.c $(ex_header_deps)
expclb.$(obj-suffix): expclb.c $(ex_header_deps)
expcns.$(obj-suffix): expcns.c $(ex_header_deps)
expcon.$(obj-suffix): expcon.c $(ex_header_deps)
expconn.$(obj-suffix): expconn.c $(ex_header_deps)
expcor.$(obj-suffix): expcor.c $(ex_header_deps)
expcset.$(obj-suffix): expcset.c $(ex_header_deps)
expcss.$(obj-suffix): expcss.c $(ex_header_deps)
expean.$(obj-suffix): expean.c $(ex_header_deps)
expeat.$(obj-suffix): expeat.c $(ex_header_deps)
expecpp.$(obj-suffix): expecpp.c $(ex_header_deps)
expelb.$(obj-suffix): expelb.c $(ex_header_deps)
expelc.$(obj-suffix): expelc.c $(ex_header_deps)
expem.$(obj-suffix): expem.c $(ex_header_deps)
expenm.$(obj-suffix): expenm.c $(ex_header_deps)
expev.$(obj-suffix): expev.c $(ex_header_deps)
expfrm.$(obj-suffix): expfrm.c $(ex_header_deps)
expgv.$(obj-suffix): expgv.c $(ex_header_deps)
expidm.$(obj-suffix): expidm.c $(ex_header_deps)
expinf.$(obj-suffix): expinf.c $(ex_header_deps)
expini.$(obj-suffix): expini.c $(ex_header_deps)
expinix.$(obj-suffix): expinix.c $(ex_header_deps)
expmap.$(obj-suffix): expmap.c $(ex_header_deps)
expmp.$(obj-suffix): expmp.c $(ex_header_deps)
expnam.$(obj-suffix): expnam.c $(ex_header_deps)
expnams.$(obj-suffix): expnams.c $(ex_header_deps)
expncor.$(obj-suffix): expncor.c $(ex_header_deps)
expnm.$(obj-suffix): expnm.c $(ex_header_deps)
expnmap.$(obj-suffix): expnmap.c $(ex_header_deps)
expnnm.$(obj-suffix): expnnm.c $(ex_header_deps)
expnnv.$(obj-suffix): expnnv.c $(ex_header_deps)
expnoatt.$(obj-suffix): expnoatt.c $(ex_header_deps)
expnp.$(obj-suffix): expnp.c $(ex_header_deps)
expns.$(obj-suffix): expns.c $(ex_header_deps)
expnsd.$(obj-suffix): expnsd.c $(ex_header_deps)
expnstt.$(obj-suffix): expnstt.c $(ex_header_deps)
expnsv.$(obj-suffix): expnsv.c $(ex_header_deps)
expnv.$(obj-suffix): expnv.c $(ex_header_deps)
expnvar.$(obj-suffix): expnvar.c $(ex_header_deps)
expnvv.$(obj-suffix): expnvv.c $(ex_header_deps)
expoatt.$(obj-suffix): expoatt.c $(ex_header_deps)
expoea.$(obj-suffix): expoea.c $(ex_header_deps)
expp.$(obj-suffix): expp.c $(ex_header_deps)
exppa.$(obj-suffix): exppa.c $(ex_header_deps)
exppem.$(obj-suffix): exppem.c $(ex_header_deps)
exppn.$(obj-suffix): exppn.c $(ex_header_deps)
exppsetd.$(obj-suffix): exppsetd.c $(ex_header_deps)
expqa.$(obj-suffix): expqa.c $(ex_header_deps)
expset.$(obj-suffix): expset.c $(ex_header_deps)
expsetd.$(obj-suffix): expsetd.c $(ex_header_deps)
expsetp.$(obj-suffix): expsetp.c $(ex_header_deps)
expsp.$(obj-suffix): expsp.c $(ex_header_deps)
expss.$(obj-suffix): expss.c $(ex_header_deps)
expssd.$(obj-suffix): expssd.c $(ex_header_deps)
expsstt.$(obj-suffix): expsstt.c $(ex_header_deps)
expssv.$(obj-suffix): expssv.c $(ex_header_deps)
exptim.$(obj-suffix): exptim.c $(ex_header_deps)
exptt.$(obj-suffix): exptt.c $(ex_header_deps)
expvan.$(obj-suffix): expvan.c $(ex_header_deps)
expvar.$(obj-suffix): expvar.c $(ex_header_deps)
expvarnam.$(obj-suffix): expvarnam.c $(ex_header_deps)
expvarnams.$(obj-suffix): expvarnams.c $(ex_header_deps)
expvarparam.$(obj-suffix): expvarparam.c $(ex_header_deps)
expvartab.$(obj-suffix): expvartab.c $(ex_header_deps)
expvnm.$(obj-suffix): expvnm.c $(ex_header_deps)
expvp.$(obj-suffix): expvp.c $(ex_header_deps)
expvpa.$(obj-suffix): expvpa.c $(ex_header_deps)
expvpax.$(obj-suffix): expvpax.c $(ex_header_deps)
expvpc.$(obj-suffix): expvpc.c $(ex_header_deps)
expvtt.$(obj-suffix): expvtt.c $(ex_header_deps)
expvv.$(obj-suffix): expvv.c $(ex_header_deps)
exupda.$(obj-suffix): exupda.c $(ex_header_deps)

