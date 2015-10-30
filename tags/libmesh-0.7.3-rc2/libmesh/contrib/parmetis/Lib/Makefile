include ../../../Make.common

libmesh_INCLUDE += -I.


SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := $(EXTERNAL_LIBDIR)/libparmetis$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/libparmetis$(shared_libext)
endif


ifeq ($(enable-parmetis),no)
  all::
	@echo "<<< parmetis is disabled, not building >>>"
else
    ifeq (x$(MPI_IMPL),x)
      all::
	@echo "<<< no MPI found, not building parmetis >>>"
    else
      all:: $(target)
    endif
endif

ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/libparmetis$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/libparmetis$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

$(EXTERNAL_LIBDIR)/libparmetis$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(libmesh_CC) $(libmesh_CSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)

clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/libparmetis.*

distclean:
	@$(MAKE) clobber
	@rm -f ../../lib/*/libparmetis.*

echo:
	@echo $(libmesh_INCLUDE)


# Dependencies for parmetis files.  Every c file depends on the following headers...
parmetis_header_deps := defs.h macros.h parmetis.h proto.h rename.h stdheaders.h struct.h parmetislib.h

# Script to generate the following lines automatically, unfortunately
# I don't think there is a nice Makfile way to specify all of these in 
# one line?
# for i in *.c; do echo "`basename $i c`\$(obj-suffix): $i \$(parmetis_header_deps)"; done

akwayfm.$(obj-suffix): akwayfm.c $(parmetis_header_deps)
ametis.$(obj-suffix): ametis.c $(parmetis_header_deps)
backcompat.$(obj-suffix): backcompat.c $(parmetis_header_deps)
balancemylink.$(obj-suffix): balancemylink.c $(parmetis_header_deps)
coarsen.$(obj-suffix): coarsen.c $(parmetis_header_deps)
comm.$(obj-suffix): comm.c $(parmetis_header_deps)
csrmatch.$(obj-suffix): csrmatch.c $(parmetis_header_deps)
debug.$(obj-suffix): debug.c $(parmetis_header_deps)
diffutil.$(obj-suffix): diffutil.c $(parmetis_header_deps)
fpqueue.$(obj-suffix): fpqueue.c $(parmetis_header_deps)
frename.$(obj-suffix): frename.c $(parmetis_header_deps)
gkmetis.$(obj-suffix): gkmetis.c $(parmetis_header_deps)
grsetup.$(obj-suffix): grsetup.c $(parmetis_header_deps)
iidxsort.$(obj-suffix): iidxsort.c $(parmetis_header_deps)
iintsort.$(obj-suffix): iintsort.c $(parmetis_header_deps)
ikeysort.$(obj-suffix): ikeysort.c $(parmetis_header_deps)
ikeyvalsort.$(obj-suffix): ikeyvalsort.c $(parmetis_header_deps)
initbalance.$(obj-suffix): initbalance.c $(parmetis_header_deps)
initmsection.$(obj-suffix): initmsection.c $(parmetis_header_deps)
initpart.$(obj-suffix): initpart.c $(parmetis_header_deps)
kmetis.$(obj-suffix): kmetis.c $(parmetis_header_deps)
kwaybalance.$(obj-suffix): kwaybalance.c $(parmetis_header_deps)
kwayfm.$(obj-suffix): kwayfm.c $(parmetis_header_deps)
kwayrefine.$(obj-suffix): kwayrefine.c $(parmetis_header_deps)
lmatch.$(obj-suffix): lmatch.c $(parmetis_header_deps)
match.$(obj-suffix): match.c $(parmetis_header_deps)
mdiffusion.$(obj-suffix): mdiffusion.c $(parmetis_header_deps)
memory.$(obj-suffix): memory.c $(parmetis_header_deps)
mesh.$(obj-suffix): mesh.c $(parmetis_header_deps)
mmetis.$(obj-suffix): mmetis.c $(parmetis_header_deps)
move.$(obj-suffix): move.c $(parmetis_header_deps)
node_refine.$(obj-suffix): node_refine.c $(parmetis_header_deps)
ometis.$(obj-suffix): ometis.c $(parmetis_header_deps)
order.$(obj-suffix): order.c $(parmetis_header_deps)
pspases.$(obj-suffix): pspases.c $(parmetis_header_deps)
redomylink.$(obj-suffix): redomylink.c $(parmetis_header_deps)
remap.$(obj-suffix): remap.c $(parmetis_header_deps)
rmetis.$(obj-suffix): rmetis.c $(parmetis_header_deps)
selectq.$(obj-suffix): selectq.c $(parmetis_header_deps)
serial.$(obj-suffix): serial.c $(parmetis_header_deps)
setup.$(obj-suffix): setup.c $(parmetis_header_deps)
stat.$(obj-suffix): stat.c $(parmetis_header_deps)
timer.$(obj-suffix): timer.c $(parmetis_header_deps)
util.$(obj-suffix): util.c $(parmetis_header_deps)
wave.$(obj-suffix): wave.c $(parmetis_header_deps)
weird.$(obj-suffix): weird.c $(parmetis_header_deps)
xyzpart.$(obj-suffix): xyzpart.c $(parmetis_header_deps)
