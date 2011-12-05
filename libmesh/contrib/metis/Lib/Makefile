include ../../../Make.common

INCLUDE = -I. 

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


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
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/libmetis.*

distclean:
	@$(MAKE) clobber
	@rm -f ../../lib/*/libmetis.*


# Dependencies for metis files.  Every c file depends on the following headers...
metis_header_deps := defs.h macros.h proto.h rename.h stdheaders.h struct.h metis.h

# Script to generate the following lines automatically, unfortunately
# I don't think there is a nice Makfile way to specify all of these in 
# one line?
# for i in *.c; do echo "`basename $i c`\$(obj-suffix): $i \$(metis_header_deps)"; done

balance.$(obj-suffix): balance.c $(metis_header_deps)
bucketsort.$(obj-suffix): bucketsort.c $(metis_header_deps)
ccgraph.$(obj-suffix): ccgraph.c $(metis_header_deps)
checkgraph.$(obj-suffix): checkgraph.c $(metis_header_deps)
coarsen.$(obj-suffix): coarsen.c $(metis_header_deps)
compress.$(obj-suffix): compress.c $(metis_header_deps)
debug.$(obj-suffix): debug.c $(metis_header_deps)
estmem.$(obj-suffix): estmem.c $(metis_header_deps)
fm.$(obj-suffix): fm.c $(metis_header_deps)
fortran.$(obj-suffix): fortran.c $(metis_header_deps)
frename.$(obj-suffix): frename.c $(metis_header_deps)
graph.$(obj-suffix): graph.c $(metis_header_deps)
initpart.$(obj-suffix): initpart.c $(metis_header_deps)
kmetis.$(obj-suffix): kmetis.c $(metis_header_deps)
kvmetis.$(obj-suffix): kvmetis.c $(metis_header_deps)
kwayfm.$(obj-suffix): kwayfm.c $(metis_header_deps)
kwayrefine.$(obj-suffix): kwayrefine.c $(metis_header_deps)
kwayvolfm.$(obj-suffix): kwayvolfm.c $(metis_header_deps)
kwayvolrefine.$(obj-suffix): kwayvolrefine.c $(metis_header_deps)
match.$(obj-suffix): match.c $(metis_header_deps)
mbalance.$(obj-suffix): mbalance.c $(metis_header_deps)
mbalance2.$(obj-suffix): mbalance2.c $(metis_header_deps)
mcoarsen.$(obj-suffix): mcoarsen.c $(metis_header_deps)
memory.$(obj-suffix): memory.c $(metis_header_deps)
mesh.$(obj-suffix): mesh.c $(metis_header_deps)
meshpart.$(obj-suffix): meshpart.c $(metis_header_deps)
mfm.$(obj-suffix): mfm.c $(metis_header_deps)
mfm2.$(obj-suffix): mfm2.c $(metis_header_deps)
mincover.$(obj-suffix): mincover.c $(metis_header_deps)
minitpart.$(obj-suffix): minitpart.c $(metis_header_deps)
minitpart2.$(obj-suffix): minitpart2.c $(metis_header_deps)
mkmetis.$(obj-suffix): mkmetis.c $(metis_header_deps)
mkwayfmh.$(obj-suffix): mkwayfmh.c $(metis_header_deps)
mkwayrefine.$(obj-suffix): mkwayrefine.c $(metis_header_deps)
mmatch.$(obj-suffix): mmatch.c $(metis_header_deps)
mmd.$(obj-suffix): mmd.c $(metis_header_deps)
mpmetis.$(obj-suffix): mpmetis.c $(metis_header_deps)
mrefine.$(obj-suffix): mrefine.c $(metis_header_deps)
mrefine2.$(obj-suffix): mrefine2.c $(metis_header_deps)
mutil.$(obj-suffix): mutil.c $(metis_header_deps)
myqsort.$(obj-suffix): myqsort.c $(metis_header_deps)
ometis.$(obj-suffix): ometis.c $(metis_header_deps)
parmetis.$(obj-suffix): parmetis.c $(metis_header_deps)
pmetis.$(obj-suffix): pmetis.c $(metis_header_deps)
pqueue.$(obj-suffix): pqueue.c $(metis_header_deps)
refine.$(obj-suffix): refine.c $(metis_header_deps)
separator.$(obj-suffix): separator.c $(metis_header_deps)
sfm.$(obj-suffix): sfm.c $(metis_header_deps)
srefine.$(obj-suffix): srefine.c $(metis_header_deps)
stat.$(obj-suffix): stat.c $(metis_header_deps)
subdomains.$(obj-suffix): subdomains.c $(metis_header_deps)
timing.$(obj-suffix): timing.c $(metis_header_deps)
util.$(obj-suffix): util.c $(metis_header_deps)



