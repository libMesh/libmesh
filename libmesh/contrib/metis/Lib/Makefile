include ../../../Make.common

INCLUDE = -I.

LD = $(CC) -L. 


# OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o memory.o \
#        pmetis.o pqueue.o refine.o util.o timing.o debug.o \
#        bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
#        kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
#        mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
#        myqsort.o compress.o parmetis.o estmem.o \
#        mpmetis.o mcoarsen.o mmatch.o minitpart.o mbalance.o \
#        mrefine.o mutil.o mfm.o mkmetis.o mkwayrefine.o mkwayfmh.o \
#        mrefine2.o minitpart2.o mbalance2.o mfm2.o \
#        kvmetis.o kwayvolrefine.o kwayvolfm.o subdomains.o


SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


target := ../../lib/$(hosttype)_$(METHOD)/libmetis.a

ifeq ($(enable-shared),yes)
  target := ../../lib/$(hosttype)_$(METHOD)/libmetis.so
endif


all:: $(target)

../../lib/$(hosttype)_$(METHOD)/libmetis.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../../lib/$(hosttype)_$(METHOD)/libmetis.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(CC) -shared -o $@ $(OBJS)

clean:
	@rm -f $(OBJS) *~

clobber:
	@rm -f *.o *.go *.pgo *.syo
	@rm -f ../../lib/$(hosttype)_$(METHOD)/libmetis.*

distclean:
	@$(MAKE) clobber
	@rm -f ../../lib/*/libmetis.*
