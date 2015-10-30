include ../../../Make.common

INCLUDE = $(MPI_INCLUDE) -I.


SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := ../../lib/$(hosttype)_$(METHOD)/libparmetis.a

ifeq ($(enable-shared),yes)
  target := ../../lib/$(hosttype)_$(METHOD)/libparmetis.so
endif


#ifeq (x$(MPI_IMPL),x)
ifeq ($(enable-parmetis),no)
  all::
	@echo "<<< parmetis is disabled, not building >>>"
else
  all:: $(target)
endif

../../lib/$(hosttype)_$(METHOD)/libparmetis.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../../lib/$(hosttype)_$(METHOD)/libparmetis.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)

clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f ../../lib/$(hosttype)_$(METHOD)/libparmetis.*

distclean:
	@$(MAKE) clobber
	@rm -f ../../lib/*/libparmetis.*

echo:
	@echo $(INCLUDE)