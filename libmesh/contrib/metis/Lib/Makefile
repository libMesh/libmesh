include ../../../Make.common

INCLUDE = -I. 

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := ../../lib/$(hosttype)_$(METHOD)/libmetis.a

ifeq ($(enable-shared),yes)
  target := ../../lib/$(hosttype)_$(METHOD)/libmetis.so
endif

ifeq ($(enable-metis),no)
  all::
	@echo "<<< metis is disabled, not building >>>"
else
  all:: $(target)
endif

../../lib/$(hosttype)_$(METHOD)/libmetis.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../../lib/$(hosttype)_$(METHOD)/libmetis.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)

clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f ../../lib/$(hosttype)_$(METHOD)/libmetis.*

distclean:
	@$(MAKE) clobber
	@rm -f ../../lib/*/libmetis.*
