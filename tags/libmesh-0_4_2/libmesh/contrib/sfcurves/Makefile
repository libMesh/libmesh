include ../../Make.common

INCLUDES = -I.

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))



target := ../lib/$(hosttype)_$(METHOD)/libsfcurves.a

ifeq ($(enable-shared),yes)
  target := ../lib/$(hosttype)_$(METHOD)/libsfcurves.so
endif


ifeq ($(enable-sfcurves),no)
  all::
	@echo "<<< sfcurves is disabled, not building >>>"
else
  all:: $(target)
endif

../lib/$(hosttype)_$(METHOD)/libsfcurves.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../lib/$(hosttype)_$(METHOD)/libsfcurves.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)


.PHONY: clean clobber distclean



clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f ../lib/$(hosttype)_$(METHOD)/libsfcurves.*

distclean:
	@$(MAKE) clobber
	@rm -f ../lib/*/libsfcurves.*