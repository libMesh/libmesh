include ../../Make.common

INCLUDES = -I.

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))



target := ../lib/$(hosttype)_$(METHOD)/libsfcurves.a

ifeq ($(enable-shared),yes)
  target := ../lib/$(hosttype)_$(METHOD)/libsfcurves.so
endif


all:: $(target)

../lib/$(hosttype)_$(METHOD)/libsfcurves.a: $(OBJS)
	@echo "Linking "$@
	$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../lib/$(hosttype)_$(METHOD)/libsfcurves.so: $(OBJS)
	@echo "Linking "$@
	$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(CC) -shared -o $@ $(OBJS)

clean:
	@rm -f $(OBJS) *~

clobber:
	@rm -f *.o *.go *.pgo *.syo
	@rm -f ../lib/$(hosttype)_$(METHOD)/libsfcurves.*

distclean:
	@$(MAKE) clobber
	@rm -f ../lib/*/libsfcurves.*