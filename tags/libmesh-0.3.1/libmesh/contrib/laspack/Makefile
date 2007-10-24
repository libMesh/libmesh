include ../../Make.common

INCLUDES = -I.

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))



target := ../lib/$(hosttype)_$(METHOD)/liblaspack.a

ifeq ($(enable-shared),yes)
  target := ../lib/$(hosttype)_$(METHOD)/liblaspack.so
endif


all:: $(target)

../lib/$(hosttype)_$(METHOD)/liblaspack.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../lib/$(hosttype)_$(METHOD)/liblaspack.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS)


.PHONY: clean clobber distclean



clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f ../lib/$(hosttype)_$(METHOD)/liblaspack.*

distclean:
	@$(MAKE) clobber
	@rm -f ../lib/*/liblaspack.*