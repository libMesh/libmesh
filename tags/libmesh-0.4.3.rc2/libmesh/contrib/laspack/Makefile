#
# let Make.common know that we are the Laspack
# Makefile, so that it can decide which compilation
# rules to tell
#
included-from-laspack := yes

include ../../Make.common

INCLUDES = -I.

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))



target := ../lib/$(hosttype)_$(METHOD)/liblaspack.a

ifeq ($(enable-shared),yes)
  target := ../lib/$(hosttype)_$(METHOD)/liblaspack.so
endif


ifeq ($(enable-laspack),no)
  all::
	@echo "<<< laspack is disabled, not building >>>"
else
  all:: $(target)
endif

../lib/$(hosttype)_$(METHOD)/liblaspack.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

ifeq ($(enable-complex),yes)
#
# Complex arithmetic: use C++ 
#
../lib/$(hosttype)_$(METHOD)/liblaspack.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)

else
#
# Real arithmetic: use C
#
../lib/$(hosttype)_$(METHOD)/liblaspack.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)

endif


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
