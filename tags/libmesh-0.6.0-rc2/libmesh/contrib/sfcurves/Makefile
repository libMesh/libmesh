include ../../Make.common

INCLUDES = -I.

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))



target := $(EXTERNAL_LIBDIR)/libsfcurves$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/libsfcurves$(shared_libext)
endif


ifeq ($(enable-sfcurves),no)
  all::
	@echo "<<< sfcurves is disabled, not building >>>"
else
  all:: $(target)
endif

ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/libsfcurves$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/libsfcurves$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

$(EXTERNAL_LIBDIR)/libsfcurves$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)


.PHONY: clean clobber distclean



clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/libsfcurves.*

distclean:
	@$(MAKE) clobber
	@rm -f ../lib/*/libsfcurves.*
