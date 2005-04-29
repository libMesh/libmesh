include ../../Make.common

INCLUDE = -I.

# Package name and additional flags required for library compilation
pkgname  := triangle
CFLAGS   += -DTRILIBRARY

# (Relative) Path to location of shared libs
libpath := ../lib

# Only one source file for triangle
SRC     := triangle.c
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := $(libpath)/$(hosttype)_$(METHOD)/lib$(pkgname).a

ifeq ($(enable-shared),yes)
  target := $(libpath)/$(hosttype)_$(METHOD)/lib$(pkgname).so
endif

ifeq ($(enable-triangle),no)
  all::
	@echo "<<< $(pkgname) is disabled, not building >>>"
else
  all:: $(target)
endif

$(libpath)/$(hosttype)_$(METHOD)/lib$(pkgname).a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(libpath)/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

$(libpath)/$(hosttype)_$(METHOD)/lib$(pkgname).so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(libpath)/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)

clean:
	@rm -f $(OBJS) showme triangle tricall *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(libpath)/$(hosttype)_$(METHOD)/lib$(pkgname).*

distclean:
	@$(MAKE) clobber
	@rm -f $(libpath)/*/lib$(pkgname).*
