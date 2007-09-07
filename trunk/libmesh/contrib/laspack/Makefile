#
# let Make.common know that we are the Laspack
# Makefile, so that it can decide which compilation
# rules to tell
#
included-from-laspack := yes

include ../../Make.common

libmesh_INCLUDES = -I.

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))



target := $(EXTERNAL_LIBDIR)/liblaspack$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/liblaspack$(shared_libext)
endif


ifeq ($(enable-laspack),no)
  all::
	@echo "<<< laspack is disabled, not building >>>"
else
  all:: $(target)
endif

ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/liblaspack$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/liblaspack$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

ifeq ($(enable-complex),yes)
#
# Complex arithmetic: use C++ 
#
$(EXTERNAL_LIBDIR)/liblaspack$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(libmesh_CXX) $(libmesh_CXXSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)

else
#
# Real arithmetic: use C
#
$(EXTERNAL_LIBDIR)/liblaspack$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(libmesh_CC) $(libmesh_CSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)

endif


.PHONY: clean clobber distclean



clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f $(EXTERNAL_LIBDIR)/liblaspack.*

distclean:
	@$(MAKE) clobber
	@rm -f ../lib/*/liblaspack.*
