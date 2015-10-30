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

# Dependencies for Laspack source files
eigenval.$(obj-suffix):\
    eigenval.c\
    copyrght.h\
    eigenval.h\
    elcmp.h\
    errhandl.h\
    itersolv.h\
    laspack_config.h\
    lastypes.h\
    operats.h\
    precond.h\
    qmatrix.h\
    qvector.h\
    rtc.h
errhandl.$(obj-suffix):\
    errhandl.c\
    copyrght.h\
    errhandl.h\
    laspack_config.h
factor.$(obj-suffix):\
    factor.c\
    copyrght.h\
    elcmp.h\
    errhandl.h\
    factor.h\
    laspack_config.h\
    lastypes.h\
    qmatrix.h\
    qvector.h
itersolv.$(obj-suffix):\
    itersolv.c\
    copyrght.h\
    eigenval.h\
    elcmp.h\
    errhandl.h\
    itersolv.h\
    laspack_config.h\
    lastypes.h\
    operats.h\
    precond.h\
    qmatrix.h\
    qvector.h\
    rtc.h
matrix.$(obj-suffix):\
    matrix.c\
    copyrght.h\
    elcmp.h\
    errhandl.h\
    laspack_config.h\
    lastypes.h\
    matrix.h
mlsolv.$(obj-suffix):\
    mlsolv.c\
    copyrght.h\
    eigenval.h\
    elcmp.h\
    errhandl.h\
    itersolv.h\
    laspack_config.h\
    lastypes.h\
    mlsolv.h\
    operats.h\
    precond.h\
    qmatrix.h\
    qvector.h\
    rtc.h
operats.$(obj-suffix):\
    operats.c\
    copyrght.h\
    elcmp.h\
    errhandl.h\
    laspack_config.h\
    lastypes.h\
    operats.h\
    qmatrix.h\
    qvector.h
precond.$(obj-suffix):\
    precond.c\
    copyrght.h\
    eigenval.h\
    elcmp.h\
    errhandl.h\
    factor.h\
    itersolv.h\
    laspack_config.h\
    lastypes.h\
    operats.h\
    precond.h\
    qmatrix.h\
    qvector.h
qmatrix.$(obj-suffix):\
    qmatrix.c\
    copyrght.h\
    elcmp.h\
    errhandl.h\
    laspack_config.h\
    lastypes.h\
    qmatrix.h\
    qvector.h
qvector.$(obj-suffix):\
    qvector.c\
    copyrght.h\
    elcmp.h\
    errhandl.h\
    laspack_config.h\
    lastypes.h\
    qvector.h
rtc.$(obj-suffix):\
    rtc.c\
    copyrght.h\
    eigenval.h\
    elcmp.h\
    errhandl.h\
    itersolv.h\
    laspack_config.h\
    lastypes.h\
    operats.h\
    precond.h\
    qmatrix.h\
    qvector.h\
    rtc.h
