include ../../../Make.common

# If enable-mpi is set, that means
# petsc has been disabled, and we need
# to include the mpi headers with an
# additional -I.
# Otherwise there will be one too many
# -I's.  (This is a hack)
#
# If both petsc and mpi are explicitly
# disabled, better not try to build anything!
ifeq ($(enable-mpi),yes)
  INCLUDE = -I$(MPI_INCLUDE) -I.
else
  INCLUDE = $(MPI_INCLUDE) -I.
endif


SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %.$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean



target := ../../lib/$(hosttype)_$(METHOD)/libparmetis.a

ifeq ($(enable-shared),yes)
  target := ../../lib/$(hosttype)_$(METHOD)/libparmetis.so
endif


ifeq ($(enable-mpi)_$(enable-petsc) , no_no)
  # Do nothing
  all::
	@echo "<<< Cannot build Parmetis without either petsc or mpi >>>"

else

all:: $(target)

../../lib/$(hosttype)_$(METHOD)/libparmetis.a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../../lib/$(hosttype)_$(METHOD)/libparmetis.so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)
endif

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
	@echo $(MPI_INCLUDE)