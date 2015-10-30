include ../../Make.common

INCLUDES = -I.

OBJS = hilbert.o morton.o cmp.o

%.o : %.c
	@echo "Compiling "$<"..."
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@


target := ../lib/$(HOSTTYPE)/libsfcurves.a

ifeq ($(enable-shared),yes)
  target := ../lib/$(HOSTTYPE)/libsfcurves.so
endif


all:: $(target)

../lib/$(HOSTTYPE)/libsfcurves.a: $(OBJS)
	@echo "Linking "$@
	$(shell mkdir -p ../lib/$(HOSTTYPE))
	@$(AR) rv $@ $(OBJS)

../lib/$(HOSTTYPE)/libsfcurves.so: $(OBJS)
	@echo "Linking "$@
	$(shell mkdir -p ../lib/$(HOSTTYPE))
	@$(CC) -shared -o $@ $(OBJS)

clean:
	rm -f *.o *~

clobber:
	make clean
	rm -f ../lib/$(HOSTTYPE)/libsfcurves.*
