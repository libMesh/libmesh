include ../../Make.common

INCLUDES = -I.

OBJS = hilbert.o morton.o cmp.o

%.o : %.c
	@echo "Compiling "$<"..."
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@


target := ../lib/$(hosttype)/libsfcurves.a

ifeq ($(enable-shared),yes)
  target := ../lib/$(hosttype)/libsfcurves.so
endif


all:: $(target)

../lib/$(hosttype)/libsfcurves.a: $(OBJS)
	@echo "Linking "$@
	$(shell mkdir -p ../lib/$(hosttype))
	@$(AR) rv $@ $(OBJS)

../lib/$(hosttype)/libsfcurves.so: $(OBJS)
	@echo "Linking "$@
	$(shell mkdir -p ../lib/$(hosttype))
	@$(CC) -shared -o $@ $(OBJS)

clean:
	@rm -f *.o *~

clobber:
	@$(MAKE) clean
	@rm -f ../lib/$(hosttype)/libsfcurves.*
