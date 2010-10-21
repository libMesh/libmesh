#
# set appropriate the following variables in your environment (if necessary):
#
# CC		for the C compiler
# CPP		for the C++ compiler
# CFLAGS	for options of the C compiler
# CPPFLAGS	for options of the C++ compiler
# PFLAGS	for options of the Pascal compiler
# FFLAGS	for options of the FORTRAN compiler
# LDFLAGS	for linker options
#

#
# ARCH_EXT can be used in order to install libraries in different directories
# depending on the computer architecture,
# e.g. $(HOME)/lib/sunos for ARCH_EXT = '/sunos'
#
#ARCH_EXT	=

#
# set the path for the root of the include directories here,
# e.g. /usr/local/include
#
INCROOT		= ..
#
# set the destination directories for the library and include files
#
LIBDEST		= $(HOME)/lib$(ARCH_EXT)
INCDEST		= $(HOME)/include

#
# the following text was created automaticaly. You should change it carefully.
#

SHELL		= /bin/sh

LIBNAME		= laspack

LIBRARY		= lib$(LIBNAME).a

HDRS		= copyrght.h \
		eigenval.h \
		elcmp.h \
		errhandl.h \
		factor.h \
		itersolv.h \
		lastypes.h \
		matrix.h \
		mlsolv.h \
		operats.h \
		precond.h \
		qmatrix.h \
		rtc.h \
		qvector.h \
		version.h

EXTHDRS		=

SRCS		= eigenval.c \
		errhandl.c \
		factor.c \
		itersolv.c \
		matrix.c \
		mlsolv.c \
		operats.c \
		precond.c \
		qmatrix.c \
		rtc.c \
		vector.c

OBJS		= eigenval.o \
		errhandl.o \
		factor.o \
		itersolv.o \
		matrix.o \
		mlsolv.o \
		operats.o \
		precond.o \
		qmatrix.o \
		rtc.o \
		vector.o

LIBS		=

COMPFLAGS	=  

# compiler options passed throuth enviroment variables
#CFLAGS		=
#PFLAGS		=
#FFLAGS		=
#CXXFLAGS	=

LIBLOCAL	= /usr/local/lib
INCLOCAL	= /usr/local/include

INSTALL		= mv

ARFLAGS		= cru

LINTLIBS	=

LINTFLAGS	= -u -I$(INCROOT) $(CFLAGS)

MAKEFILE	= makefile

PRINT		= pr

PRINTFLAGS	=

LP		= lp

LPFLAGS		= 

all:		$(LIBRARY)

clean:;		@rm -rf $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags
		@if [ -f compllist ]; then rm -f compllist; fi
		@if [ -f cleanlist ]; then rm -f cleanlist; fi
		@find . -type f -print > compllist
		@sed -n \
			-e "/~/ w cleanlist" \
			-e '/%/ w cleanlist' \
			-e '/.bak/ w cleanlist' \
			-e '/.obj/ w cleanlist' \
			-e '/.exe/ w cleanlist' \
			-e '/.aux/ w cleanlist' \
			-e '/.blg/ w cleanlist' \
			-e '/.dvi/ w cleanlist' \
			-e '/.glo/ w cleanlist' \
			-e '/.idx/ w cleanlist' \
			-e '/.ilg/ w cleanlist' \
			-e '/.ind/ w cleanlist' \
			-e '/.lof/ w cleanlist' \
			-e '/.log/ w cleanlist' \
			-e '/.lot/ w cleanlist' \
			-e '/.toc/ w cleanlist' \
			compllist
		@rm -f `cat cleanlist`
		@rm -f compllist
		@rm -f cleanlist

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(LIBDEST)
		@if [ $(LIBDEST) != . ]; then rm -f $(LIBDEST)/$(LIBRARY); fi
		@if [ $(LIBDEST) != . ]; then $(INSTALL) -f $(LIBRARY) $(LIBDEST); fi
		@echo Installing header files in $(INCDEST)/$(LIBNAME)
		@rm -rf $(INCDEST)/$(LIBNAME).old
		@if [ -d $(INCDEST)/$(LIBNAME) ]; then \
			mv $(INCDEST)/$(LIBNAME) $(INCDEST)/$(LIBNAME).old; \
		fi
		@mkdir $(INCDEST)/$(LIBNAME)
		@chmod 755 $(INCDEST)/$(LIBNAME)
		@cp *.h $(INCDEST)/$(LIBNAME)
		@chmod 644 $(INCDEST)/$(LIBNAME)/*

install-local:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(LIBLOCAL)
		@rm -f $(LIBLOCAL)/$(LIBRARY).old
		@if [ -f $(LIBLOCAL)/$(LIBRARY) ]; then \
			mv $(LIBLOCAL)/$(LIBRARY) $(LIBLOCAL)/$(LIBRARY).old; \
		fi
		@$(INSTALL) $(LIBRARY) $(LIBLOCAL)
		@chmod 755 $(LIBLOCAL)/$(LIBRARY)
		@echo Installing header files in $(INCLOCAL)/$(LIBNAME)
		@rm -rf $(INCLOCAL)/$(LIBNAME).old
		@if [ -d $(INCLOCAL)/$(LIBNAME) ]; then \
			mv $(INCLOCAL)/$(LIBNAME) $(INCLOCAL)/$(LIBNAME).old; \
		fi
		@mkdir $(INCLOCAL)/$(LIBNAME)
		@chmod 755 $(INCLOCAL)/$(LIBNAME)
		@cp *.h $(INCLOCAL)/$(LIBNAME)
		@chmod 755 $(INCLOCAL)/$(LIBNAME)/*

lint:		$(LINTLIBS) $(HDRS) $(EXTHDRS) $(SRCS)
		@$(LINT) $(LINTFLAGS) $(LINTLIBS) $(SRCS)

print:;		@$(PRINT) $(PRINTFLAGS) $(HDRS) $(SRCS) | $(LP) $(LPFLAGS)

tags:           $(HDRS) $(SRCS) 
		@ctags $(HDRS) $(SRCS)

touch:;		@touch $(HDRS) $(SRCS) $(MAKEFILE)

update:		$(LIBDEST)/$(LIBRARY)

d2u:;		@d2u $(HDRS) $(SRCS)

c:;		@$(MAKE) -f $(MAKEFILE) clean
cl:;		@$(MAKE) -f $(MAKEFILE) clobber
i:;             @$(MAKE) -f $(MAKEFILE) install
il:;		@$(MAKE) -f $(MAKEFILE) install-local
l:;		@$(MAKE) -f $(MAKEFILE) lint
t:;		@$(MAKE) -f $(MAKEFILE) touch
u:;		@$(MAKE) -f $(MAKEFILE) update 

$(LIBRARY):     $(OBJS) $(MAKEFILE)
		@echo "Loading $(LIBRARY) ..."
		@ar $(ARFLAGS) $(LIBRARY) $(OBJS)
		@ranlib $(LIBRARY)

$(LIBDEST)/$(LIBRARY):  $(HDRS) $(EXTHDRS) $(SRCS) $(LIBS) 
		@$(MAKE) -f $(MAKEFILE) install

.SUFFIXES: .c .cc  .p .f .o

.c.o:;		$(CC) -I$(INCROOT) $(CFLAGS) $(COMPFLAGS) -c $<
.cc.o:;		$(CPP) -I$(INCROOT) $(CPPFLAGS) $(COMPFLAGS) -c $<
.p.o:;		pc $(PFLAGS) $(COMPFLAGS) -c $<
.f.o:;		f77 $(FFLAGS) $(COMPFLAGS) -c $<
