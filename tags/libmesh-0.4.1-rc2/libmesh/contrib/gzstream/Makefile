# ============================================================================
# gzstream, C++ iostream classes wrapping the zlib compression library.
# Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# ============================================================================
# 
# File          : Makefile
# Revision      : $Revision: 1.1 $
# Revision_date : $Date: 2003-09-16 00:41:29 $
# Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
# 
# ============================================================================

# ----------------------------------------------------------------------------
# adapt these settings to your need:
# add '-DGZSTREAM_NAMESPACE=name' to CPPFLAGS to place the classes
# in its own namespace. Note, this macro needs to be set while creating
# the library as well while compiling applications based on it.
# As an alternative, gzstream.C and gzstream.h can be edited.
# ----------------------------------------------------------------------------



# The preceeding text came from the authors of this package.  The
# following allows it to link properly with our libraries.
include ../../Make.common

INCLUDES = -I.
SRC     := $(wildcard *.C)
OBJS   	:= $(patsubst %.C, %.$(obj-suffix), $(SRC))

library_name := libgzstream



#
# The remaining lines will usually work with any simple
# contributed project.  For example sfcurves follows the
# same procedure.
#
target       := ../lib/$(hosttype)_$(METHOD)/$(library_name).a


ifeq ($(enable-shared),yes)	
  target := ../lib/$(hosttype)_$(METHOD)/$(library_name).so
endif

all:: $(target)

../lib/$(hosttype)_$(METHOD)/$(library_name).a: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(AR) rv $@ $(OBJS)

../lib/$(hosttype)_$(METHOD)/$(library_name).so: $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p ../lib/$(hosttype)_$(METHOD))
	@$(CC) $(CSHAREDFLAG) -o $@ $(OBJS) $(LDFLAGS)

.PHONY: clean clobber distclean



clean:
	@rm -f $(OBJS) *~

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f ../lib/$(hosttype)_$(METHOD)/$(library_name).*

distclean:
	@$(MAKE) clobber
	@rm -f ../lib/*/$(library_name).*

