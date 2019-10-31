check_PROGRAMS  =
check_PROGRAMS += message_tag_unit
check_PROGRAMS += packed_range_unit
check_PROGRAMS += parallel_sync_unit
check_PROGRAMS += parallel_unit

AM_CPPFLAGS  =
AM_CPPFLAGS += -I$(top_srcdir)/src/algorithms/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parallel/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #passmess_version.h

AM_LDFLAGS = $(top_builddir)/src/libpassmess.la

# Sources for these tests
message_tag_unit_SOURCES = message_tag_unit.C
packed_range_unit_SOURCES = packed_range_unit.C
parallel_sync_unit_SOURCES = parallel_sync_unit.C
parallel_unit_SOURCES = parallel_unit.C

#Define tests to actually be run
TESTS  =
TESTS += message_tag_unit
TESTS += packed_range_unit
TESTS += parallel_sync_unit
TESTS += parallel_unit

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
