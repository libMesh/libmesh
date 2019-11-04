
AM_CPPFLAGS  =
AM_CPPFLAGS += -I$(top_srcdir)/src/algorithms/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parallel/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #passmess_version.h

#Define tests to actually be run
TESTS  =

if BUILD_DBG_MODE
  message_tag_unit_dbg_SOURCES = message_tag_unit.C
  packed_range_unit_dbg_SOURCES = packed_range_unit.C
  parallel_sync_unit_dbg_SOURCES = parallel_sync_unit.C
  parallel_unit_dbg_SOURCES = parallel_unit.C
  message_tag_unit_dbg_LDFLAGS = $(top_builddir)/src/libpassmess_dbg.la
  packed_range_unit_dbg_LDFLAGS = $(top_builddir)/src/libpassmess_dbg.la
  parallel_sync_unit_dbg_LDFLAGS = $(top_builddir)/src/libpassmess_dbg.la
  parallel_unit_dbg_LDFLAGS = $(top_builddir)/src/libpassmess_dbg.la

  TESTS += message_tag_unit-dbg
  TESTS += packed_range_unit-dbg
  TESTS += parallel_sync_unit-dbg
  TESTS += parallel_unit-dbg
endif

if BUILD_DEVEL_MODE
  message_tag_unit_devel_SOURCES = message_tag_unit.C
  packed_range_unit_devel_SOURCES = packed_range_unit.C
  parallel_sync_unit_devel_SOURCES = parallel_sync_unit.C
  parallel_unit_devel_SOURCES = parallel_unit.C
  message_tag_unit_devel_LDFLAGS = $(top_builddir)/src/libpassmess_devel.la
  packed_range_unit_devel_LDFLAGS = $(top_builddir)/src/libpassmess_devel.la
  parallel_sync_unit_devel_LDFLAGS = $(top_builddir)/src/libpassmess_devel.la
  parallel_unit_devel_LDFLAGS = $(top_builddir)/src/libpassmess_devel.la

  TESTS += message_tag_unit-devel
  TESTS += packed_range_unit-devel
  TESTS += parallel_sync_unit-devel
  TESTS += parallel_unit-devel
endif

if BUILD_OPT_MODE
  message_tag_unit_opt_SOURCES = message_tag_unit.C
  packed_range_unit_opt_SOURCES = packed_range_unit.C
  parallel_sync_unit_opt_SOURCES = parallel_sync_unit.C
  parallel_unit_opt_SOURCES = parallel_unit.C
  message_tag_unit_opt_LDFLAGS = $(top_builddir)/src/libpassmess_opt.la
  packed_range_unit_opt_LDFLAGS = $(top_builddir)/src/libpassmess_opt.la
  parallel_sync_unit_opt_LDFLAGS = $(top_builddir)/src/libpassmess_opt.la
  parallel_unit_opt_LDFLAGS = $(top_builddir)/src/libpassmess_opt.la

  TESTS += message_tag_unit-opt
  TESTS += packed_range_unit-opt
  TESTS += parallel_sync_unit-opt
  TESTS += parallel_unit-opt
endif

check_PROGRAMS = $(TESTS)

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
