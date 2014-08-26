# if we are building the development package, there is
# nothing to do here, and fp_opcode_add.inc will get
# generated. however, if we are building the release
# version, we need to link it.

EXTRA_DIST = fp_opcode_add.inc.release

if FPARSER_RELEASE
    # handle linking source files
    BUILT_SOURCES  = fp_opcode_add.inc
    CLEANFILES = $(BUILT_SOURCES)

fp_opcode_add.inc: $(srcdir)/fp_opcode_add.inc.release
	$(AM_V_at)$(RM) $@
	$(AM_V_GEN)$(LN_S) -f $< $@

endif
