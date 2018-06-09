pkg_cppflags = -I$(srcdir)/fpoptimizer


FEATURE_FLAGS =
FEATURE_FLAGS += -DFP_ENABLE_EVAL
#FEATURE_FLAGS += -DFP_SUPPORT_TR1_MATH_FUNCS

#FEATURE_FLAGS += -DFP_USE_THREAD_SAFE_EVAL
#FEATURE_FLAGS += -DFP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA
#FEATURE_FLAGS += -DFP_NO_EVALUATION_CHECKS
#FEATURE_FLAGS += -D_GLIBCXX_DEBUG
#FEATURE_FLAGS += -DFP_DISABLE_SHORTCUT_LOGICAL_EVALUATION
FEATURE_FLAGS += -DFP_SUPPORT_FLOAT_TYPE
FEATURE_FLAGS += -DFP_SUPPORT_LONG_DOUBLE_TYPE
FEATURE_FLAGS += -DFP_SUPPORT_LONG_INT_TYPE
#FEATURE_FLAGS += -DFP_SUPPORT_MPFR_FLOAT_TYPE
#FEATURE_FLAGS += -DFP_SUPPORT_GMP_INT_TYPE
FEATURE_FLAGS += -DFP_SUPPORT_COMPLEX_DOUBLE_TYPE
FEATURE_FLAGS += -DFP_SUPPORT_COMPLEX_FLOAT_TYPE
FEATURE_FLAGS += -DFP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
#FEATURE_FLAGS += -DFP_USE_STRTOLD

if FPARSER_SUPPORT_DEBUGGING
  FEATURE_FLAGS += -DFUNCTIONPARSER_SUPPORT_DEBUGGING
endif

# We now explicitly require compiler support for std::move
FEATURE_FLAGS += -DFP_SUPPORT_CXX11_MOVE

EXTRA_DIST      = # none, append below
BUILT_SOURCES   = # none, append below
CLEANFILES      = # none, append below
DISTCLEANFILES  = # none, append below
noinst_PROGRAMS = # none, append below

# generally we do not install contributed headers, but for fparser.hh we
# must because it can optionally be included directly by our
# parsed_function.h
includedir = $(prefix)/include/libmesh

include_HEADERS = fparser.hh fparser_ad.hh


EXTRA_DIST += util/bytecoderules.dat util/bytecoderules_header.txt


# fpoptimizer/grammar_data.cc is a generated source, and this rule specifies how we build it.
fpoptimizer/grammar_data.cc: \
		util/tree_grammar_parser \
		fpoptimizer/treerules.dat
	@echo "Parsing Grammar from "$(srcdir)/fpoptimizer/treerules.dat" using "$<" ..."
	@util/tree_grammar_parser < $(srcdir)/fpoptimizer/treerules.dat > $@

EXTRA_DIST += fpoptimizer/treerules.dat


FPOPTIMIZER_CC_FILES=\
	    lib/crc32.hh \
	    lib/autoptr.hh \
	    lib/functional.hh \
	    fpoptimizer/hash.hh \
	    fpoptimizer/codetree.hh \
	    fpoptimizer/grammar.hh \
	    fpoptimizer/consts.hh \
	    fpoptimizer/optimize.hh \
	    fpoptimizer/opcodename.hh \
	    fpoptimizer/opcodename.cc \
	    fpoptimizer/bytecodesynth.hh \
	    fpoptimizer/bytecodesynth.cc \
	    fpoptimizer/valuerange.hh \
	    fpoptimizer/rangeestimation.hh \
	    fpoptimizer/constantfolding.hh \
	    fpoptimizer/logic_boolgroups.hh \
	    fpoptimizer/logic_collections.hh \
	    fpoptimizer/logic_ifoperations.hh \
	    fpoptimizer/logic_powoperations.hh \
	    fpoptimizer/logic_comparisons.hh \
	    fpoptimizer/codetree.cc \
	    fpoptimizer/debug.cc \
	    fpoptimizer/grammar.cc \
	    fpoptimizer/grammar_data.cc \
	    fpoptimizer/optimize.cc \
	    fpoptimizer/optimize_match.cc \
	    fpoptimizer/optimize_synth.cc \
	    fpoptimizer/optimize_debug.cc \
	    fpoptimizer/hash.cc \
	    fpoptimizer/makebytecode.cc \
	    fpoptimizer/readbytecode.cc \
	    fpoptimizer/constantfolding.cc \
	    fpoptimizer/valuerange.cc \
	    fpoptimizer/rangeestimation.cc \
	    fpoptimizer/transformations.cc \
	    fpoptimizer/cse.cc \
	    fpoptimizer/optimize_main.cc

EXTRA_DIST += $(FPOPTIMIZER_CC_FILES) fpoptimizer/fpoptimizer_header.txt fpoptimizer/fpoptimizer_footer.txt


pkg_sources = fparser.cc fparser_ad.cc Faddeeva.cc Faddeeva.hh
pkg_sources += lib/sha1.cpp
if FPARSER_SUPPORT_JIT
  pkg_cppflags += -DFPARSER_JIT_COMPILER="\"$(CXX) $(CXXFLAGS)\""
endif

pkg_cppflags += $(FEATURE_FLAGS)


# util/tree_grammar_parser is a utility which is required to build
# fpoptimizer/grammar_data.cc.  But this file itself is only needed
# when we are compiling with optimization support enabled.
#
# Note that we always need this rule, so 'make dist' will work
util_tree_grammar_parser_CXXFLAGS  = # nothing fancy - don't use our compile flags for this utility code lest it be horribly slow
util_tree_grammar_parser_SOURCES   = \
                                 util/tree_grammar_parser.yy fpoptimizer/opcodename.cc fpoptimizer/opcodename.hh \
                                 mpfr/MpfrFloat.hh \
                                 mpfr/GmpInt.hh \
                                 fparser_mpfr.hh \
                                 fparser.hh \
                                 fparser.cc \
                                 fpoptimizer/logic_comparisons.hh \
                                 fpoptimizer/constantfolding.hh \
                                 fpoptimizer/logic_ifoperations.hh \
                                 fpoptimizer/instantiate.hh \
                                 fpoptimizer/instantiate_for_ad.hh \
                                 fpoptimizer/opcodename.hh \
                                 fpoptimizer/consts.hh \
                                 fpoptimizer/optimize.hh \
                                 fpoptimizer/bytecodesynth.hh \
                                 fpoptimizer/valuerange.hh \
                                 fpoptimizer/codetree.hh \
                                 fpoptimizer/logic_collections.hh \
                                 fpoptimizer/hash.hh \
                                 fpoptimizer/grammar.hh \
                                 fpoptimizer/rangeestimation.hh \
                                 fpoptimizer/logic_boolgroups.hh \
                                 fpoptimizer/logic_powoperations.hh \
                                 fpconfig.hh \
                                 lib/functional.hh \
                                 lib/crc32.hh \
                                 lib/autoptr.hh \
                                 lib/sha1.h \
                                 lib/sha1.cpp \
                                 extrasrc/fpaux.hh \
                                 extrasrc/fptypes.hh \
                                 extrasrc/fp_opcode_add.inc \
                                 fparser_gmpint.hh extrasrc/fp_identifier_parser.inc

                                 # util/cpp_compress.hh

  util_bytecoderules_parser_CXXFLAGS = # nothing fancy - don't use our compile flags for this utility code lest it be horribly slow
  util_bytecoderules_parser_SOURCES  = util/bytecoderules_parser.cc

  # 'util/make_function_name_parser' is used to create extrasrc/fp_identifier_parser.inc,
  # but it seems to only define a snippet.  So for now we use the distributed version
  # of that file rather than treating it as a built source.
  # util_make_function_name_parser_SOURCES  = util/make_function_name_parser.cc util/cpp_compress.cc

#  util_cpp_compress_CXXFLAGS = # nothing fancy - don't use our compile flags for this utility code lest it be horribly slow
#  util_cpp_compress_SOURCES  = util/cpp_compress.hh util/cpp_compress.cc util/cpp_compress_main.cc


# when doing 'make clean' we need to remove the generated sources
CLEANFILES += fpoptimizer/grammar_data.cc
# in case they weren't conditionally cleaned:
DISTCLEANFILES += util/bytecoderules_parser
# DISTCLEANFILES += util/cpp_compress
DISTCLEANFILES += util/tree_grammar_parser


# extrasrc/fp_opcode_add.inc is a generated source, and this rule specifies how we build it.
extrasrc/fp_opcode_add.inc: \
		util/bytecoderules_parser \
		util/bytecoderules.dat \
		util/bytecoderules_header.txt
	@echo "Creating system-specific Bytecode in "$@" ..."
	@$(MKDIR_P) extrasrc
	@cat $(srcdir)/util/bytecoderules_header.txt > $@
	@util/bytecoderules_parser \
		< $(srcdir)/util/bytecoderules.dat \
		>> $@

BUILT_SOURCES += extrasrc/fp_opcode_add.inc

fparser.cc: extrasrc/fp_opcode_add.inc

if FPARSER_DEVEL

  BUILT_SOURCES += fpoptimizer/grammar_data.cc

  # Build up utility programs used to define parsing, not needed
  # for installation
  noinst_PROGRAMS           += util/tree_grammar_parser
  noinst_PROGRAMS           += util/bytecoderules_parser
  # noinst_PROGRAMS           += util/cpp_compress
  # noinst_PROGRAMS           += util/make_function_name_parser

  pkg_sources += $(FPOPTIMIZER_CC_FILES)

  # when doing 'make clean' we need to remove the generated sources
  CLEANFILES += $(BUILT_SOURCES)

else !FPARSER_DEVEL
  SUBDIRS = extrasrc

  pkg_sources += fpoptimizer.cc
endif

include $(top_srcdir)/contrib/Make.common
