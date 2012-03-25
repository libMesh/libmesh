include $(top_builddir)/Make.common

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

AUTOMAKE_OPTIONS = subdir-objects 

INCLUDES = -I$(srcdir)/include  -I$(srcdir)/fpoptimizer
EXTRA_DIST = # none, append below

# Build up utility programs used to define parsing, not needed
# for installation
noinst_PROGRAMS                   = util/bytecoderules_parser
util_bytecoderules_parser_SOURCES = util/bytecoderules_parser.cc

noinst_PROGRAMS           += util/cpp_compress
util_cpp_compress_SOURCES  = util/cpp_compress.hh util/cpp_compress.cc util/cpp_compress_main.cc

noinst_PROGRAMS +=  util/tree_grammar_parser
util_tree_grammar_parser_SOURCES = \
                                 util/tree_grammar_parser.yy fpoptimizer/opcodename.cc fpoptimizer/opcodename.hh \
                                 mpfr/MpfrFloat.hh \
                                 mpfr/GmpInt.hh \
                                 fparser_mpfr.hh \
                                 fparser.hh \
                                 fpoptimizer/logic_comparisons.hh \
                                 fpoptimizer/constantfolding.hh \
                                 fpoptimizer/logic_ifoperations.hh \
                                 fpoptimizer/instantiate.hh \
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
                                 util/cpp_compress.hh \
                                 fpconfig.hh \
                                 lib/functional.hh \
                                 lib/crc32.hh \
                                 lib/autoptr.hh \
                                 extrasrc/fpaux.hh \
                                 extrasrc/fptypes.hh \
                                 fparser_gmpint.hh

BUILT_SOURCES = extrasrc/fp_opcode_add.inc fpoptimizer/grammar_data.cc
CLEANFILES = $(BUILT_SOURCES)

extrasrc/fp_opcode_add.inc: \
		util/bytecoderules_parser \
		util/bytecoderules.dat \
		util/bytecoderules_header.txt \
		util/cpp_compress
	@echo "Creating Bytecode from "$<"..."
	@mkdir -p extrasrc
	@cat $(srcdir)/util/bytecoderules_header.txt > $@
	@util/bytecoderules_parser \
		< $(srcdir)/util/bytecoderules.dat \
		| util/cpp_compress \
		>> $@

EXTRA_DIST += util/bytecoderules.dat util/bytecoderules_header.txt

fpoptimizer/grammar_data.cc: \
		util/tree_grammar_parser \
		fpoptimizer/treerules.dat
	@echo "Parsing Grammar "$<"..."
	@util/tree_grammar_parser < $(srcdir)/fpoptimizer/treerules.dat > $@

EXTRA_DIST += fpoptimizer/treerules.dat

FEATURE_FLAGS =
FEATURE_FLAGS += -DFP_ENABLE_EVAL
#FEATURE_FLAGS += -DFP_SUPPORT_TR1_MATH_FUNCS

# # Mac OSX has trouble linking the optimizer stuff, so we disable
# # it on those systems.
# ifeq ($(findstring darwin,$(hostos)),darwin)
#   FEATURE_FLAGS += -DFP_NO_SUPPORT_OPTIMIZER
# endif

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
FEATURE_FLAGS += -DFP_USE_STRTOLD

FEATURE_FLAGS += -DFUNCTIONPARSER_SUPPORT_DEBUGGING

# generally we do not install contributed headers, but for fparser.hh we
# must because it can optionally be included directly by our
# parsed_function.h
includedir = $(prefix)/include/contrib

include_HEADERS       =  fparser.hh

noinst_LTLIBRARIES    = libfparser.la
libfparser_la_CPPFLAGS = $(AM_CPPFLAGS) $(FEATURE_FLAGS)
libfparser_la_CXXFLAGS = $(AM_CXXFLAGS)
libfparser_la_SOURCES  = \
                $(builddir)/extrasrc/fp_opcode_add.inc \
                fparser.cc \
		fpoptimizer/grammar_data.cc \
		fpoptimizer/optimize_main.cc \
		fpoptimizer/readbytecode.cc \
		fpoptimizer/makebytecode.cc \
		fpoptimizer/codetree.cc \
		fpoptimizer/grammar.cc \
		fpoptimizer/optimize.cc \
		fpoptimizer/optimize_match.cc \
		fpoptimizer/optimize_synth.cc \
		fpoptimizer/optimize_debug.cc \
		fpoptimizer/constantfolding.cc \
		fpoptimizer/valuerange.cc \
		fpoptimizer/rangeestimation.cc \
		fpoptimizer/opcodename.cc \
		fpoptimizer/bytecodesynth.cc \
		fpoptimizer/transformations.cc \
		fpoptimizer/cse.cc \
		fpoptimizer/debug.cc \
		fpoptimizer/hash.cc \
		fpoptimizer/consts.hh
