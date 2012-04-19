include $(top_builddir)/Make.common.@METHOD@

FEATURE_FLAGS =
FEATURE_FLAGS += -DFP_ENABLE_EVAL
#FEATURE_FLAGS += -DFP_SUPPORT_TR1_MATH_FUNCS

# some systems have trouble with the fparser optimizer
if FPARSER_NO_SUPPORT_OPTIMIZER
  FEATURE_FLAGS += -DFP_NO_SUPPORT_OPTIMIZER
endif

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


AM_CPPFLAGS = $(libmesh_CPPFLAGS) $(FEATURE_FLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

AUTOMAKE_OPTIONS = subdir-objects 

INCLUDES = -I$(srcdir)/fpoptimizer

EXTRA_DIST    = # none, append below
BUILT_SOURCES = # none, append below



# Build up utility programs used to define parsing, not needed
# for installation
noinst_PROGRAMS                    = util/bytecoderules_parser
util_bytecoderules_parser_CXXFLAGS = # nothing fancy - don't use our compile flags for this utility code lest it be horribly slow
util_bytecoderules_parser_SOURCES  = util/bytecoderules_parser.cc

# 'util/make_function_name_parser' is used to create extrasrc/fp_identifier_parser.inc,
# but it seems to only define a snippet.  So for now we use the distributed version
# of that file rather than treating it as a built source.
# noinst_PROGRAMS                        += util/make_function_name_parser
# util_make_function_name_parser_SOURCES  = util/make_function_name_parser.cc util/cpp_compress.cc 

noinst_PROGRAMS           += util/cpp_compress
util_cpp_compress_CXXFLAGS = -O2 # nothing fancy - don't use our compile flags for this utility code lest it be horribly slow
util_cpp_compress_SOURCES  = util/cpp_compress.hh util/cpp_compress.cc util/cpp_compress_main.cc


# util/tree_grammar_parser is a utility which is required to build 
# fpoptimizer/grammar_data.cc.  But this file itself is only needed
# when we are compiling with optimization support enabled.
# 
# using rule-specific CXXFLAGS here will compile target-specific objects,
# e.g. fpoptimizer/util_tree_grammar_parser-opcodename.o, which gets around
# the issue of 'object created both with libtool and without' that would
# otherwise occur when we build objects directly for our temporary utility
# programs and when constructing the final library.
if FPARSER_SUPPORT_OPTIMIZER
  noinst_PROGRAMS                 +=  util/tree_grammar_parser
endif
util_tree_grammar_parser_CXXFLAGS  = # nothing fancy - don't use our compile flags for this utility code lest it be horribly slow
util_tree_grammar_parser_SOURCES   = \
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
                                 fparser_gmpint.hh extrasrc/fp_identifier_parser.inc

# Eventually I think we want to build this generated source as well,
# but as of v4.4.3 the 'util/make_function_name_parser' only creates
# a snippet.
# extrasrc/fp_identifier_parser.inc: \
# 		util/make_function_name_parser
# 	@echo "Creating "$@
# 	@$(MKDIR_P) extrasrc
# 	@util/make_function_name_parser > $@



# extrasrc/fp_opcode_add.inc is a generated source, and this rule specifies how we build it.
BUILT_SOURCES += extrasrc/fp_opcode_add.inc
extrasrc/fp_opcode_add.inc: \
		util/bytecoderules_parser \
		util/bytecoderules.dat \
		util/bytecoderules_header.txt \
		util/cpp_compress
	@echo "Creating system-specific Bytecode in "$@" ..."
	@$(MKDIR_P) extrasrc
	@cat $(srcdir)/util/bytecoderules_header.txt > $@
	@util/bytecoderules_parser \
		< $(srcdir)/util/bytecoderules.dat \
		| util/cpp_compress \
		>> $@

EXTRA_DIST += util/bytecoderules.dat util/bytecoderules_header.txt



# fpoptimizer/grammar_data.cc is a generated source, and this rule specifies how we build it.
if FPARSER_SUPPORT_OPTIMIZER
  BUILT_SOURCES += fpoptimizer/grammar_data.cc
endif
fpoptimizer/grammar_data.cc: \
		util/tree_grammar_parser \
		fpoptimizer/treerules.dat
	@echo "Parsing Grammar from "$(srcdir)/fpoptimizer/treerules.dat" using "$<" ..."
	@util/tree_grammar_parser < $(srcdir)/fpoptimizer/treerules.dat > $@

EXTRA_DIST += fpoptimizer/treerules.dat



# generally we do not install contributed headers, but for fparser.hh we
# must because it can optionally be included directly by our
# parsed_function.h
includedir = $(prefix)/include/contrib

include_HEADERS       =  fparser.hh

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

# # fpoptimizer.cc is a generated source, and this rule specifies how we build it.
# # this rule is expanded from the upstream makefile because we need to account for
# # VPATH builds, in which the FPOPTIMIZER_CC_FILES are scattered around the source
# # and build trees.  So I have introduced the if/else and file_to_include construct.
# if FPARSER_SUPPORT_OPTIMIZER
#   BUILT_SOURCES += fpoptimizer.cc
# endif
# fpoptimizer.cc: fpoptimizer/fpoptimizer_header.txt \
# 		fpoptimizer/fpoptimizer_footer.txt \
# 	        extrasrc/fp_opcode_add.inc \
# 		$(FPOPTIMIZER_CC_FILES) \
# 		util/cpp_compress
# 	@echo "Generating source file "$@"..."
# 	@rm -f fpoptimizer.cc
# 	@cat $(srcdir)/fpoptimizer/fpoptimizer_header.txt  > $@
# 	@for file in $(FPOPTIMIZER_CC_FILES); do \
#                 file_to_include="" ;\
# 	        if (test -f $(srcdir)/$$file); then \
#                   file_to_include=$(srcdir)/$$file ; \
#                 elif (test -f $$file); then \
#                   file_to_include=$$file ; \
#                 else \
#                   echo "ERROR: Cannot find $$file in source or build tree!" ; \
# 	          exit 1 ; \
#                 fi ; \
# 		echo "#line 1 \"$$file_to_include\""; \
# 		$(SED) 's@^(#include \".*)@// line removed for fpoptimizer.cc: \\1@' < "$$file_to_include"; \
# 		echo; \
# 	done | $(SED) 's@BEGIN_EXPLICIT_INSTANTATION.*@@;s@.*END_EXPLICIT_INSTANTATION@@' \
# 	     | util/cpp_compress "lnxyceti" >> $@
# 	@cat $(srcdir)/fpoptimizer/fpoptimizer_footer.txt >> $@


# On OSX, compiling the $(FPOPTIMIZER_CC_FILES) files separately and linking
# them into libfparser.la seems to get around the multiple definition problem.
# note the upstream approach is to pack all the sources into fpoptimizer.cc 
# and compress it.
noinst_LTLIBRARIES    = libfparser.la
libfparser_la_SOURCES = fparser.cc
if FPARSER_SUPPORT_OPTIMIZER 
  libfparser_la_SOURCES += $(FPOPTIMIZER_CC_FILES) #fpoptimizer.cc
endif



# when doing 'make clean' we need to remove the generated sources
CLEANFILES = extrasrc/fp_opcode_add.inc fpoptimizer/grammar_data.cc fpoptimizer.cc
DISTCLEANFILES = util/tree_grammar_parser # in case it wasn't conditionally cleaned