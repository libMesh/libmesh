
include ../../Make.common

# Package name and additional flags required for library compilation
pkgname  := fparser

.PHONY: clean clobber distclean

target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext)

ifeq ($(enable-shared),yes)
  target := $(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext)
endif

ifeq ($(enable-fparser),no)
  all::
        @echo "<<< $(pkgname) is disabled, not building >>>"
else
  all:: $(target)
#  all:: $(target) testbed speedtest functioninfo
endif



# Modified original Makefile:

#===========================================================================
# This Makefile uses quite heavily GNU Make extensions, so it's probably
# hopeless to try to use it with other Make programs which do not have the
# same extensions.
#
# Also requires: rm, grep, sed and g++ (regardless of what CXX and LD are).
# The optimizer code generator requires bison.
#===========================================================================

RELEASE_VERSION=4.4.3

# The FP_FEATURE_FLAGS is set by run_full_release_testing.sh, but can be
# used otherwise as well.
ifeq ($(FP_FEATURE_FLAGS),)
FEATURE_FLAGS =
FEATURE_FLAGS += -DFP_ENABLE_EVAL
#FEATURE_FLAGS += -DFP_SUPPORT_TR1_MATH_FUNCS
#FEATURE_FLAGS += -DFP_NO_SUPPORT_OPTIMIZER
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
else
FEATURE_FLAGS = $(FP_FEATURE_FLAGS)
endif

FEATURE_FLAGS += -DFUNCTIONPARSER_SUPPORT_DEBUGGING

fp_CPPFLAGS=$(FEATURE_FLAGS) $(libmesh_CPPFLAGS)
fp_CXXFLAGS=$(libmesh_CXXFLAGS)
fp_LDFLAGS=$(libmesh_LDFLAGS)

ifneq (,$(findstring -DFP_SUPPORT_MPFR_FLOAT_TYPE,$(FEATURE_FLAGS)))
fp_LDFLAGS += -lgmp -lmpfr
ADDITIONAL_SRCS = mpfr/MpfrFloat.cc
ifneq (,$(findstring -DFP_SUPPORT_GMP_INT_TYPE,$(FEATURE_FLAGS)))
ADDITIONAL_SRCS += mpfr/GmpInt.cc
endif
else
ifneq (,$(findstring -DFP_SUPPORT_GMP_INT_TYPE,$(FEATURE_FLAGS)))
fp_LDFLAGS += -lgmp
ADDITIONAL_SRCS = mpfr/GmpInt.cc
endif
endif

ifneq (,$(findstring -DFP_USE_THREAD_SAFE_EVAL,$(FEATURE_FLAGS)))
BOOST_THREAD_LIB = -lboost_thread-mt
else
ifneq (,$(findstring -DFP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA,$(FEATURE_FLAGS)))
BOOST_THREAD_LIB = -lboost_thread-mt
endif
endif

fp_CPPFLAGS += -I$(LIBMESH_DIR)/contrib/fparser

SRC = 	fparser.cc \
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
		$(ADDITIONAL_SRCS)

OBJS   	:= $(patsubst %.cc, %.$(obj-suffix), $(SRC))

testbed: testbed.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) $(BOOST_THREAD_LIB) -lfparser

fpoptimizer.$(obj-suffix): fpoptimizer.cc

testbed_release: testbed.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) $(BOOST_THREAD_LIB) -lfparser

speedtest: util/speedtest.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

speedtest_release: util/speedtest.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

examples/example: examples/example.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

examples/example2: examples/example2.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

ftest: util/ftest.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

powi_speedtest: util/powi_speedtest.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

koe: koe.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

functioninfo: util/functioninfo.$(obj-suffix) $(target)
	$(CXX) -o $@ $< $(fp_LDFLAGS) -lfparser

fpoptimizer/grammar_data.cc: \
		util/tree_grammar_parser \
		fpoptimizer/treerules.dat
	util/tree_grammar_parser < fpoptimizer/treerules.dat > $@

extrasrc/fp_opcode_add.inc: \
		util/bytecoderules_parser \
		util/bytecoderules.dat \
		util/bytecoderules_header.txt \
		util/cpp_compress
	cat util/bytecoderules_header.txt > $@
	util/bytecoderules_parser \
		< util/bytecoderules.dat \
		| util/cpp_compress \
		>> $@

tests/make_tests: tests/make_tests.$(obj-suffix) \
                  util/cpp_compress.$(obj-suffix)
	$(CXX) -o $@ $^ $(fp_LDFLAGS)

testbed_tests.inc: tests/make_tests
	tests/make_tests tests/*/* -o $@

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

fpoptimizer.cc: fpoptimizer/fpoptimizer_header.txt \
		fpoptimizer/fpoptimizer_footer.txt \
		$(FPOPTIMIZER_CC_FILES) \
		util/cpp_compress
	rm -f fpoptimizer.cc
	cat fpoptimizer/fpoptimizer_header.txt  > $@
	for file in $(FPOPTIMIZER_CC_FILES); do \
		echo "#line 1 \"$$file\""; \
		sed -r "s@^(#include \".*)@// line removed for fpoptimizer.cc: \\1@" < "$$file"; \
		echo; \
	done | sed 's@BEGIN_EXPLICIT_INSTANTATION.*@@;s@.*END_EXPLICIT_INSTANTATION@@' \
	     | util/cpp_compress "lnxyceti" >> $@
	#     >> $@
	cat fpoptimizer/fpoptimizer_footer.txt >> $@

util/tree_grammar_parser: \
		util/tree_grammar_parser.$(obj-suffix) \
		fpoptimizer/opcodename.$(obj-suffix)
	$(CXX) -o $@ $^ $(fp_LDFLAGS)

util/tree_grammar_parser.cc: \
		util/tree_grammar_parser.y
	bison --output=$@ $<
	sed -i 's/ *$$//' $@

util/cpp_compress: \
		util/cpp_compress.$(obj-suffix) util/cpp_compress_main.$(obj-suffix)
	$(CXX) -o $@ $^ $(fp_LDFLAGS)

util/bytecoderules_parser: util/bytecoderules_parser.$(obj-suffix)
	$(CXX) -o $@ $^ $(fp_LDFLAGS)


util/version_changer: util/version_changer.cc
	g++ -O3 $^ -s -o $@ $(fp_LDFLAGS) $(fp_CXXFLAGS) $(fp_CPPFLAGS)

util/make_function_name_parser: util/make_function_name_parser.cc util/cpp_compress.$(obj-suffix)
	g++ -O3 $^ -s -o $@ $(fp_LDFLAGS) $(fp_CXXFLAGS) $(fp_CPPFLAGS)

util/powi_opt: \
		util/powi_opt.$(obj-suffix) \
		fpoptimizer/hash.$(obj-suffix) \
		fpoptimizer/constantfolding.$(obj-suffix) \
		fpoptimizer/codetree.$(obj-suffix) \
		fpoptimizer/valuerange.$(obj-suffix) \
		fpoptimizer/rangeestimation.$(obj-suffix)
	g++ -O3 $^ -s -o $@ $(fp_LDFLAGS) $(fp_CXXFLAGS) $(fp_CPPFLAGS)

util/create_testrules_for_optimization_rules: \
		util/create_testrules_for_optimization_rules.cc \
		fpoptimizer/grammar_data.$(obj-suffix) \
		fpoptimizer/opcodename.$(obj-suffix) \
		fpoptimizer/grammar.$(obj-suffix)
	g++ -O3 $^ -s -o $@ $(fp_LDFLAGS) $(fp_CXXFLAGS) $(fp_CPPFLAGS)

fpoptimizer_tests.sh: util/create_testrules_for_optimization_rules
	./$< > $@
	chmod +x $@

set_version_string: util/version_changer
	util/version_changer $(RELEASE_VERSION) fparser.cc \
		fparser.hh fparser_mpfr.hh fparser_gmpint.hh fpconfig.hh \
		fpoptimizer.cc extrasrc/fptypes.hh extrasrc/fpaux.hh \
		extrasrc/fp_opcode_add.inc \
		fpoptimizer/fpoptimizer_header.txt \
		util/bytecoderules_header.txt \
		docs/fparser.html webpage/index.html

pack: set_version_string distro_pack devel_pack

distro_pack: $(RELEASE_PACK_FILES)
	zip -9 fparser$(RELEASE_VERSION).zip $(RELEASE_PACK_FILES)
	# Use KZIP&ZIPMIX (advsys.net/ken), if possible, to create a smaller zip file
	if which kzip; then \
	  rm -rf fparser-$(RELEASE_VERSION);\
	  mkdir fparser-$(RELEASE_VERSION); \
	  tar cf - $(RELEASE_PACK_FILES) | tar -x -v -C fparser-$(RELEASE_VERSION) -f -; \
	  for s in -b0 -b128 -b256 -b512 -b1024 \
	           -rn -rn -rn -rn -rn -rn -rn -rn \
	           -rn -rn -rn -rn -rn -rn -rn -rn; do \
	    (cd fparser-$(RELEASE_VERSION); \
	    kzip -r -y "$$s" ../fparser$(RELEASE_VERSION)-tmp.zip * );\
	    DeflOpt ../fparser$(RELEASE_VERSION)-tmp.zip; \
	    zipmix -y fparser$(RELEASE_VERSION).zip \
	    	      fparser$(RELEASE_VERSION)-tmp.zip \
	    	      fparser$(RELEASE_VERSION)-tmp2.zip; \
	    if [ -f fparser$(RELEASE_VERSION)-tmp2.zip ]; then \
	      mv -f fparser$(RELEASE_VERSION)-tmp2.zip fparser$(RELEASE_VERSION).zip; \
	    fi; \
	    ls -al fparser$(RELEASE_VERSION)*.zip; \
	  done; \
	  rm -f fparser$(RELEASE_VERSION)-tmp.zip; \
	fi

devel_pack:
	tar --exclude='*~' \
		--transform="s|^|fparser_$(RELEASE_VERSION)_devel/|" \
		-cjvf fparser$(RELEASE_VERSION)_devel.tar.bz2 \
		Makefile examples/example.cc examples/example2.cc fparser.cc \
		fparser.hh fparser_mpfr.hh fparser_gmpint.hh \
		fpconfig.hh extrasrc/fptypes.hh extrasrc/fpaux.hh \
		extrasrc/fp_opcode_add.inc \
		extrasrc/fp_identifier_parser.inc \
		testbed_tests.inc \
		util/speedtest.cc testbed.cc \
		tests/*.cc tests/*.txt tests/*/* \
		util/*.cc util/*.hh util/*.dat util/*.txt util/*.y \
		docs/fparser.html docs/style.css docs/lgpl.txt docs/gpl.txt \
		fpoptimizer/*.hh fpoptimizer/*.cc \
		fpoptimizer/*.dat \
		fpoptimizer/*.txt \
		lib/*.hh \
		mpfr/MpfrFloat.hh mpfr/MpfrFloat.cc \
		mpfr/GmpInt.hh mpfr/GmpInt.cc \
		run_full_release_testing.sh \
		util/functioninfo.cc

clean:
	rm -f	testbed testbed_release \
		speedtest speedtest_release \
		functioninfo \
		examples/example examples/example2 ftest powi_speedtest \
		util/tree_grammar_parser \
		tests/make_tests \
		util/bytecoderules_parser \
		util/cpp_compress \
		util/make_function_name_parser \
		examples/*.$(obj-suffix) \
		fpoptimizer/*.$(obj-suffix) \
		tests/*.$(obj-suffix) \
		mpfr/*.$(obj-suffix) \
		util/*.$(obj-suffix) \
		*.$(obj-suffix) \
		.dep \
		util/tree_grammar_parser.output

release_clean:
	rm -f testbed_release speedtest_release \
		testbed.$(obj-suffix) fparser.$(obj-suffix) fpoptimizer.$(obj-suffix)

distclean: clean
	rm -f	*~
	@rm -f $(libpath)/*/lib$(pkgname).*

TESTBED_TEST_FILES = $(wildcard tests/*/*)
testbed_tests.inc: $(TESTBED_TEST_FILES)

.dep:
	echo -n '' > .dep
	- g++ -MM -MG $(fp_CPPFLAGS) $(wildcard *.cc) | sed "s|^\(.*\).o:|\1.$(obj-suffix):|" >> .dep
	- g++ -MM $(fp_CPPFLAGS) $(wildcard examples/*.cc) | sed "s|^\(.*\).o:|examples/\1.$(obj-suffix):|" >> .dep
	- g++ -MM $(fp_CPPFLAGS) $(wildcard fpoptimizer/*.cc) | sed "s|^\(.*\).o:|fpoptimizer/\1.$(obj-suffix):|" >> .dep
	- g++ -MM $(fp_CPPFLAGS) $(wildcard tests/*.cc) | sed "s|^\(.*\).o:|tests/\1.$(obj-suffix):|" >> .dep
	- g++ -MM $(fp_CPPFLAGS) $(wildcard util/*.cc) | sed "s|^\(.*\).o:|util/\1.$(obj-suffix):|" >> .dep
	- g++ -MM $(fp_CPPFLAGS) $(wildcard mpfr/*.cc) | sed "s|^\(.*\).o:|mpfr/\1.$(obj-suffix):|" >> .dep
#	- g++ -MM $(fp_CPPFLAGS) $(wildcard lib/*.cc) | sed "s|^\(.*\).o:|lib/\1.$(obj-suffix):|" >> .dep
	sed -i "s@`pwd`/@@" .dep

-include .dep

# More libMesh additions:

# Build static library
ifeq ($(findstring darwin,$(hostos)),darwin)
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@libtool -static -o $@ $(OBJS)
else
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(static_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(AR) rv $@ $(OBJS)
endif

# Build dynamic library
$(EXTERNAL_LIBDIR)/lib$(pkgname)$(shared_libext): $(OBJS)
	@echo "Linking "$@
	@$(shell mkdir -p $(EXTERNAL_LIBDIR))
	@$(libmesh_CXX) $(libmesh_CXXSHAREDFLAG) -o $@ $(OBJS) $(libmesh_LDFLAGS)

clobber: clean

%.$(obj-suffix) : %.cc
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(fp_CPPFLAGS) $(fp_CXXFLAGS) -c $< -o $@

