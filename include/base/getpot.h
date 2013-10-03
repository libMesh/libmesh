//  -*- c++ -*-
//  GetPot Version libMeshHPCT_fork-1.2                        Apr/14/2010
//  Based on "getpot-1.1.1.tgz" version from SourceForge
//
//  New code (C) 2009-2013 Roy Stogner, Karl Schulz
//
//  GetPot Version 1.0                                        Sept/13/2002
//
//  WEBSITE: http://getpot.sourceforge.net
//
//  This library is  free software; you can redistribute  it and/or modify
//  it  under  the terms  of  the GNU  Lesser  General  Public License  as
//  published by the  Free Software Foundation; either version  2.1 of the
//  License, or (at your option) any later version.
//
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//
//  You  should have  received a  copy of  the GNU  Lesser  General Public
//  License along  with this library; if  not, write to  the Free Software
//  Foundation, Inc.,  59 Temple Place,  Suite 330, Boston,  MA 02111-1307
//  USA
//
//  (C) 2001-2002 Frank R. Schaefer
//==========================================================================
#ifndef LIBMESH_GETPOT_H
#define LIBMESH_GETPOT_H

#if defined(WIN32) || defined(SOLARIS_RAW) || (__GNUC__ == 2) || defined(__HP_aCC)
#define strtok_r(a, b, c) strtok(a, b)
#endif // WINDOWS or SOLARIS or gcc 2.* or HP aCC

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream> // not every compiler distribution includes <iostream>
//                  // with <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include <cmath>
#include <cstdlib>
#include <ctime>

extern "C" {
#include <stdarg.h> // --> va_list and friends
#include <string.h> // --> strcmp, strncmp, strlen, strncpy
}

// Undefine USE_LIBMESH to avoid libMesh-specific code

#define USE_LIBMESH 1

#ifdef USE_LIBMESH

#include "libmesh/libmesh_common.h"

// We need a mutex to keep const operations thread-safe in the
// presence of mutable containers.  Right now GetPot supports a
// Threads::scoped_mutex wrapper around TBB, and we're assuming that
// users aren't doing any threaded GetPot usage when TBB threads are
// disabled.
#if !defined(GETPOT_DISABLE_MUTEX)
  #include "libmesh/threads.h"
  #define SCOPED_MUTEX  libMesh::Threads::spin_mutex::scoped_lock lock(_getpot_mtx)
#else
  #define SCOPED_MUTEX
#endif

#define getpot_cerr libMesh::err
#define getpot_error() libmesh_error()
#define getpot_file_error(filename) libmesh_file_error(filename)
#define getpot_cast_int libMesh::libmesh_cast_int

#else // USE_LIBMESH

// Currently threaded GetPot use is only supported via libMesh Threads
#define SCOPED_MUTEX

#define getpot_cerr std::cerr
#define getpot_error() throw std::runtime_error(std::string("GetPot Error"))
#define getpot_file_error(filename) getpot_error()
#define getpot_cast_int static_cast

#endif



typedef  std::vector<std::string>  STRING_VECTOR;

#define victorate(TYPE, VARIABLE, ITERATOR)                        \
  std::vector<TYPE>::const_iterator ITERATOR = (VARIABLE).begin(); \
  for(; (ITERATOR) != (VARIABLE).end(); (ITERATOR)++)

// We allow GETPOT_NAMESPACE to be defined before this file is
// included; if libraries using two different versions of GetPot might
// be linked together, the result may be unsafe unless they're put in
// different namespaces.
#ifdef GETPOT_NAMESPACE
namespace GETPOT_NAMESPACE {
#endif

class GetPot {
    //--------
    inline void _basic_initialization();
public:
    // (*) constructors, destructor, assignment operator -----------------------
    inline GetPot();
    inline GetPot(const GetPot&);
    inline GetPot(const int argc_, const char* const* argv_,
		  const char* FieldSeparator=0x0);
    inline GetPot(const char* FileName,
		  const char* CommentStart=0x0, const char* CommentEnd=0x0,
		  const char* FieldSeparator=0x0);
    inline GetPot(const std::string& FileName,
		  const std::string& CommentStart   = std::string("#"),
                  const std::string& CommentEnd     = std::string("\n"),
		  const std::string& FieldSeparator = std::string(" \t\n"));
    inline ~GetPot();
    inline GetPot& operator=(const GetPot&);

    // Re-initialization methods
    inline void parse_command_line(const int argc_, const char * const* argv_,
                                   const char* FieldSeparator =0x0);
    inline void parse_input_file(const std::string& FileName,
                                 const std::string& CommentStart=std::string("#"),
                                 const std::string& CommentEnd=std::string("\n"),
                                 const std::string& FieldSeparator=std::string(" \t\n"));

    // (*) absorbing contents of another GetPot object
    inline void            absorb(const GetPot& Other);
    //     -- for ufo detection: recording requested arguments, options etc.
    inline void            clear_requests();
    inline void            disable_request_recording() { request_recording_f = false; }
    inline void            enable_request_recording()  { request_recording_f = true; }

    // (*) direct access to command line arguments -----------------------------
    inline const char*     operator[](unsigned Idx) const;
    template <typename T>
    inline T               get(unsigned Idx, const T&    Default) const;
    inline const char*     get(unsigned Idx, const char* Default) const;
    inline unsigned        size() const;

    // (*) flags ---------------------------------------------------------------
    inline bool            options_contain(const char* FlagList) const;
    inline bool            argument_contains(unsigned Idx, const char* FlagList) const;

    // (*) variables -----------------------------------------------------------
    //     -- check for a variable
    inline bool            have_variable(const char* VarName) const;
    inline bool            have_variable(const std::string& VarName) const;
    //     -- scalar values
    template<typename T>
    inline T               operator()(const char* VarName, const T&    Default) const;
    template<typename T>
    inline T               operator()(const std::string& VarName, const T&    Default) const;
    inline const char*     operator()(const char* VarName, const char* Default) const;
    inline const char*     operator()(const std::string& VarName, const char* Default) const;
    //     -- vectors
    template<typename T>
    inline T               operator()(const char* VarName, const T&    Default, unsigned Idx) const;
    template<typename T>
    inline T               operator()(const std::string& VarName, const T&    Default, unsigned Idx) const;
    inline const char*     operator()(const char* VarName, const char* Default, unsigned Idx) const;
    inline const char*     operator()(const std::string& VarName, const char* Default, unsigned Idx) const;
    // (*) access varibles, but error out if not present -----------------------
    //     -- scalar values
    template<typename T>
    inline T               get_value_no_default(const char* VarName, const T& Default) const;
    template<typename T>
    inline T               get_value_no_default(const std::string& VarName, const T& Default) const;
    inline const char*     get_value_no_default(const char* VarName, const char* Default) const;
    inline const char*     get_value_no_default(const std::string& VarName, const char* Default) const;
    //     -- vectors
    template<typename T>
    inline T               get_value_no_default(const char* VarName, const T&    Default, unsigned Idx) const;
    template<typename T>
    inline T               get_value_no_default(const std::string& VarName, const T&    Default, unsigned Idx) const;
    inline const char*     get_value_no_default(const char* VarName, const char* Default, unsigned Idx) const;
    inline const char*     get_value_no_default(const std::string& VarName, const char* Default, unsigned Idx) const;

    //     -- setting variables
    //                  i) from outside of GetPot (considering prefix etc.)
    //                  ii) from inside, use '_set_variable()' below
    template<typename T>
    inline void            set(const char* VarName, const T& Value, const bool Requested = true);
    template<typename T>
    inline void            set(const std::string& VarName, const T& Value, const bool Requested = true);
    inline void            set(const char* VarName, const char* Value, const bool Requested = true);
    inline void            set(const std::string& VarName, const char* Value, const bool Requested = true);

    inline unsigned        vector_variable_size(const char* VarName) const;
    inline unsigned        vector_variable_size(const std::string& VarName) const;
    inline STRING_VECTOR   get_variable_names() const;
    inline STRING_VECTOR   get_section_names() const;
    inline
    std::set<std::string>  get_overridden_variables() const;

    // (*) cursor oriented functions -------------------------------------------
    inline void            set_prefix(const char* Prefix) { prefix = std::string(Prefix); }
    inline bool            search_failed() const { return search_failed_f; }

    //     -- enable/disable search for an option in loop
    inline void            disable_loop() { search_loop_f = false; }
    inline void            enable_loop()  { search_loop_f = true; }

    //     -- reset cursor to position '1'
    inline void            reset_cursor();
    inline void            init_multiple_occurrence();

    //     -- search for a certain option and set cursor to position
    inline bool            search(const char* option);
    inline bool            search(const std::string& option);
    inline bool            search(unsigned No, const char* P, ...);
    //     -- get argument at cursor++
    template<typename T>
    inline T               next(const T&    Default);
    inline const char*     next(const char* Default);
    //     -- search for option and get argument at cursor++
    template<typename T>
    inline T               follow(const T&    Default, const char* Option);
    inline const char*     follow(const char* Default, const char* Option);
    //     -- search for one of the given options and get argument that follows it
    template<typename T>
    inline T               follow(const T&    Default, unsigned No, const char* Option, ...);
    inline const char*     follow(const char* Default, unsigned No, const char* Option, ...);
    //     -- directly followed arguments
    template<typename T>
    inline T               direct_follow(const T&    Default, const char* Option);
    inline const char*     direct_follow(const char* Default, const char* Option);

    // (*) nominus arguments ---------------------------------------------------
    inline void            reset_nominus_cursor();
    inline STRING_VECTOR   nominus_vector() const;
    inline unsigned        nominus_size() const {
	return getpot_cast_int<unsigned>(idx_nominus.size());
    }
    inline const char*     next_nominus();

    // (*) unidentified flying objects -----------------------------------------
    inline STRING_VECTOR   unidentified_arguments(unsigned Number, const char* Known, ...) const;
    inline STRING_VECTOR   unidentified_arguments(const std::set<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_arguments(const std::vector<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_arguments() const;

    inline STRING_VECTOR   unidentified_options(unsigned Number, const char* Known, ...) const;
    inline STRING_VECTOR   unidentified_options(const std::set<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_options(const std::vector<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_options() const;

    inline std::string     unidentified_flags(const char* Known,
					     int ArgumentNumber /* =-1 */) const;

    inline STRING_VECTOR   unidentified_variables(unsigned Number, const char* Known, ...) const;
    inline STRING_VECTOR   unidentified_variables(const std::set<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_variables(const std::vector<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_variables() const;

    inline STRING_VECTOR   unidentified_sections(unsigned Number, const char* Known, ...) const;
    inline STRING_VECTOR   unidentified_sections(const std::set<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_sections(const std::vector<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_sections() const;

    inline STRING_VECTOR   unidentified_nominuses(unsigned Number, const char* Known, ...) const;
    inline STRING_VECTOR   unidentified_nominuses(const std::set<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_nominuses(const std::vector<std::string>& Knowns) const;
    inline STRING_VECTOR   unidentified_nominuses() const;

    // (*) output --------------------------------------------------------------
    // Print everything
    inline int print(std::ostream &out_stream = std::cout) const;
    // Print everything after skipping skip_count arguments, with a
    // custom prefix.  skip_count defaults to 1 to handle the common
    // "executable input_file" command line case.
    inline int print(const char *custom_prefix,
		     std::ostream &out_stream = std::cout,
		     unsigned int skip_count=1) const;

private:
    // (*) Type Declaration ----------------------------------------------------
    struct variable {
	//-----------
	// Variable to be specified on the command line or in input files.
	// (i.e. of the form var='12 312 341')

	// -- constructors, destructors, assignment operator
	variable();
	variable(const variable&);
	variable(const char* Name, const char* Value, const char* FieldSeparator);
	~variable();
	variable& operator=(const variable& Other);

	void      take(const char* Value, const char* FieldSeparator);

	// -- get a specific element in the string vector
	//    (return 0 if not present)
	const std::string*  get_element(unsigned Idx) const;

	// -- data memebers
	std::string       name;      // identifier of variable
	STRING_VECTOR     value;     // value of variable stored in vector
	std::string       original;  // value of variable as given on command line
    };

    // (*) member variables --------------------------------------------------------------
    std::string           prefix;          // prefix automatically added in queries
    std::string           section;         // (for dollar bracket parsing)
    STRING_VECTOR         section_list;    // list of all parsed sections
    //     -- argument vector
    STRING_VECTOR         argv;            // vector of command line arguments stored as strings
    unsigned              cursor;          // cursor for argv
    bool                  search_loop_f;   // shall search start at beginning after
    //                                     // reaching end of arg array ?
    bool                  search_failed_f; // flag indicating a failed search() operation
    //                                     // (e.g. next() functions react with 'missed')
    std::set<std::string> overridden_vars; // vector of variables that were supplied more than once during parsing

    //     --  nominus vector
    int                   nominus_cursor;  // cursor for nominus_pointers
    std::vector<unsigned> idx_nominus;     // indecies of 'no minus' arguments

    //     -- variables
    //       (arguments of the form "variable=value")
    std::vector<variable> variables;

    //     -- comment delimiters
    std::string           _comment_start;
    std::string           _comment_end;

    //     -- field separator (separating elements of a vector)
    std::string           _field_separator;

    //     -- helper functor for creating sets of C-style strings
    struct ltstr {
      bool operator()(const char* s1, const char* s2) const {
        return strcmp(s1, s2) < 0;
      }
    };

    //     -- we have some mutable non-thread-safe members, but we
    //        want to be able to call const member functions from
    //        multiple threads at once, so we'll wrap access to
    //        mutable objects in a mutex.
#if !defined(GETPOT_DISABLE_MUTEX)
    mutable libMesh::Threads::spin_mutex _getpot_mtx;
#endif

    //     -- some functions return a char pointer to a string created on the fly.
    //        this container makes them 'available' until the getpot object is destroyed.
    //        user codes are recommended to instead request std::string values.
    //        We use char* here because c_str() results are only
    //        guaranteed to remain valid until a non-const string
    //        method is called
    mutable std::set<const char*, ltstr> _internal_string_container;

    //     -- some functions return a char pointer to a temporarily existing string
    //        this function adds them to our container
    const char*    _internal_managed_copy(const std::string& Arg) const;

    //     -- keeping track about arguments that are requested, so that the UFO detection
    //        can be simplified
    mutable std::set<std::string> _requested_arguments;
    mutable std::set<std::string> _requested_variables;
    mutable std::set<std::string> _requested_sections;

    bool            request_recording_f;   // speed: request recording can be turned off

    //     -- if an argument is requested record it and the 'tag' the section branch to which
    //        it belongs. Caution: both functions mark the sections as 'tagged'.
    //        These are "const" functions but they do modify the
    //        mutable _requested_* members
    void                      _record_argument_request(const std::string& Arg) const;
    void                      _record_variable_request(const std::string& Arg) const;

    // (*) helper functions ----------------------------------------------------
    //                  set variable from inside GetPot (no prefix considered)
    inline void               _set_variable(const std::string& VarName,
					    const std::string& Value,
					    const bool Requested);

    //     -- produce three basic data vectors:
    //          - argument vector
    //          - nominus vector
    //          - variable dictionary
    inline void               _parse_argument_vector(const STRING_VECTOR& ARGV);

    //     -- helpers for argument list processing
    //        * search for a variable in 'variables' array
    inline const variable*    _find_variable(const char*) const;
    //        * search (and record request) for a variable in 'variables' array
    inline const variable*    _request_variable(const char*) const;
    //        * support finding directly followed arguments
    inline const char*        _match_starting_string(const char* StartString);
    //        * support search for flags in a specific argument
    inline bool               _check_flags(const std::string& Str, const char* FlagList) const;
    //        * type conversion if possible
    template<typename T>
    inline T                  _convert_to_type(const std::string& String, const T& Default) const;
    inline std::string        _convert_to_type(const std::string& String, const char* Default) const;
    template<typename T>
    inline T                  _convert_to_type_no_default(const char* VarName, const std::string& String, const T& Default) const;
    inline std::string        _convert_to_type_no_default(const char* VarName, const std::string& String, const char* Default) const;
    //        * prefix extraction
    const std::string         _get_remaining_string(const std::string& String,
						     const std::string& Start) const;
    //        * search for a specific string
    inline bool               _search_string_vector(const STRING_VECTOR& Vec,
						     const std::string& Str) const;

    //     -- helpers to parse input file
    //        create an argument vector based on data found in an input file, i.e.:
    //           1) delete comments (in between '_comment_start' '_comment_end')
    //           2) contract assignment expressions, such as
    //                   my-variable   =    '007 J. B.'
    //             into
    //                   my-variable='007 J. B.'
    //           3) interprete sections like '[../my-section]' etc.
    inline void               _skip_whitespace(std::istream& istr);
    inline const std::string  _get_next_token(std::istream& istr);
    inline const std::string  _get_string(std::istream& istr);
    inline const std::string  _get_until_closing_bracket(std::istream& istr);
    inline const std::string  _get_until_closing_square_bracket(std::istream& istr);

    inline STRING_VECTOR      _read_in_stream(std::istream& istr);
    inline STRING_VECTOR      _read_in_file(const std::string& FileName);
    inline std::string        _process_section_label(const std::string& Section,
						      STRING_VECTOR& section_stack);

    //      -- dollar bracket expressions
    std::string               _DBE_expand_string(const std::string& str);
    std::string               _DBE_expand(const std::string& str);
    const GetPot::variable*   _DBE_get_variable(const std::string& str);
    STRING_VECTOR             _DBE_get_expr_list(const std::string& str, const unsigned ExpectedNumber);

    template <typename T>
    static std::string _convert_from_type(const T& Value) {
      std::ostringstream out_string;
      out_string << Value;
      return out_string.str();
    }

    static STRING_VECTOR _get_section_tree(const std::string& FullPath) {
	// -- cuts a variable name into a tree of sub-sections. this is requested for recording
	//    requested sections when dealing with 'ufo' detection.
	STRING_VECTOR   result;
	for (std::size_t pos = 0; pos != FullPath.size(); ++pos) {
	    if( FullPath[pos] == '/' ) {
		result.push_back(FullPath.substr(0,pos));
	    }
	}

	return result;
    }
};


///////////////////////////////////////////////////////////////////////////////
// (*) constructors, destructor, assignment operator
//.............................................................................
//
inline void
GetPot::_basic_initialization()
{
    cursor = 0;              nominus_cursor = -1;
    search_failed_f = true;  search_loop_f = true;
    prefix = "";             section = "";

    // automatic request recording for later ufo detection
    request_recording_f = true;

    // comment start and end strings
    _comment_start = std::string("#");
    _comment_end   = std::string("\n");

    // default: separate vector elements by whitespaces
    _field_separator = " \t\n";
}

inline
GetPot::GetPot() :
  prefix(),
  section(),
  section_list(),
  argv(),
  cursor(),
  search_loop_f(),
  search_failed_f(),
  nominus_cursor(),
  idx_nominus(),
  variables(),
  _comment_start(),
  _comment_end(),
  _field_separator(),
#if !defined(GETPOT_DISABLE_MUTEX)
  _getpot_mtx(),
#endif
  _internal_string_container(),
  _requested_arguments(),
  _requested_variables(),
  _requested_sections(),
  request_recording_f()
{
    _basic_initialization();
}

inline
GetPot::GetPot(const int argc_, const char * const * argv_,
	       const char* FieldSeparator /* =0x0 */) :
    // leave 'char**' non-const to honor less capable compilers ...
  prefix(),
  section(),
  section_list(),
  argv(),
  cursor(),
  search_loop_f(),
  search_failed_f(),
  nominus_cursor(),
  idx_nominus(),
  variables(),
  _comment_start(),
  _comment_end(),
  _field_separator(),
#if !defined(GETPOT_DISABLE_MUTEX)
  _getpot_mtx(),
#endif
  _internal_string_container(),
  _requested_arguments(),
  _requested_variables(),
  _requested_sections(),
  request_recording_f()
{
    this->parse_command_line(argc_, argv_, FieldSeparator);
}


inline void
GetPot::parse_command_line(const int argc_, const char * const * argv_,
                           const char* FieldSeparator /* =0x0 */)
    // leave 'char**' non-const to honor less capable compilers ...
{
    _basic_initialization();

    // if specified -> overwrite default string
    if( FieldSeparator ) _field_separator = std::string(FieldSeparator);

    // -- make an internal copy of the argument list:
    STRING_VECTOR _apriori_argv;
    // -- for the sake of clarity: we do want to include the first
    //    argument of the first parsing source in the argument vector!
    //    it will not be a nominus argument, though. This gives us a
    //    minimum vector size of one which facilitates error checking
    //    in many functions. Also the user will be able to retrieve
    //    the name of his application or input file by "get[0]"
    _apriori_argv.push_back(std::string(argv_[0]));
    int i=1;
    for(; i<argc_; i++) {
	std::string tmp(argv_[i]);   // recall the problem with temporaries,
	_apriori_argv.push_back(tmp);       // reference counting in arguement lists ...
    }
    _parse_argument_vector(_apriori_argv);
}


inline
GetPot::GetPot(const char* FileName,
	       const char* CommentStart  /* = 0x0 */, const char* CommentEnd /* = 0x0 */,
	       const char* FieldSeparator/* = 0x0 */) :
  prefix(),
  section(),
  section_list(),
  argv(),
  cursor(),
  search_loop_f(),
  search_failed_f(),
  nominus_cursor(),
  idx_nominus(),
  variables(),
  _comment_start(),
  _comment_end(),
  _field_separator(),
#if !defined(GETPOT_DISABLE_MUTEX)
  _getpot_mtx(),
#endif
  _internal_string_container(),
  _requested_arguments(),
  _requested_variables(),
  _requested_sections(),
  request_recording_f()
{
  const std::string& StrCommentStart   = CommentStart   ? CommentStart   : std::string("#");
  const std::string& StrCommentEnd     = CommentEnd     ? CommentEnd     : std::string("\n");
  const std::string& StrFieldSeparator = FieldSeparator ? FieldSeparator : std::string(" \t\n");
  this->parse_input_file(FileName, StrCommentStart, StrCommentEnd, StrFieldSeparator);
}



inline
GetPot::GetPot(const std::string& FileName,
               const std::string& CommentStart,
               const std::string& CommentEnd,
               const std::string& FieldSeparator) :
  prefix(),
  section(),
  section_list(),
  argv(),
  cursor(),
  search_loop_f(),
  search_failed_f(),
  nominus_cursor(),
  idx_nominus(),
  variables(),
  _comment_start(),
  _comment_end(),
  _field_separator(),
#if !defined(GETPOT_DISABLE_MUTEX)
  _getpot_mtx(),
#endif
  _internal_string_container(),
  _requested_arguments(),
  _requested_variables(),
  _requested_sections(),
  request_recording_f()
{
    this->parse_input_file(FileName, CommentStart, CommentEnd, FieldSeparator);
}



inline void
GetPot::parse_input_file(const std::string& FileName,
                         const std::string& CommentStart,
                         const std::string& CommentEnd,
                         const std::string& FieldSeparator)
{
    _basic_initialization();

    // overwrite default strings
    _comment_start = std::string(CommentStart);
    _comment_end   = std::string(CommentEnd);
    _field_separator = FieldSeparator;

    STRING_VECTOR _apriori_argv;
    // -- the first element of the argument vector stores the name of
    //    the first parsing source; however, this element is not
    //    parsed for variable assignments or nominuses.
    //
    //    Regardless, we don't add more than one name to the argument
    //    vector.
    _apriori_argv.push_back(FileName);

    STRING_VECTOR args = _read_in_file(FileName.c_str());
    _apriori_argv.insert(_apriori_argv.begin()+1, args.begin(), args.end());
    _parse_argument_vector(_apriori_argv);
}

inline
GetPot::GetPot(const GetPot& Other) :
  prefix(Other.prefix),
  section(Other.section),
  section_list(Other.section_list),
  argv(Other.argv),
  cursor(Other.cursor),
  search_loop_f(Other.search_loop_f),
  search_failed_f(Other.search_failed_f),
  overridden_vars(),
  nominus_cursor(Other.nominus_cursor),
  idx_nominus(Other.idx_nominus),
  variables(Other.variables),
  _comment_start(Other._comment_start),
  _comment_end(Other._comment_end),
  _field_separator(Other._field_separator),
#if !defined(GETPOT_DISABLE_MUTEX)
  _getpot_mtx(Other._getpot_mtx),
#endif
  _internal_string_container(),
  _requested_arguments(Other._requested_arguments),
  _requested_variables(Other._requested_variables),
  _requested_sections(Other._requested_sections),
  request_recording_f(Other.request_recording_f)
{
    std::set<const char*,ltstr>::const_iterator it =
      Other._internal_string_container.begin();

    const std::set<const char*,ltstr>::const_iterator end =
      Other._internal_string_container.end();

    for (; it != end; ++it) {
        const char* otherstr = *it;
        char* newcopy = new char[strlen(otherstr)+1];
        strncpy(newcopy, otherstr, strlen(otherstr)+1);
        this->_internal_string_container.insert(newcopy);
    }
}

inline
GetPot::~GetPot()
{
    // // may be some return strings had to be created, delete now !
    std::set<const char*, ltstr>::const_iterator        it = _internal_string_container.begin();
    const std::set<const char*, ltstr>::const_iterator end = _internal_string_container.end();
    for(; it != end; ++it)
        delete [] *it;
}

inline GetPot&
GetPot::operator=(const GetPot& Other)
{
    if (&Other == this) return *this;

    prefix               = Other.prefix;
    section              = Other.section;
    section_list         = Other.section_list;
    argv                 = Other.argv;
    cursor               = Other.cursor;
    search_loop_f        = Other.search_loop_f;
    search_failed_f      = Other.search_failed_f;
    nominus_cursor       = Other.nominus_cursor;
    overridden_vars      = Other.overridden_vars;
    idx_nominus          = Other.idx_nominus;
    variables            = Other.variables;
    _comment_start       = Other._comment_start;
    _comment_end         = Other._comment_end;
    _field_separator     = Other._field_separator;
#if !defined(GETPOT_DISABLE_MUTEX)
    _getpot_mtx          = Other._getpot_mtx;
#endif
    _requested_arguments = Other._requested_arguments;
    _requested_variables = Other._requested_variables;
    _requested_sections  = Other._requested_sections;
    request_recording_f  = Other.request_recording_f;

    std::set<const char*, ltstr>::const_iterator        my_it =
      _internal_string_container.begin();
    const std::set<const char*, ltstr>::const_iterator my_end =
      _internal_string_container.end();

    for(; my_it != my_end; ++my_it)
        delete [] *my_it;

    _internal_string_container.clear();

    std::set<const char*,ltstr>::const_iterator it =
      Other._internal_string_container.begin();
    const std::set<const char*,ltstr>::const_iterator end =
      Other._internal_string_container.end();

    for (; it != end; ++it) {
        const char* otherstr = *it;
        char* newcopy = new char[strlen(otherstr)+1];
        strncpy(newcopy, otherstr, strlen(otherstr)+1);
        this->_internal_string_container.insert(newcopy);
    }

    return *this;
}


inline void
GetPot::absorb(const GetPot& Other)
{
    if (&Other == this) return;

    // variables that are not influenced by absorption:
    //               _comment_start
    //               _comment_end
    //               cursor
    //               nominus_cursor
    //               search_failed
    //               idx_nominus
    //               search_loop_f
    argv      = Other.argv;
    variables = Other.variables;

    if( request_recording_f ) {
        // Get a lock before touching anything mutable
        SCOPED_MUTEX;

	_requested_arguments.insert(Other._requested_arguments.begin(), Other._requested_arguments.end());
	_requested_variables.insert(Other._requested_variables.begin(), Other._requested_variables.end());
	_requested_sections.insert(Other._requested_sections.begin(), Other._requested_sections.end());
    }

}

inline void
GetPot::clear_requests()
{
    // Get a lock before touching anything mutable
    SCOPED_MUTEX;

    _requested_arguments.clear();
    _requested_variables.clear();
    _requested_sections.clear();
}

inline void
GetPot::_parse_argument_vector(const STRING_VECTOR& ARGV)
{
    if( ARGV.empty() ) return;

    // build internal databases:
    //   1) array with no-minus arguments (usually used as filenames)
    //   2) variable assignments:
    //             'variable name' '=' number | string
    STRING_VECTOR                 section_stack;
    STRING_VECTOR::const_iterator it = ARGV.begin();


    section = "";

    // -- do not parse the first argument, so that this parsing source
    // name is not interpreted a s a nominus or so.  If we already
    // have parsed arguments, don't bother adding another parsing
    // source name
    if (argv.empty())
      argv.push_back(*it);
    ++it;

    // -- loop over remaining arguments
    for(; it != ARGV.end(); ++it) {
	std::string arg = *it;

	if( arg.length() == 0 ) continue;

	// -- [section] labels and [include file] directives
	if( arg.length() > 1 && arg[0] == '[' && arg[arg.length()-1] == ']' ) {

            // Is this an include file directive?
            std::size_t include_pos = arg.find("include ", 1);
            if (include_pos != std::string::npos) {

	        const std::string includefile =
                  _DBE_expand_string(arg.substr(9, arg.length()-9-include_pos));

                this->parse_input_file
                  (includefile, _comment_start, _comment_end, _field_separator);
            } else {

	        // (*) sections are considered 'requested arguments'
	        if( request_recording_f ) {
                    // Get a lock before touching anything mutable
                    SCOPED_MUTEX;

                    _requested_arguments.insert(arg);
                }

	        const std::string Name = _DBE_expand_string(arg.substr(1, arg.length()-2));
	        section = _process_section_label(Name, section_stack);
	        // new section --> append to list of sections
	        if( find(section_list.begin(), section_list.end(), section) == section_list.end() )
		    if( section.length() != 0 ) section_list.push_back(section);
	        argv.push_back(arg);
            }
	}
	else {
	    arg = section + _DBE_expand_string(arg);
	    argv.push_back(arg);
	}

	// -- separate array for nominus arguments
	if( arg[0] != '-' )
	    idx_nominus.push_back(getpot_cast_int<unsigned>(argv.size()-1));

	// -- variables: does arg contain a '=' operator ?
	const std::size_t equals_pos = arg.find_first_of('=');
	if( equals_pos != std::string::npos ) {
	    // (*) record for later ufo detection
	    //     arguments carriying variables are always treated as 'requested' arguments.
	    //     unrequested variables have to be detected with the ufo-variable
	    //     detection routine.
	    if( request_recording_f ) {
		// Get a lock before touching anything mutable
		SCOPED_MUTEX;

		_requested_arguments.insert(arg);
            }

	    // => arg (from start to '=') = Name of variable
	    //        (from '=' to end)   = value of variable
	    _set_variable(arg.substr(0,equals_pos),
		          arg.substr(equals_pos+1), false);
	}
    }
}


inline STRING_VECTOR
GetPot::_read_in_file(const std::string& FileName)
{
    std::ifstream  i(FileName.c_str());

    // if( ! i ) return STRING_VECTOR();

    if (!i)
      libmesh_file_error(FileName);

    // argv[0] == the filename of the file that was read in
    return _read_in_stream(i);
}

inline STRING_VECTOR
GetPot::_read_in_stream(std::istream& istr)
{
    STRING_VECTOR  brute_tokens;
    while(istr) {
	_skip_whitespace(istr);
	const std::string Token = _get_next_token(istr);
        // Allow 'keyword =' to parse with an empty string as value.
        // Only break at EOF.
// 	if( Token.length() == 0 || Token[0] == EOF) break;
	if( Token[0] == EOF) break;
	brute_tokens.push_back(Token);
    }

    // -- reduce expressions of token1'='token2 to a single
    //    string 'token1=token2'
    // -- copy everything into 'argv'
    // -- arguments preceded by something like '[' name ']' (section)
    //    produce a second copy of each argument with a prefix '[name]argument'
    unsigned i1 = 0;
    unsigned i2 = 1;
    unsigned i3 = 2;

    STRING_VECTOR  arglist;
    while( i1 < brute_tokens.size() ) {
	// 1) concatenate 'abcdef' '=' 'efgasdef' to 'abcdef=efgasdef'
	// note: java.lang.String: substring(a,b) = from a to b-1
	//        C++ string:      substr(a,b)    = from a to a + b
        std::string result;
	if( i2 < brute_tokens.size() && brute_tokens[i2] == "=" ) {
	    if( i3 >= brute_tokens.size() )
              result = brute_tokens[i1] + brute_tokens[i2];
	    else
              result = brute_tokens[i1] + brute_tokens[i2] + brute_tokens[i3];
	    i1 = i3+1; i2 = i3+2; i3 = i3+3;
	}
        else if ( i2 < brute_tokens.size() &&
		  brute_tokens[i2].length() > 0 &&
		  brute_tokens[i2][0] == '=' ) {
          // This case should not be hit if '=' at the beginning of a word
          //   is always separated into its own word
          result = brute_tokens[i1] + brute_tokens[i2];
          i1 = i3; i2 = i3+1; i3 = i3+2;
        }
        else if ( i2 < brute_tokens.size() && brute_tokens[i1][brute_tokens[i1].size()-1] == '=' ) {
          result = brute_tokens[i1] + brute_tokens[i2];
          i1 = i3; i2 = i3+1; i3 = i3+2;
        }
	else {
          result = brute_tokens[i1];
          i1=i2; i2=i3; i3++;
	}
        // Now strip out any comment
        size_t comment_start_loc = result.find(_comment_start, 0);
        if (comment_start_loc != std::string::npos)
        {
          result = result.substr(0, comment_start_loc);
        }
        arglist.push_back(result);
    }
    return arglist;
}

inline void
GetPot::_skip_whitespace(std::istream& istr)
    // find next non-whitespace while deleting comments
{
    int tmp = istr.get();
    do {
	// -- search a non whitespace
	while( isspace(tmp) ) {
	    tmp = istr.get();
	    if( ! istr ) return;
	}

	// -- look if characters match the comment starter string
	const std::istream::pos_type  Pos = istr.tellg();
	unsigned    i=0;
	for(; i<_comment_start.length() ; i++) {
	    if( tmp != _comment_start[i] ) {
		istr.seekg(Pos);
		// -- one step more backwards, since 'tmp' already at non-whitespace
		istr.unget();
		return;
	    }

// RHS: Why is this here?  It breaks on empty comments
//	    tmp = istr.get();
//	    if( ! istr ) { istr.unget(); return; }
	}
	// 'tmp' contains last character of _comment_starter

	// -- comment starter found -> search for comment ender
	unsigned match_no=0;
	while(1+1 == 2) {
	    tmp = istr.get();
	    if( ! istr ) { istr.unget(); return; }

	    if( tmp == _comment_end[match_no] ) {
		match_no++;
		if( match_no == _comment_end.length() ) {
		    istr.unget();
		    break; // shuffle more whitespace, end of comment found
		}
	    }
	    else
		match_no = 0;
	}

	tmp = istr.get();

    } while( istr );
    istr.unget();
}

inline const std::string
GetPot::_get_next_token(std::istream& istr)
    // get next concatenates string token. consider quotes that embrace
    // whitespaces
{
    std::string token;
    int    tmp = 0;
    int    last_letter = 0;
    while(1+1 == 2) {
	last_letter = tmp; tmp = istr.get();
        if( tmp == '=' )
        {
          // Always break at '='.
          // This separates '=' at the beginning of a word into its own word.
          token += getpot_cast_int<char>(tmp);
          return token;
        }
        else if( tmp == EOF
	    || ((tmp == ' ' || tmp == '\t' || tmp == '\n') && last_letter != '\\') ) {
	    return token;
	}
	else if( tmp == '\'' && last_letter != '\\' ) {
	    // QUOTES: un-backslashed quotes => it's a string
	    token += _get_string(istr);
	    continue;
	}
	else if( tmp == '{' && last_letter == '$') {
	    token += '{' + _get_until_closing_bracket(istr);
	    continue;
	}
	else if( tmp == '[') {
	    token += '[' + _get_until_closing_square_bracket(istr);
	    continue;
	}
	else if( tmp == '$' && last_letter == '\\') {
	    token += getpot_cast_int<char>(tmp); tmp = 0;  //  so that last_letter will become = 0, not '$';
	    continue;
	}
	else if( tmp == '\\' && last_letter != '\\')
	    continue;              // don't append un-backslashed backslashes
	token += getpot_cast_int<char>(tmp);
    }
}

inline const std::string
GetPot::_get_string(std::istream& istr)
    // parse input until next matching '
{
    std::string str;
    int    tmp = 0;
    int    last_letter = 0;
    while(1 + 1 == 2) {
	last_letter = tmp; tmp = istr.get();
	if( tmp == EOF)  return str;
	// un-backslashed quotes => it's the end of the string
	else if( tmp == '\'' && last_letter != '\\')  return str;
	else if( tmp == '\\' && last_letter != '\\')  continue; // don't append

	str += getpot_cast_int<char>(tmp);
    }
}

inline const std::string
GetPot::_get_until_closing_bracket(std::istream& istr)
    // parse input until next matching }
{
    std::string str = "";
    int    tmp = 0;
    int    last_letter = 0;
    int    brackets = 1;
    while(1 + 1 == 2) {
	last_letter = tmp; tmp = istr.get();
	if( tmp == EOF) return str;
	else if( tmp == '{' && last_letter == '$') brackets += 1;
	else if( tmp == '}') {
	    brackets -= 1;
	    // un-backslashed brackets => it's the end of the string
	    if( brackets == 0) return str + '}';
	    else if( tmp == '\\' && last_letter != '\\')
		continue;  // do not append an unbackslashed backslash
	}
	str += getpot_cast_int<char>(tmp);
    }
}


inline const std::string
GetPot::_get_until_closing_square_bracket(std::istream& istr)
    // parse input until next matching ]
{
    std::string str = "";
    int    tmp = 0;
    int    brackets = 1;
    while(1 + 1 == 2) {
	tmp = istr.get();
	if( tmp == EOF) return str;
	else if( tmp == '[') {
            brackets += 1;
        }
	else if( tmp == ']') {
	    brackets -= 1;
	    if( brackets == 0) return str + ']';
	}
	str += getpot_cast_int<char>(tmp);
    }
}

inline std::string
GetPot::_process_section_label(const std::string& Section,
				STRING_VECTOR& section_stack)
{
    std::string sname = Section;
    //  1) subsection of actual section ('./' prefix)
    if( sname.length() >= 2 && sname.substr(0, 2) == "./" ) {
	sname = sname.substr(2);
    }
    //  2) subsection of parent section ('../' prefix)
    else if( sname.length() >= 3 && sname.substr(0, 3) == "../" ) {
	do {
	    if( section_stack.end() != section_stack.begin() )
		section_stack.pop_back();
	    sname = sname.substr(3);
	} while( sname.substr(0, 3) == "../" );
    }
    // 3) subsection of the root-section
    else {
	section_stack.erase(section_stack.begin(), section_stack.end());
	// [] => back to root section
    }

    if( sname != "" ) {
	// parse section name for 'slashes'
	unsigned i=0;
	while( i < sname.length() ) {
	    if( sname[i] == '/' ) {
		section_stack.push_back(sname.substr(0,i));
		if( i+1 < sname.length()-1 )
		    sname = sname.substr(i+1);
		i = 0;
	    }
	    else
		i++;
	}
	section_stack.push_back(sname);
    }
    std::string section_label = "";
    if( !section_stack.empty() ) {
	victorate(std::string, section_stack, it)
	    section_label += *it + "/";
    }
    return section_label;
}

// Use C++ istream/ostream to handle most type conversions.
template <typename T>
inline T
GetPot::_convert_to_type(const std::string& String, const T& Default) const
{
  std::istringstream in_string(String);
  T retval;
  in_string >> retval;
  if (in_string.fail())
    retval = Default;
  return retval;
}

// copy string - operator>> would have stopped upon seeing whitespace!
template <>
inline std::string
GetPot::_convert_to_type(const std::string& String, const std::string&) const
{
    return String;
}

// copy string
inline std::string
GetPot::_convert_to_type(const std::string& String, const char*) const
{
    return String;
}

// be more liberal than std C++ in what we interpret as a boolean
template<>
inline bool
GetPot::_convert_to_type<bool>(const std::string& String, const bool& Default) const
{
  std::string newstring(String);
  //std::transform(newstring.begin(), newstring.end(), newstring.begin(), std::toupper);
  for (unsigned int i=0; i<newstring.length(); ++i)
  {
    newstring[i]=getpot_cast_int<char>(toupper(newstring[i]));
  }

  // "true"/"True"/"TRUE" should work
  if (newstring.find("TRUE")!=std::string::npos)  return true;
  if (newstring.find("FALSE")!=std::string::npos) return false;

  // And if we don't find that, let's search for an integer and use C unsigned
  // int->bool conversion before giving up; i.e. a user could specify "0" for
  // false or "1" for true
  std::istringstream in_string(String);
  unsigned int retval;
  in_string >> retval;
  if (in_string.fail())
    return Default;

  return retval;
}

// Use C++ istream/ostream to handle most type conversions.
template <typename T>
inline T
GetPot::_convert_to_type_no_default(const char* VarName, const std::string& String, const T&) const
{
  std::istringstream in_string(String);
  T retval;
  in_string >> retval;
  if (in_string.fail())
  {
    getpot_cerr <<"ERROR: Input value for variable "<<VarName<<" is of the wrong type."<<std::endl;
    getpot_cerr <<"       value = "<<String<<" expected type = "<<typeid(T).name()<<std::endl;
    getpot_error();
  }
  return retval;
}

// copy string - operator>> would have stopped upon seeing whitespace!
template <>
inline std::string
GetPot::_convert_to_type_no_default(const char*, const std::string& String, const std::string&) const
{
    return String;
}

// copy string
inline std::string
GetPot::_convert_to_type_no_default(const char*, const std::string& String, const char*) const
{
    return String;
}

// be more liberal than std C++ in what we interpret as a boolean
template<>
inline bool
GetPot::_convert_to_type_no_default<bool>(const char* VarName, const std::string& String, const bool&) const
{
  std::string newstring(String);
  //std::transform(newstring.begin(), newstring.end(), newstring.begin(), std::toupper);
  for (unsigned int i=0; i<newstring.length(); ++i)
  {
    newstring[i]=getpot_cast_int<char>(toupper(newstring[i]));
  }

  // "true"/"True"/"TRUE" should work
  if (newstring.find("TRUE")!=std::string::npos)  return true;
  if (newstring.find("FALSE")!=std::string::npos) return false;

  // And if we don't find that, let's search for an integer and use C unsigned
  // int->bool conversion before giving up; i.e. a user could specify "0" for
  // false or "1" for true
  std::istringstream in_string(String);
  unsigned int retval;
  in_string >> retval;
  if (in_string.fail())
  {
    getpot_cerr <<"ERROR: Input value for variable "<<VarName<<" is of the wrong type."<<std::endl;
    getpot_cerr <<"       value = "<<String<<" expected type = "<<typeid(bool).name()<<std::endl;
    getpot_error();
  }

  return retval;
}

inline const char*
GetPot::_internal_managed_copy(const std::string& Arg) const
{
    const char* arg = Arg.c_str();

    // Get a lock before touching anything mutable
    SCOPED_MUTEX;

    // See if there's already an identical string saved
    std::set<const char*,ltstr>::const_iterator it =
        _internal_string_container.find(arg);

    // If so, return it
    if (it != _internal_string_container.end())
        return *it;

    // Otherwise, create a new one
    char* newcopy = new char[strlen(arg)+1];
    strncpy(newcopy, arg, strlen(arg)+1);
    _internal_string_container.insert(newcopy);
    return newcopy;
}

//////////////////////////////////////////////////////////////////////////////
// (*) cursor oriented functions
//.............................................................................
inline const std::string
GetPot::_get_remaining_string(const std::string& String, const std::string& Start) const
    // Checks if 'String' begins with 'Start' and returns the remaining String.
    // Returns None if String does not begin with Start.
{
    if( Start == "" ) return String;
    // note: java.lang.String: substring(a,b) = from a to b-1
    //        C++ string:      substr(a,b)    = from a to a + b
    if( String.find(Start) == 0 ) return String.substr(Start.length());
    else                          return "";
}

//     -- search for a certain argument and set cursor to position
inline bool
GetPot::search(const std::string &Option)
{
  return search(Option.c_str());
}

//     -- search for a certain argument and set cursor to position
inline bool
GetPot::search(const char* Option)
{
    unsigned           OldCursor = cursor;
    const std::string  SearchTerm = prefix + Option;

    // (*) record requested arguments for later ufo detection
    _record_argument_request(SearchTerm);

    if( OldCursor >= argv.size() )
      OldCursor = getpot_cast_int<unsigned>(argv.size() - 1);
    search_failed_f = true;

    // (*) first loop from cursor position until end
    unsigned  c = cursor;
    for(; c < argv.size(); c++) {
	if( argv[c] == SearchTerm )
	{ cursor = c; search_failed_f = false; return true; }
    }
    if( ! search_loop_f ) return false;

    // (*) second loop from 0 to old cursor position
    for(c = 1; c < OldCursor; c++) {
	if( argv[c] == SearchTerm )
	{ cursor = c; search_failed_f = false; return true; }
    }
    // in case nothing is found the cursor stays where it was
    return false;
}


inline bool
GetPot::search(unsigned No, const char* P, ...)
{
    // (*) recording the requested arguments happens in subroutine 'search'
    if( No == 0 ) return false;

    // search for the first argument
    if( search(P) == true ) return true;

    // start interpreting variable argument list
    va_list ap;
    va_start(ap, P);
    unsigned i = 1;
    for(; i < No; i++) {
	char* Opt = va_arg(ap, char *);
	// (*) search records itself for later ufo detection
	if( search(Opt) == true ) break;
    }

    if( i < No ) {
	i++;
	// loop was left before end of array --> hit but
	// make sure that the rest of the search terms is marked
	// as requested.
	for(; i < No; i++) {
	    char* Opt = va_arg(ap, char *);
	    // (*) record requested arguments for later ufo detection
	    _record_argument_request(Opt);
	}
	va_end(ap);
	return true;
    }

    va_end(ap);
    // loop was left normally --> no hit
    return false;
}

inline void
GetPot::reset_cursor()
{ search_failed_f = false; cursor = 0; }

inline void
GetPot::init_multiple_occurrence()
{ disable_loop(); reset_cursor(); }
///////////////////////////////////////////////////////////////////////////////
// (*) direct access to command line arguments
//.............................................................................
//
inline const char*
GetPot::operator[](unsigned idx) const
{ return idx<argv.size() ? argv[idx].c_str() : 0; }

template <typename T>
inline T
GetPot::get(unsigned int Idx, const T& Default) const
{
    if( Idx >= argv.size() ) return Default;
    return _convert_to_type(argv[Idx], Default);
}

inline const char*
GetPot::get(unsigned int Idx, const char* Default) const
{
    if( Idx >= argv.size() ) return Default;
    return argv[Idx].c_str();
}

inline unsigned
GetPot::size() const
{ return getpot_cast_int<unsigned>(argv.size()); }


//     -- next() function group
template <typename T>
inline T
GetPot::next(const T& Default)
{
    if( search_failed_f ) return Default;
    cursor++;
    if( cursor >= argv.size() )
    { cursor = getpot_cast_int<unsigned>(argv.size()); return Default; }

    // (*) record requested argument for later ufo detection
    _record_argument_request(argv[cursor]);

    const std::string Remain = _get_remaining_string(argv[cursor], prefix);

    return Remain != "" ? _convert_to_type(Remain, Default) : Default;
}

inline const char*
GetPot::next(const char* Default)
{
    return _internal_managed_copy(next(std::string(Default)));
}

//     -- follow() function group
//        distinct option to be searched for
template <typename T>
inline T
GetPot::follow(const T& Default, const char* Option)
{
    // (*) record requested of argument is entirely handled in 'search()' and 'next()'
    if( search(Option) == false ) return Default;
    return next(Default);
}

inline const char*
GetPot::follow(const char* Default, const char* Option)
{
    return _internal_managed_copy(follow(std::string(Default), Option));
}

//     -- second follow() function group
//        multiple option to be searched for
template <typename T>
inline T
GetPot::follow(const T& Default, unsigned int No, const char* P, ...)
{
    // (*) record requested of argument is entirely handled in 'search()' and 'next()'
    if( No == 0 ) return Default;
    if( search(P) == true ) return next(Default);

    va_list ap;
    va_start(ap, P);
    unsigned i=1;
    for(; i<No; i++) {
	char* Opt = va_arg(ap, char *);
	if( search(Opt) == true ) {
	    va_end(ap);
	    return next(Default);
	}
    }
    va_end(ap);
    return Default;
}

inline const char*
GetPot::follow(const char* Default, unsigned No, const char* P, ...)
{
    // (*) record requested of argument is entirely handled in 'search()' and 'next()'
    if( No == 0 ) return Default;
    if( search(P) == true ) return next(Default);

    va_list ap;
    va_start(ap, P);
    unsigned i=1;
    for(; i<No; i++) {
	char* Opt = va_arg(ap, char *);
	if( search(Opt) == true ) {
	    va_end(ap);
	    return next(Default);
	}
    }
    va_end(ap);
    return Default;
}

///////////////////////////////////////////////////////////////////////////////
// (*) directly connected options
//.............................................................................
//
template <typename T>
inline T
GetPot::direct_follow(const T& Default, const char* Option)
{
    const char* FollowStr = _match_starting_string(Option);

    // (*) record requested of argument for later ufo-detection
    _record_argument_request(std::string(Option) + FollowStr);

    if( FollowStr == 0 )                    return Default;
    if( ++cursor >= argv.size() )
      cursor = getpot_cast_int<unsigned>(argv.size());
    return _convert_to_type(FollowStr, Default);
}

inline const char*
GetPot::direct_follow(const char* Default, const char* Option)
{
    return _internal_managed_copy(direct_follow(std::string(Default), Option));
}

inline const char*
GetPot::_match_starting_string(const char* StartString)
    // pointer  to the place where the string after
    //          the match inside the found argument starts.
    // 0        no argument matches the starting string.
{
    const unsigned N =
      getpot_cast_int<unsigned>(strlen(StartString));
    unsigned       OldCursor = cursor;

    if( OldCursor >= argv.size() )
      OldCursor = getpot_cast_int<unsigned>(argv.size() - 1);
    search_failed_f = true;

    // (*) first loop from cursor position until end
    unsigned c = cursor;
    for(; c < argv.size(); c++) {
	if( strncmp(StartString, argv[c].c_str(), N) == 0)
	{ cursor = c; search_failed_f = false; return &(argv[c].c_str()[N]); }
    }

    if( ! search_loop_f ) return NULL;

    // (*) second loop from 0 to old cursor position
    for(c = 1; c < OldCursor; c++) {
	if( strncmp(StartString, argv[c].c_str(), N) == 0)
	{ cursor = c; search_failed_f = false; return &(argv[c].c_str()[N]); }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// (*) search for flags
//.............................................................................
//
inline bool
GetPot::options_contain(const char* FlagList) const
{
    // go through all arguments that start with a '-' (but not '--')
    std::string str;
    STRING_VECTOR::const_iterator it = argv.begin();
    for(; it != argv.end(); ++it) {
	str = _get_remaining_string(*it, prefix);

	if( str.length() >= 2 && str[0] == '-' && str[1] != '-' )
	    if( _check_flags(str, FlagList) ) return true;
    }
    return false;
}

inline bool
GetPot::argument_contains(unsigned Idx, const char* FlagList) const
{
    if( Idx >= argv.size() ) return false;

    // (*) record requested of argument for later ufo-detection
    //     an argument that is checked for flags is considered to be 'requested'
    _record_argument_request(argv[Idx]);

    if( prefix == "" )
	// search argument for any flag in flag list
	return _check_flags(argv[Idx], FlagList);

    // if a prefix is set, then the argument index is the index
    //   inside the 'namespace'
    // => only check list of arguments that start with prefix
    unsigned no_matches = 0;
    unsigned i=0;
    for(; i<argv.size(); i++) {
	const std::string Remain = _get_remaining_string(argv[i], prefix);
	if( Remain != "") {
	    no_matches += 1;
	    if( no_matches == Idx)
		return _check_flags(Remain, FlagList);
	}
    }
    // no argument in this namespace
    return false;
}

inline bool
GetPot::_check_flags(const std::string& Str, const char* FlagList) const
{
    const char* p=FlagList;
    for(; *p != '\0' ; p++)
	if( Str.find(*p) != std::string::npos ) return true; // found something
    return false;
}

///////////////////////////////////////////////////////////////////////////////
// (*) nominus arguments
inline STRING_VECTOR
GetPot::nominus_vector() const
    // return vector of nominus arguments
{
    STRING_VECTOR nv;
    std::vector<unsigned>::const_iterator it = idx_nominus.begin();
    for(; it != idx_nominus.end(); ++it) {
	nv.push_back(argv[*it]);

	// (*) record for later ufo-detection
	//     when a nominus vector is requested, the entire set of nominus arguments are
	//     tagged as 'requested'
	_record_argument_request(argv[*it]);
    }
    return nv;
}

inline const char*
GetPot::next_nominus()
{
    if( nominus_cursor < int(idx_nominus.size()) - 1 ) {
	const std::string Tmp = argv[idx_nominus[++nominus_cursor]];

	// (*) record for later ufo-detection
	_record_argument_request(Tmp);

	return Tmp.c_str();
    }
    return 0;
}

inline void
GetPot::reset_nominus_cursor()
{ nominus_cursor = -1; }
///////////////////////////////////////////////////////////////////////////////
// (*) variables
//.............................................................................
//
inline bool
GetPot::have_variable(const char* VarName) const
{
    const variable* sv = _request_variable(VarName);
    if (sv == 0) return false;
    return true;
}

inline bool
GetPot::have_variable(const std::string& VarName) const
{
    return have_variable(VarName.c_str());
}

template <typename T>
inline T
GetPot::operator()(const char* VarName, const T& Default) const
{
    // (*) recording of requested variables happens in '_request_variable()'
    const variable*  sv = _request_variable(VarName);
    if( sv == 0 ) return Default;
    return _convert_to_type(sv->original, Default);
}

template <typename T>
inline T
GetPot::operator()(const std::string& VarName, const T& Default) const
{
    return operator()(VarName.c_str(), Default);
}

inline const char*
GetPot::operator()(const char* VarName, const char* Default) const
{
    return _internal_managed_copy(operator()(VarName, std::string(Default)));
}

inline const char*
GetPot::operator()(const std::string& VarName, const char* Default) const
{
    return operator()(VarName.c_str(), Default);
}

template <typename T>
inline T
GetPot::operator()(const char* VarName, const T& Default, unsigned int Idx) const
{
    // (*) recording of requested variables happens in '_request_variable()'
    const variable* sv = _request_variable(VarName);
    if( sv == 0 ) return Default;
    const std::string*  element = sv->get_element(Idx);
    if( element == 0 ) return Default;
    return _convert_to_type(*element, Default);
}

template <typename T>
inline T
GetPot::operator()(const std::string& VarName, const T& Default, unsigned int Idx) const
{
    return operator()(VarName.c_str(), Default, Idx);
}

inline const char*
GetPot::operator()(const char* VarName, const char* Default, unsigned int Idx) const
{
    return _internal_managed_copy(operator()(VarName, std::string(Default), Idx));
}

inline const char*
GetPot::operator()(const std::string& VarName, const char* Default, unsigned int Idx) const
{
    return operator()(VarName.c_str(), Default, Idx);
}

template <typename T>
inline T
GetPot::get_value_no_default(const char* VarName, const T& Default ) const
{
    // (*) recording of requested variables happens in '_request_variable()'
    const variable*  sv = _request_variable(VarName);
    if( sv == 0 )
    {
      getpot_cerr << "ERROR: cannot find variable "<<VarName<<std::endl;
      getpot_error();
    }
    return _convert_to_type_no_default(VarName, sv->original, Default);
}

template <typename T>
inline T
GetPot::get_value_no_default(const std::string& VarName, const T& Default ) const
{
    return get_value_no_default(VarName.c_str(),Default);
}

inline const char*
GetPot::get_value_no_default(const char* VarName, const char* Default) const
{
    return _internal_managed_copy(get_value_no_default(VarName, Default));
}

inline const char*
GetPot::get_value_no_default(const std::string& VarName, const char* Default ) const
{
    return get_value_no_default(VarName.c_str(),Default);
}

template <typename T>
inline T
GetPot::get_value_no_default(const char* VarName, const T& Default, unsigned int Idx) const
{
    // (*) recording of requested variables happens in '_request_variable()'
    const variable* sv = _request_variable(VarName);
    if( sv == 0 ){
      getpot_cerr << "ERROR: cannot find variable "<<VarName<<std::endl;
      getpot_error();
    }
    const std::string*  element = sv->get_element(Idx);
    if( element == 0 ){
      getpot_cerr << "ERROR: cannot find index "<<Idx<<" of variable "<<VarName<<std::endl;
      getpot_error();
    }
    return _convert_to_type_no_default(VarName, *element, Default);
}

template <typename T>
inline T
GetPot::get_value_no_default(const std::string& VarName, const T& Default, unsigned int Idx) const
{
    return get_value_no_default(VarName.c_str(), Default, Idx);
}

inline const char*
GetPot::get_value_no_default(const char* VarName, const char* Default, unsigned int Idx) const
{
    return _internal_managed_copy(get_value_no_default(VarName, std::string(Default), Idx));
}

inline const char*
GetPot::get_value_no_default(const std::string& VarName, const char* Default, unsigned int Idx) const
{
    return get_value_no_default(VarName.c_str(), Default, Idx);
}

inline void
GetPot::_record_argument_request(const std::string& Name) const
{
    if( ! request_recording_f ) return;

    // Get a lock before touching anything mutable
    SCOPED_MUTEX;

    // (*) record requested variable for later ufo detection
    _requested_arguments.insert(Name);

    // (*) record considered section for ufo detection
    STRING_VECTOR      STree = _get_section_tree(Name);
    victorate(std::string, STree, it)
	if( _requested_sections.find(*it) == _requested_sections.end() )
	    if( section.length() != 0 ) _requested_sections.insert(*it);
}

inline void
GetPot::_record_variable_request(const std::string& Name) const
{
    if( ! request_recording_f ) return;

    // Get a lock before touching anything mutable
    SCOPED_MUTEX;

    // (*) record requested variable for later ufo detection
    _requested_variables.insert(Name);

    // (*) record considered section for ufo detection
    STRING_VECTOR      STree = _get_section_tree(Name);
    victorate(std::string, STree, it)
	if( _requested_sections.find(*it) == _requested_sections.end() )
	    if( section.length() != 0 ) _requested_sections.insert(*it);
}

// (*) following functions are to be used from 'outside', after getpot has parsed its
//     arguments => append an argument in the argument vector that reflects the addition
inline void
GetPot::_set_variable(const std::string& VarName,
		      const std::string& Value, const bool Requested /* = true */)
{
    const GetPot::variable* Var = Requested ?
		    _request_variable(VarName.c_str()) :
		    _find_variable(VarName.c_str());
    if( Var == 0 ) variables.push_back(variable(VarName.c_str(),
						Value.c_str(), _field_separator.c_str()));
    else {
      overridden_vars.insert(VarName.c_str());
      (const_cast<GetPot::variable*>(Var))->take(Value.c_str(), _field_separator.c_str());
    }
}

template <typename T>
inline void
GetPot::set(const char* VarName, const T& Value, const bool Requested /* = true */)
{
  std::ostringstream string_value;
  string_value << Value;
  _set_variable(VarName, string_value.str().c_str(), Requested);
}

template <typename T>
inline void
GetPot::set(const std::string& VarName, const T& Value, const bool Requested /* = true */)
{
    set(VarName.c_str(), Value, Requested);
}

inline void
GetPot::set(const char* VarName, const char* Value, const bool Requested /* = true */)
{
  _set_variable(VarName, Value, Requested);
}

inline void
GetPot::set(const std::string& VarName, const char* Value, const bool Requested /* = true */)
{
    set(VarName.c_str(), Value, Requested);
}

inline unsigned
GetPot::vector_variable_size(const char* VarName) const
{
    const variable*  sv = _request_variable(VarName);
    if( sv == 0 ) return 0;
    return (unsigned)(sv->value.size());
}

inline unsigned
GetPot::vector_variable_size(const std::string& VarName) const
{
    return vector_variable_size(VarName.c_str());
}

inline STRING_VECTOR
GetPot::get_variable_names() const
{
    STRING_VECTOR  result;
    std::vector<GetPot::variable>::const_iterator it = variables.begin();
    for(; it != variables.end(); ++it) {
	const std::string Tmp = _get_remaining_string((*it).name, prefix);
	if( Tmp != "" ) result.push_back(Tmp);
    }
    return result;
}

inline STRING_VECTOR
GetPot::get_section_names() const
{ return section_list; }

inline std::set<std::string>
GetPot::get_overridden_variables() const
{ return overridden_vars; }

inline const GetPot::variable*
GetPot::_find_variable(const char* VarName) const
{
    const std::string Name = prefix + VarName;

    std::vector<variable>::const_iterator it = variables.begin();
    for(; it != variables.end(); ++it) {
	if( (*it).name == Name ) return &(*it);
    }
    return 0;
}

inline const GetPot::variable*
GetPot::_request_variable(const char* VarName) const
{
    // (*) record requested variable for later ufo detection
    this->_record_variable_request(VarName);

    return this->_find_variable(VarName);
}

///////////////////////////////////////////////////////////////////////////////
// (*) ouput (basically for debugging reasons
//.............................................................................
//
inline int
GetPot::print(std::ostream &out_stream) const
{
    out_stream << "argc = " << argv.size() << std::endl;
    STRING_VECTOR::const_iterator it = argv.begin();
    for(; it != argv.end(); ++it)
	out_stream << *it << std::endl;
    out_stream << std::endl;
    return 1;
}

// PECOS/HPCT Addition - add option to prepend output with a delimiter
// while also disabling argc print and skipping first print (the name
// of the input file)
//
// PECOS Development Team: (ks. 4/16/09)

inline int
GetPot::print(const char* custom_prefix, std::ostream &out_stream, unsigned int skip_count) const
{
    STRING_VECTOR::const_iterator it = argv.begin();
    it += skip_count;
    for(; it != argv.end(); ++it)
      {
	out_stream << custom_prefix;
        out_stream << *it << std::endl;
      }
    out_stream << std::endl;
    return 1;
}


// (*) dollar bracket expressions (DBEs) ------------------------------------
//
//     1) Entry Function: _DBE_expand_string()
//        Takes a string such as
//
//          "${+ ${x} ${y}}   Subject-${& ${section} ${subsection}}:   ${title}"
//
//        calls _DBE_expand() for each of the expressions
//
//           ${+ ${x} ${y}}
//           ${& ${section} ${subsection}}
//           ${Title}
//
//        and returns the string
//
//          "4711 Subject-1.01:   Mit den Clowns kamen die Schwaene"
//
//        assuming that
//            x          = "4699"
//            y          = "12"
//            section    = "1."
//            subsection = "01"
//            title      = "Mit den Clowns kamen die Schwaene"
//
//      2) _DBE_expand():
//
//           checks for the command, i.e. the 'sign' that follows '${'
//           divides the argument list into sub-expressions using
//           _DBE_get_expr_list()
//
//           ${+ ${x} ${y}}                 -> "${x}"  "${y}"
//           ${& ${section} ${subsection}}  -> "${section}" "${subsection}"
//           ${Title}                       -> Nothing, variable expansion
//
//      3) _DBE_expression_list():
//
//           builds a vector of unbracketed whitespace separated strings, i.e.
//
//           "  ${Number}.a ${: Das Marmorbild} AB-${& Author= ${Eichendorf}-1870}"
//
//           is split into a vector
//
//              [0] ${Number}.a
//              [1] ${: Das Marmorbild}
//              [2] AB-${& Author= ${Eichendorf}}-1870
//
//           Each sub-expression is expanded using expand().
//---------------------------------------------------------------------------
inline std::string
GetPot::_DBE_expand_string(const std::string& str)
{
    // Parses for closing operators '${ }' and expands them letting
    // white spaces and other letters as they are.
    std::string   new_string = "";
    unsigned open_brackets = 0;
    unsigned first = 0;
    unsigned i = 0;
    for(;  i<str.size(); i++) {
	if( i < str.size() - 2 && str.substr(i, 2) == "${" ) {
	    if( open_brackets == 0 ) first = i+2;
	    open_brackets++;
	}
	else if( str[i] == '}' && open_brackets > 0) {
	    open_brackets -= 1;
	    if( open_brackets == 0 ) {
		const std::string Replacement = _DBE_expand(str.substr(first, i - first));
		new_string += Replacement;
	    }
	}
	else if( open_brackets == 0 )
	    new_string += str[i];
    }
    return new_string;
}

inline STRING_VECTOR
GetPot::_DBE_get_expr_list(const std::string& str_, const unsigned ExpectedNumber)
    // ensures that the resulting vector has the expected number
    // of arguments, but they may contain an error message
{
    std::string str = str_;
    // Separates expressions by non-bracketed whitespaces, expands them
    // and puts them into a list.

    unsigned i=0;
    // (1) eat initial whitespaces
    for(; i < str.size(); i++)
	if( ! isspace(str[i]) ) break;

    STRING_VECTOR   expr_list;
    unsigned         open_brackets = 0;
    std::vector<unsigned> start_idx;
    unsigned         start_new_string = i;
    unsigned         l = (unsigned)(str.size());

    // (2) search for ${ } expressions ...
    while( i < l ) {
	const char letter = str[i];
	// whitespace -> end of expression
	if( isspace(letter) && open_brackets == 0) {
	    expr_list.push_back(str.substr(start_new_string, i - start_new_string));
	    bool no_breakout_f = true;
	    for(i++; i < l ; i++) {
		if( ! isspace(str[i]) )
		{ no_breakout_f = false; start_new_string = i; break; }
	    }
	    if( no_breakout_f ) {
		// end of expression list
		if( expr_list.size() < ExpectedNumber ) {
		    const std::string   pre_tmp("<< ${ }: missing arguments>>");
		    STRING_VECTOR tmp(ExpectedNumber - expr_list.size(), pre_tmp);
		    expr_list.insert(expr_list.end(), tmp.begin(), tmp.end());
		}
		return expr_list;
	    }
	}

	// dollar-bracket expression
	if( str.length() >= i+2 && str.substr(i, 2) == "${" ) {
	    open_brackets++;
	    start_idx.push_back(i+2);
	}
	else if( letter == '}' && open_brackets > 0) {
	    int start = start_idx[start_idx.size()-1];
	    start_idx.pop_back();
	    const std::string Replacement = _DBE_expand(str.substr(start, i-start));
	    if( start - 3 < (int)0)
		str = Replacement + str.substr(i+1);
	    else
		str = str.substr(0, start-2) + Replacement + str.substr(i+1);
	    l = (int)(str.size());
	    i = start + (int)(Replacement.size()) - 3;
	    open_brackets--;
	}
	i++;
    }

    // end of expression list
    expr_list.push_back(str.substr(start_new_string, i-start_new_string));

    if( expr_list.size() < ExpectedNumber ) {
	const std::string   pre_tmp("<< ${ }: missing arguments>>");
	STRING_VECTOR tmp(ExpectedNumber - expr_list.size(), pre_tmp);
	expr_list.insert(expr_list.end(), tmp.begin(), tmp.end());
    }

    return expr_list;
}

inline const GetPot::variable*
GetPot::_DBE_get_variable(const std::string& VarName)
{
    static GetPot::variable ev;
    std::string secure_Prefix = prefix;

    prefix = section;
    // (1) first search in currently active section
    const GetPot::variable* var = _request_variable(VarName.c_str());
    if( var != 0 ) { prefix = secure_Prefix; return var; }

    // (2) search in root name space
    prefix = "";
    var = _request_variable(VarName.c_str());
    if( var != 0 ) { prefix = secure_Prefix; return var; }

    prefix = secure_Prefix;

    // error occured => variable name == ""
    ev.original = "<<${ } variable '";
    ev.original += VarName + "' undefined>>";
    return &ev;
}

inline std::string
GetPot::_DBE_expand(const std::string& expr)
{
    // ${: } pure text
    if( expr[0] == ':' )
	return expr.substr(1);

    // ${& expr expr ... } text concatination
    else if( expr[0] == '&' ) {
	const STRING_VECTOR A = _DBE_get_expr_list(expr.substr(1), 1);

	STRING_VECTOR::const_iterator it = A.begin();
	std::string result = *it++;
	for(; it != A.end(); ++it) result += *it;

	return result;
    }

    // ${<-> expr expr expr} text replacement
    else if( expr.length() >= 3 && expr.substr(0, 3) == "<->" ) {
	STRING_VECTOR A = _DBE_get_expr_list(expr.substr(3), 3);
	size_t tmp = 0;
	const size_t L = A[1].length();
	while( (tmp = A[0].find(A[1])) != std::string::npos ) {
	    A[0].replace(tmp, L, A[2]);
	}
	return A[0];
    }

    // ${=func [expr...] } function evaluation
    else if( expr.length() >= 2 &&
	     expr.substr(0, 1) == "=" &&
	     expr.substr(0, 2) != "==" ) {
	size_t funcnamestart = expr.find_first_not_of(" \t", 1);
	if (funcnamestart != std::string::npos) {
            size_t funcnameend = expr.find_first_of(" \t",funcnamestart);
	    std::string funcname = expr.substr(funcnamestart,
					       funcnameend-funcnamestart);
	    if (funcname == "log") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::log(arg));
	    }
	    else if (funcname == "log10") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::log10(arg));
	    }
	    else if (funcname == "sin") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::sin(arg));
	    }
	    else if (funcname == "cos") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::cos(arg));
	    }
	    else if (funcname == "tan") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::tan(arg));
	    }
	    else if (funcname == "asin") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::asin(arg));
	    }
	    else if (funcname == "acos") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::acos(arg));
	    }
	    else if (funcname == "atan") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::atan(arg));
	    }
	    else if (funcname == "atan2") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 2);
	        double arg1 = _convert_to_type(A[0], 0.0);
	        double arg2 = _convert_to_type(A[1], 0.0);
	        return _convert_from_type(std::atan2(arg1, arg2));
	    }
	    else if (funcname == "sinh") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::sinh(arg));
	    }
	    else if (funcname == "cosh") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::cosh(arg));
	    }
	    else if (funcname == "tanh") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::tanh(arg));
	    }
	    else if (funcname == "sqrt") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::sqrt(arg));
	    }
	    else if (funcname == "abs") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::abs(arg));
	    }
	    else if (funcname == "max") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        STRING_VECTOR::const_iterator it = A.begin();
	        double result = _convert_to_type(*it++, 0.0);
	        for(; it != A.end(); ++it)
	            result = std::max(result, _convert_to_type(*it, 0.0));
	        return _convert_from_type(result);
	    }
	    else if (funcname == "min") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        STRING_VECTOR::const_iterator it = A.begin();
	        double result = _convert_to_type(*it++, 0.0);
	        for(; it != A.end(); ++it)
	            result = std::min(result, _convert_to_type(*it, 0.0));
	        return _convert_from_type(result);
	    }
	    else if (funcname == "ceil") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::ceil(arg));
	    }
	    else if (funcname == "floor") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        double arg = _convert_to_type(A[0], 0.0);
	        return _convert_from_type(std::floor(arg));
	    }
	    else if (funcname == "fmod") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 2);
	        double arg1 = _convert_to_type(A[0], 0.0);
	        double arg2 = _convert_to_type(A[1], 0.0);
	        return _convert_from_type(std::fmod(arg1, arg2));
	    }
	    else if (funcname == "srand") {
	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
	        unsigned int arg = _convert_to_type(A[0], 0u);
		std::srand(arg);
	        return A[0];
	    }
	    // ${=rand range} with default range==RAND_MAX
	    else if (funcname == "rand") {
                if (funcnameend >= expr.length() ||
                    expr.find_first_not_of(" \t", funcnameend) == std::string::npos)
	          return _convert_from_type(std::rand());

	        STRING_VECTOR A =
                  _DBE_get_expr_list(expr.substr(funcnameend), 1);
		unsigned int range = _convert_to_type(A[0],0u);
		if (!range)
	          return _convert_from_type(0);
		const unsigned int x = (RAND_MAX + 1u) / range;
		const unsigned int y = x * range;
		unsigned int returnval;
		do {
			returnval = rand();
		} while(returnval >= y);
		return _convert_from_type(returnval / x);
	    }
	    else if (funcname == "time") {
	        return _convert_from_type(std::time(NULL));
	    }
	}
    }

    // ${+ ...}, ${- ...}, ${* ...}, ${/ ...} expressions
    else if( expr[0] == '+' ) {
	STRING_VECTOR A = _DBE_get_expr_list(expr.substr(1), 2);
	STRING_VECTOR::const_iterator it = A.begin();
	double result = _convert_to_type(*it++, 0.0);
	for(; it != A.end(); ++it)
	    result += _convert_to_type(*it, 0.0);

	return _convert_from_type(result);
    }
    else if( expr[0] == '-' ) {
	STRING_VECTOR A = _DBE_get_expr_list(expr.substr(1), 2);
	STRING_VECTOR::const_iterator it = A.begin();
	double result = _convert_to_type(*it++, 0.0);
	for(; it != A.end(); ++it)
	    result -= _convert_to_type(*it, 0.0);

	return _convert_from_type(result);
    }
    else if( expr[0] == '*' ) {
	STRING_VECTOR A = _DBE_get_expr_list(expr.substr(1), 2);
	STRING_VECTOR::const_iterator it = A.begin();
	double result = _convert_to_type(*it++, 0.0);
	for(; it != A.end(); ++it)
	    result *= _convert_to_type(*it, 0.0);

	return _convert_from_type(result);
    }
    else if( expr[0] == '/' ) {
	STRING_VECTOR A = _DBE_get_expr_list(expr.substr(1), 2);
	STRING_VECTOR::const_iterator it = A.begin();
	double result = _convert_to_type(*it++, 0.0);
	if( result == 0 ) {
            return "0.0";
        }
	for(; it != A.end(); ++it) {
	    const double Q = _convert_to_type(*it, 0.0);
	    result /= Q;
	}
	return _convert_from_type(result);
    }

    // ${^ ... } power expressions
    else if( expr[0] == '^' ) {
	STRING_VECTOR A = _DBE_get_expr_list(expr.substr(1), 2);
	STRING_VECTOR::const_iterator it = A.begin();
	double result = _convert_to_type(*it++, 0.0);
	for(; it != A.end(); ++it)
	    result = pow(result, _convert_to_type(*it, 0.0));
	return _convert_from_type(result);
    }

    // ${==  } ${<=  } ${>= } comparisons (return the number of the first 'match'
    else if( expr.length() >= 2 &&
	     ( expr.substr(0,2) == "==" || expr.substr(0,2) == ">=" ||
	       expr.substr(0,2) == "<=" || expr[0] == '>'           || expr[0] == '<')) {
	// differentiate between two and one sign operators
	unsigned op = 0;
	enum { EQ, GEQ, LEQ, GT, LT };
	if      ( expr.substr(0, 2) == "==" ) op = EQ;
	else if ( expr.substr(0, 2) == ">=" ) op = GEQ;
	else if ( expr.substr(0, 2) == "<=" ) op = LEQ;
	else if ( expr[0] == '>' )            op = GT;
	else    /*                     "<" */ op = LT;

	STRING_VECTOR a;
	if ( op == GT || op == LT ) a = _DBE_get_expr_list(expr.substr(1), 2);
	else                        a = _DBE_get_expr_list(expr.substr(2), 2);

	std::string   x_orig = a[0];
	double   x = _convert_to_type(x_orig, 1e37);
	unsigned i = 1;

	STRING_VECTOR::const_iterator y_orig = a.begin();
	for(y_orig++; y_orig != a.end(); ++y_orig) {
	    double y = _convert_to_type(*y_orig, 1e37);

	    // set the strings as reference if one wasn't a number
	    if ( x == 1e37 || y == 1e37 ) {
		// it's a string comparison
		if( (op == EQ  && x_orig == *y_orig) || (op == GEQ && x_orig >= *y_orig) ||
		    (op == LEQ && x_orig <= *y_orig) || (op == GT  && x_orig >  *y_orig) ||
		    (op == LT  && x_orig <  *y_orig) )
		    return _convert_from_type(i);
	    }
	    else {
		// it's a number comparison
		if( (op == EQ  && x == y) || (op == GEQ && x >= y) ||
		    (op == LEQ && x <= y) || (op == GT  && x >  y) ||
		    (op == LT  && x <  y) )
		    return _convert_from_type(i);
	    }
	    i++;
	}

	// nothing fulfills the condition => return 0
	return "0";
    }
    // ${?? expr expr} select
    else if( expr.length() >= 2 && expr.substr(0, 2) == "??" ) {
	STRING_VECTOR a = _DBE_get_expr_list(expr.substr(2), 2);
	double x = _convert_to_type(a[0], 1e37);
	// last element is always the default argument
	if( x == 1e37 || x < 0 || x >= a.size() - 1 ) return a[a.size()-1];

	// round x to closest integer
	return a[int(x+0.5)];
    }
    // ${? expr expr expr} if then else conditions
    else if( expr[0] == '?' ) {
	STRING_VECTOR a = _DBE_get_expr_list(expr.substr(1), 2);
	if( _convert_to_type(a[0], 0.0) == 1.0 ) return a[1];
	else if( a.size() > 2 ) return a[2];
    }
    // ${! expr} maxro expansion
    else if( expr[0] == '!' ) {
	const GetPot::variable* Var = _DBE_get_variable(expr.substr(1));
	// error
	if( Var->name == "" ) return std::string(Var->original);

	const STRING_VECTOR A = _DBE_get_expr_list(Var->original, 2);
	return A[0];
    }
    // ${@: } - string subscription
    else if( expr.length() >= 2 && expr.substr(0,2) == "@:" ) {
	const STRING_VECTOR A = _DBE_get_expr_list(expr.substr(2), 2);
	double x = _convert_to_type(A[1], 1e37);

	// last element is always the default argument
	if( x == 1e37 || x < 0 || x >= A[0].size() - 1)
	    return "<<1st index out of range>>";

	if( A.size() > 2 ) {
	    double y = _convert_to_type(A[2], 1e37);
	    if ( y != 1e37 && y > 0 && y <= A[0].size() - 1 && y > x )
		return A[0].substr(int(x+0.5), int(y+1.5) - int(x+0.5));
	    else if( y == -1 )
		return A[0].substr(int(x+0.5));
	    return "<<2nd index out of range>>";
	}
	else {
	    char* tmp = new char[2];
	    tmp[0] = A[0][int(x+0.5)]; tmp[1] = '\0';
	    std::string result(tmp);
	    delete [] tmp;
	    return result;
	}
    }
    // ${@ } - vector subscription
    else if( expr[0] == '@' ) {
	STRING_VECTOR          A   = _DBE_get_expr_list(expr.substr(1), 2);
	const GetPot::variable* Var = _DBE_get_variable(A[0]);
	// error
	if( Var->name == "" ) {
	    // make a copy of the string if an error occured
	    // (since the error variable is a static variable inside get_variable())
	    return std::string(Var->original);
	}

	double x = _convert_to_type(A[1], 1e37);

	// last element is always the default argument
	if (x == 1e37 || x < 0 || x >= Var->value.size() )
	    return "<<1st index out of range>>";

	if ( A.size() > 2) {
	    double y = _convert_to_type(A[2], 1e37);
	    int    begin = int(x+0.5);
	    int    end = 0;
	    if ( y != 1e37 && y > 0 && y <= Var->value.size() && y > x)
		end = int(y+1.5);
	    else if( y == -1 )
		end = int(Var->value.size());
	    else
		return "<<2nd index out of range>>";

	    std::string result = *(Var->get_element(begin));
	    int i = begin+1;
	    for(; i < end; i++)
		result += std::string(" ") + *(Var->get_element(i));
	    return result;
	}
	else
	    return *(Var->get_element(int(x+0.5)));
    }

    const STRING_VECTOR    A = _DBE_get_expr_list(expr, 1);
    const GetPot::variable* B = _DBE_get_variable(A[0]);

    // make a copy of the string if an error occured
    // (since the error variable is a static variable inside get_variable())
    if( B->name == "" ) return std::string(B->original);
    // (psuggs@pobox.com mentioned to me the warning MSVC++6.0 produces
    //  with:  else return B->original (thanks))
    return B->original;
}


///////////////////////////////////////////////////////////////////////////////
// (*) unidentified flying objects
//.............................................................................
//
inline bool
GetPot::_search_string_vector(const STRING_VECTOR& VecStr, const std::string& Str) const
{
    victorate(std::string, VecStr, itk) {
	if( *itk == Str ) return true;
    }
    return false;
}

inline STRING_VECTOR
GetPot::unidentified_arguments(unsigned Number,
			       const char* KnownArgument1, ...) const
{
    std::set<std::string> known_arguments;

    // (1) create a vector of known arguments
    if( Number == 0 ) return STRING_VECTOR();

    va_list ap;
    va_start(ap, KnownArgument1);
    known_arguments.insert(std::string(KnownArgument1));
    unsigned i=1;
    for(; i<Number; i++)
	known_arguments.insert(std::string(va_arg(ap, char *)));
    va_end(ap);

    return unidentified_arguments(known_arguments);
}

inline STRING_VECTOR
GetPot::unidentified_arguments() const
{ return unidentified_arguments(_requested_arguments); }

inline STRING_VECTOR
GetPot::unidentified_arguments(const std::vector<std::string>& Knowns) const
{
    // We use set for efficiency, but want to support vector inputs for
    // backwards compatibility.
    return unidentified_arguments(std::set<std::string> (Knowns.begin(), Knowns.end()));
}

inline STRING_VECTOR
GetPot::unidentified_arguments(const std::set<std::string>& Knowns) const
{
    STRING_VECTOR ufos;
    STRING_VECTOR::const_iterator it = argv.begin();
    ++it; // forget about argv[0] (application or filename)
    for(; it != argv.end(); ++it) {
	// -- argument belongs to prefixed section ?
	const std::string arg = _get_remaining_string(*it, prefix);
	if( arg == "" ) continue;

	// -- check if in list
	if( Knowns.find(arg) == Knowns.end() )
	    ufos.push_back(*it);
    }
    return ufos;
}

inline STRING_VECTOR
GetPot::unidentified_options(unsigned Number,
			     const char* KnownOption1, ...) const
{
    std::set<std::string> known_options;

    // (1) create a vector of known arguments
    if( Number == 0 ) return STRING_VECTOR();

    va_list ap;
    va_start(ap, KnownOption1);
    known_options.insert(std::string(KnownOption1));
    unsigned i=1;
    for(; i<Number; i++)
	known_options.insert(std::string(va_arg(ap, char *)));
    va_end(ap);

    return unidentified_options(known_options);
}

inline STRING_VECTOR
GetPot::unidentified_options() const
{
    // -- every option is an argument.
    // -- the set of requested arguments contains the set of requested options.
    // -- IF the set of requested arguments contains unrequested options,
    //    THEN they were requested as 'follow' and 'next' arguments and not as real options.
    //
    // => it is not necessary to separate requested options from the list
    return unidentified_arguments(_requested_arguments);
}

inline STRING_VECTOR
GetPot::unidentified_options(const std::vector<std::string>& Knowns) const
{
    // We use set for efficiency, but want to support vector inputs for
    // backwards compatibility.
    return unidentified_options(std::set<std::string> (Knowns.begin(), Knowns.end()));
}

inline STRING_VECTOR
GetPot::unidentified_options(const std::set<std::string>& Knowns) const
{
    STRING_VECTOR ufos;
    STRING_VECTOR::const_iterator it = argv.begin();
    ++it; // forget about argv[0] (application or filename)
    for(; it != argv.end(); ++it) {
	// -- argument belongs to prefixed section ?
	const std::string arg = _get_remaining_string(*it, prefix);
	if( arg == "" ) continue;

	// is argument really an option (starting with '-') ?
	if( arg.length() < 1 || arg[0] != '-' ) continue;

	if( Knowns.find(arg) == Knowns.end() )
	    ufos.push_back(*it);
    }

    return ufos;
}

inline std::string
GetPot::unidentified_flags(const char* KnownFlagList, int ArgumentNumber=-1) const
    // Two modes:
    //  ArgumentNumber >= 0 check specific argument
    //  ArgumentNumber == -1 check all options starting with one '-'
    //                       for flags
{
    std::string         ufos;
    STRING_VECTOR known_arguments;
    std::string         KFL(KnownFlagList);

    // (2) iteration over '-' arguments (options)
    if( ArgumentNumber == -1 ) {
	STRING_VECTOR::const_iterator it = argv.begin();
	it++; // forget about argv[0] (application or filename)
	for(; it != argv.end(); ++it) {
	    // -- argument belongs to prefixed section ?
	    const std::string arg = _get_remaining_string(*it, prefix);
	    if( arg == "" ) continue;

	    // -- does arguments start with '-' (but not '--')
	    if     ( arg.length() < 2 ) continue;
	    else if( arg[0] != '-' )    continue;
	    else if( arg[1] == '-' )    continue;

	    // -- check out if flags inside option are contained in KnownFlagList
	    const char* p=arg.c_str();
	    p++; // skip starting minus
	    for(; *p != '\0' ; p++)
		if( KFL.find(*p) == std::string::npos ) ufos += *p;
	}
    }
    // (1) check specific argument
    else {
	// -- only check arguments that start with prefix
	int no_matches = 0;
	unsigned i=1;
	for(; i<argv.size(); i++) {
	    const std::string Remain = _get_remaining_string(argv[i], prefix);
	    if( Remain != "") {
		no_matches++;
		if( no_matches == ArgumentNumber) {
		    // -- the right argument number inside the section is found
		    // => check it for flags
		    const char* p = Remain.c_str();
		    p++; // skip starting minus
		    for(; *p != '\0' ; p++)
			if( KFL.find(*p) == std::string::npos ) ufos += *p;
		    return ufos;
		}
	    }
	}
    }
    return ufos;
}

inline STRING_VECTOR
GetPot::unidentified_variables(unsigned Number,
			       const char* KnownVariable1, ...) const
{
    std::set<std::string> known_variables;

    // create vector of known arguments
    if( Number == 0 ) return STRING_VECTOR();

    va_list ap;
    va_start(ap, KnownVariable1);
    known_variables.insert(std::string(KnownVariable1));
    unsigned i=1;
    for(; i<Number; i++)
	known_variables.insert(std::string(va_arg(ap, char *)));
    va_end(ap);

    return unidentified_variables(known_variables);
}

inline STRING_VECTOR
GetPot::unidentified_variables(const std::vector<std::string>& Knowns) const
{
    // We use set for efficiency, but want to support vector inputs for
    // backwards compatibility.
    return unidentified_variables(std::set<std::string> (Knowns.begin(), Knowns.end()));
}

inline STRING_VECTOR
GetPot::unidentified_variables(const std::set<std::string>& Knowns) const
{
    STRING_VECTOR ufos;

    victorate(GetPot::variable, variables, it) {
	// -- check if variable has specific prefix
	const std::string var_name = _get_remaining_string((*it).name, prefix);
	if( var_name == "" ) continue;

	// -- check if variable is known
	if( Knowns.find(var_name) == Knowns.end() )
	    ufos.push_back((*it).name);
    }
    return ufos;
}

inline STRING_VECTOR
GetPot::unidentified_variables() const
{  return unidentified_variables(_requested_variables); }


inline STRING_VECTOR
GetPot::unidentified_sections(unsigned Number,
			      const char* KnownSection1, ...) const
{
    std::set<std::string> known_sections;

    // (1) create a vector of known arguments
    if( Number == 0 ) return STRING_VECTOR();

    va_list ap;
    va_start(ap, KnownSection1);
    known_sections.insert(std::string(KnownSection1));
    unsigned i=1;
    for(; i<Number; i++) {
	std::string tmp = std::string(va_arg(ap, char *));
	if( tmp.length() == 0 ) continue;
	if( tmp[tmp.length()-1] != '/' ) tmp += '/';
	known_sections.insert(tmp);
    }
    va_end(ap);

    return unidentified_sections(known_sections);
}

inline STRING_VECTOR
GetPot::unidentified_sections() const
{ return unidentified_sections(_requested_sections); }

inline STRING_VECTOR
GetPot::unidentified_sections(const std::vector<std::string>& Knowns) const
{
    // We use set for efficiency, but want to support vector inputs for
    // backwards compatibility.
    return unidentified_sections(std::set<std::string> (Knowns.begin(), Knowns.end()));
}

inline STRING_VECTOR
GetPot::unidentified_sections(const std::set<std::string>& Knowns) const
{
    STRING_VECTOR ufos;

    victorate(std::string, section_list, it) {
	// -- check if section conform to prefix
	const std::string sec_name = _get_remaining_string(*it, prefix);
	if( sec_name == "" ) continue;

	// -- check if section is known
	if( Knowns.find(sec_name) == Knowns.end() )
	    ufos.push_back(*it);
    }

    return ufos;
}


inline STRING_VECTOR
GetPot::unidentified_nominuses(unsigned Number, const char* Known, ...) const
{
    std::set<std::string> known_nominuses;

    // create vector of known arguments
    if( Number == 0 ) return STRING_VECTOR();

    va_list ap;
    va_start(ap, Known);
    known_nominuses.insert(std::string(Known));
    unsigned i=1;
    for(; i<Number; i++) {
	std::string tmp = std::string(va_arg(ap, char *));
	if( tmp.length() == 0 ) continue;
	known_nominuses.insert(tmp);
    }
    va_end(ap);

    return unidentified_nominuses(known_nominuses);
}

inline STRING_VECTOR
GetPot::unidentified_nominuses() const {
    // -- every nominus is an argument.
    // -- the set of requested arguments contains the set of requested nominuss.
    // -- IF the set of requested arguments contains unrequested nominuss,
    //    THEN they were requested as 'follow' and 'next' arguments and not as real nominuses.
    //
    // => it is not necessary to separate requested nominus from the list

    return unidentified_nominuses(_requested_arguments);
}

inline STRING_VECTOR
GetPot::unidentified_nominuses(const std::vector<std::string>& Knowns) const
{
    // We use set for efficiency, but want to support vector inputs for
    // backwards compatibility.
    return unidentified_nominuses(std::set<std::string> (Knowns.begin(), Knowns.end()));
}

inline STRING_VECTOR
GetPot::unidentified_nominuses(const std::set<std::string>& Knowns) const
{
    STRING_VECTOR ufos;

    // (2) iterate over all arguments
    STRING_VECTOR::const_iterator it = argv.begin();
    ++it; // forget about argv[0] (application or filename)
    for(; it != argv.end(); ++it) {
	// -- check if nominus part of prefix
	const std::string arg = _get_remaining_string(*it, prefix);
	if( arg == "" )                                         continue;

	if( arg.length() < 1 )                                  continue;
	// option ? --> not a nomius
	if( arg[0] == '-' )                                     continue;
	// section ? --> not a real nominus
	if( arg[0] == '[' && arg[arg.length()-1] == ']' )       continue;
	// variable definition ? --> not a real nominus
	bool continue_f = false;
	unsigned i=0;
	for(; i<arg.length() ; i++)
	    if( arg[i] == '=' ) { continue_f = true; break; }
	if( continue_f )                                        continue;

	// real nominuses are compared with the given list
	if( Knowns.find(arg) == Knowns.end() )
	    ufos.push_back(*it);
    }
    return ufos;
}


///////////////////////////////////////////////////////////////////////////////
// (*) variable class
//.............................................................................
//
inline
GetPot::variable::variable()
 : name(),
   value(),
   original()
{}

inline
GetPot::variable::variable(const variable& Other)
{
#ifdef WIN32
    operator=(Other);
#else
    GetPot::variable::operator=(Other);
#endif
}


inline
GetPot::variable::variable(const char* Name, const char* Value, const char* FieldSeparator)
    : name(Name)
{
    // make a copy of the 'Value'
    take(Value, FieldSeparator);
}

inline const std::string*
GetPot::variable::get_element(unsigned Idx) const
{ if( Idx >= value.size() ) return 0; else return &(value[Idx]); }

inline void
GetPot::variable::take(const char* Value, const char* FieldSeparator)
{
  original = std::string(Value); // string member var
  value.clear();                 // vector<string> member var

  /*
  // separate string by white space delimiters using 'strtok'
  // thread safe usage of strtok (no static members)
  char* spt = 0;
  // make a copy of the 'Value'
  char* copy = new char[strlen(Value)+1];
  strcpy(copy, Value);
  char* follow_token = strtok_r(copy, FieldSeparator, &spt);
  while(follow_token != 0) {
    value.push_back(std::string(follow_token));
    follow_token = strtok_r(NULL, FieldSeparator, &spt);
  }

  delete [] copy;
  */

  // Don't use strtok, instead tokenize the input char "Value" using std::string operations so
  // that the results end up in the local "value" member

  // Construct std::string objects from the input char*s.  I think the only
  // FieldSeparator recognized by GetPot is whitespace?
  std::string Value_str = std::string(Value);
  std::string delimiters = std::string(FieldSeparator);

  // Skip delimiters at beginning.
  std::string::size_type lastPos = Value_str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  std::string::size_type pos     = Value_str.find_first_of(delimiters, lastPos);

  // Loop over the input string until all the tokens have been pushed back
  // into the local "value" member.
  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      value.push_back(Value_str.substr(lastPos, pos - lastPos));

      // Skip delimiters.  Note the "not_of"
      lastPos = Value_str.find_first_not_of(delimiters, pos);

      // Find next "non-delimiter"
      pos = Value_str.find_first_of(delimiters, lastPos);
    }

  // We're done, all the tokens should now be in the vector<string>

}

inline
GetPot::variable::~variable()
{}

inline GetPot::variable&
GetPot::variable::operator=(const GetPot::variable& Other)
{
    if( &Other != this) {
	name     = Other.name;
	value    = Other.value;
	original = Other.original;
    }
    return *this;
}

#ifdef GETPOT_NAMESPACE
}
#endif

#undef victorate

#endif // LIBMESH_GETPOT_H
