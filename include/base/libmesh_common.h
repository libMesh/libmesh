// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef LIBMESH_LIBMESH_COMMON_H
#define LIBMESH_LIBMESH_COMMON_H

// These flags should never be used together.  -DDEBUG means "turn on
// all the ridiculously expensive error checking"; -DNDEBUG means
// "turn off all the somewhat affordable error checking"
#if defined(DEBUG) && defined(NDEBUG)
#  error DEBUG and NDEBUG should never be defined simultaneously
#endif

// The library configuration options
#include "libmesh/libmesh_config.h"

// Use actual timestamps or constant dummies (to aid ccache)
#ifdef LIBMESH_ENABLE_TIMESTAMPS
#  define  LIBMESH_TIME __TIME__
#  define  LIBMESH_DATE __DATE__
#else
#  define  LIBMESH_TIME "notime"
#  define  LIBMESH_DATE "nodate"
#endif

// C/C++ includes everyone should know about
#include <cstdlib>
#ifdef __PGI
// BSK, Thu Feb 20 08:32:06 CST 2014 - For some reason, unless PGI gets
// <cmath> early this nonsense shows up:
// "/software/x86_64/pgi/12.9/linux86-64/12.9/include/CC/cmath", line 57: error:
//           the global scope has no "abs"
//   using _STLP_VENDOR_CSTD::abs;
// So include <cmath> as early as possible under the PGI compilers.
#  include <cmath>
#endif
#include <complex>
#include <typeinfo> // std::bad_cast
#include <type_traits> // std::decay
#include <functional> // std::less, etc
#include <string>

// Include the MPI definition
#ifdef LIBMESH_HAVE_MPI
# include "libmesh/ignore_warnings.h"
# include <mpi.h>
# include "libmesh/restore_warnings.h"
#endif

// Quad precision if we need it
#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
#include "libmesh/float128_shims.h"
#endif

// _basic_ library functionality
#include "libmesh/libmesh_base.h"
#include "libmesh/libmesh_exceptions.h"

// Proxy class for libMesh::out/err output
#include "libmesh/ostream_proxy.h"

// Make sure the libmesh_nullptr define is available for backwards
// compatibility, although we no longer use it in the library.
#include "libmesh/libmesh_nullptr.h"

// C++ headers
#include <iomanip> // setprecision, in assertion macros

namespace libMesh
{
namespace Threads
{
// For thread-safe error-messaging. Definitions in threads.h
void lock_singleton_spin_mutex();
void unlock_singleton_spin_mutex();
}

// Let's define a couple output streams - these will default
// to cout/cerr, but LibMeshInit (or the user) can also set them to
// something more sophisticated.
//
// We use a proxy class rather than references so they can be
// reseated at runtime.

extern OStreamProxy out;
extern OStreamProxy err;

// A namespace for functions used in the bodies of the macros below.
// The macros generally call these functions with __FILE__, __LINE__,
// __DATE__, and __TIME__ in the appropriate order.  These should not
// be called by users directly!  The implementations can be found in
// libmesh_common.C.
namespace MacroFunctions
{
void here(const char * file, int line, const char * date, const char * time, std::ostream & os = libMesh::err);
void stop(const char * file, int line, const char * date, const char * time);
void report_error(const char * file, int line, const char * date, const char * time, std::ostream & os = libMesh::err);
}

// Undefine any existing macros
#ifdef Real
#  undef Real
#endif

//#ifdef REAL
//#  undef REAL
//#endif

#ifdef Complex
#  undef Complex
#endif

#ifdef COMPLEX
#  undef COMPLEX
#endif

// Check to see if TOLERANCE has been defined by another
// package, if so we might want to change the name...
#ifdef TOLERANCE
DIE A HORRIBLE DEATH HERE...
#  undef TOLERANCE
#endif



// Define the type to use for real numbers

typedef LIBMESH_DEFAULT_SCALAR_TYPE Real;

// Define a corresponding tolerance.  This is what should be
// considered "good enough" when doing floating point comparisons.
// For example, v == 0 is changed to std::abs(v) < TOLERANCE.

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
static constexpr Real TOLERANCE = 2.5e-3;
# if defined (LIBMESH_DEFAULT_TRIPLE_PRECISION) || \
     defined (LIBMESH_DEFAULT_QUADRUPLE_PRECISION)
#  error Cannot define multiple precision levels
# endif
#endif

#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
static constexpr Real TOLERANCE = 1.e-8;
# if defined (LIBMESH_DEFAULT_QUADRUPLE_PRECISION)
#  error Cannot define multiple precision levels
# endif
#endif

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
static constexpr Real TOLERANCE = 1.e-11;
#endif

#if !defined (LIBMESH_DEFAULT_SINGLE_PRECISION) && \
    !defined (LIBMESH_DEFAULT_TRIPLE_PRECISION) && \
    !defined (LIBMESH_DEFAULT_QUADRUPLE_PRECISION)
static constexpr Real TOLERANCE = 1.e-6;
#endif

// Define the type to use for complex numbers
// Always use std::complex<double>, as required by Petsc?
// If your version of Petsc doesn't support
// std::complex<other_precision>, then you'd better just leave
// Real==double
typedef std::complex<Real> Complex;
typedef std::complex<Real> COMPLEX;


// Helper functions for complex/real numbers
// to clean up #ifdef LIBMESH_USE_COMPLEX_NUMBERS elsewhere
template<typename T> inline T libmesh_real(T a) { return a; }
template<typename T> inline T libmesh_imag(T /*a*/) { return 0; }
template<typename T> inline T libmesh_conj(T a) { return a; }

template<typename T>
inline T libmesh_real(std::complex<T> a) { return std::real(a); }

template<typename T>
inline T libmesh_imag(std::complex<T> a) { return std::imag(a); }

template<typename T>
inline std::complex<T> libmesh_conj(std::complex<T> a) { return std::conj(a); }

// std::isnan() is in <cmath> as of C++11.
template <typename T>
inline bool libmesh_isnan(T x) { return std::isnan(x); }

template <typename T>
inline bool libmesh_isnan(std::complex<T> a)
{ return (std::isnan(std::real(a)) || std::isnan(std::imag(a))); }

// std::isinf() is in <cmath> as of C++11.
template <typename T>
inline bool libmesh_isinf(T x) { return std::isinf(x); }

template <typename T>
inline bool libmesh_isinf(std::complex<T> a)
{ return (std::isinf(std::real(a)) || std::isinf(std::imag(a))); }

// Define the value type for unknowns in simulations.
// This is either Real or Complex, depending on how
// the library was configures
#if   defined (LIBMESH_USE_REAL_NUMBERS)
typedef Real Number;
#elif defined (LIBMESH_USE_COMPLEX_NUMBERS)
typedef Complex Number;
#else
DIE A HORRIBLE DEATH HERE...
#endif


// Define the value type for error estimates.
// Since AMR/C decisions don't have to be precise,
// we default to float for memory efficiency.
typedef float ErrorVectorReal;
#define MPI_ERRORVECTORREAL MPI_FLOAT


#ifdef LIBMESH_HAVE_MPI

/**
 * MPI Communicator used to initialize libMesh.
 */
extern MPI_Comm GLOBAL_COMM_WORLD;
#else

/**
 * Something to use with CHKERRABORT if we're just using PETSc's MPI
 * "uni" stub.
 */
extern int GLOBAL_COMM_WORLD;
#endif

// This global variable is to help us deprecate AutoPtr.  We can't
// just use libmesh_deprecated() because then you get one print out
// per template instantiation, instead of one total print out.
extern bool warned_about_auto_ptr;

// These are useful macros that behave like functions in the code.
// If you want to make sure you are accessing a section of code just
// stick a libmesh_here(); in it, for example
#define libmesh_here()                                                  \
  do {                                                                  \
    libMesh::MacroFunctions::here(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME); \
  } while (0)

// the libmesh_stop() macro will stop the code until a SIGCONT signal
// is received.  This is useful, for example, when determining the
// memory used by a given operation.  A libmesh_stop() could be
// inserted before and after a questionable operation and the delta
// memory can be obtained from a ps or top.  This macro only works for
// serial cases.
#define libmesh_stop()                                                  \
  do {                                                                  \
    libMesh::MacroFunctions::stop(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME); \
  } while (0)

// The libmesh_dbg_var() macro indicates that an argument to a function
// is used only in debug and devel modes (i.e., when NDEBUG is not defined).
#ifndef NDEBUG
#define libmesh_dbg_var(var) var
#else
#define libmesh_dbg_var(var)
#endif

// The libmesh_inf_var() macro indicates that an argument to a function
// is used only when infinite elements are enabled
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#define libmesh_inf_var(var) var
#else
#define libmesh_inf_var(var)
#endif

// The libmesh_assert() macro acts like C's assert(), but throws a
// libmesh_error() (including stack trace, etc) instead of just exiting
#ifdef NDEBUG

#define libmesh_assert_msg(asserted, msg)  ((void) 0)
#define libmesh_exceptionless_assert_msg(asserted, msg)  ((void) 0)
#define libmesh_assert_equal_to_msg(expr1,expr2, msg)  ((void) 0)
#define libmesh_assert_not_equal_to_msg(expr1,expr2, msg)  ((void) 0)
#define libmesh_assert_less_msg(expr1,expr2, msg)  ((void) 0)
#define libmesh_assert_greater_msg(expr1,expr2, msg)  ((void) 0)
#define libmesh_assert_less_equal_msg(expr1,expr2, msg)  ((void) 0)
#define libmesh_assert_greater_equal_msg(expr1,expr2, msg)  ((void) 0)

#else

#define libmesh_assertion_types(expr1,expr2)                            \
  typedef typename std::decay<decltype(expr1)>::type libmesh_type1;     \
  typedef typename std::decay<decltype(expr2)>::type libmesh_type2

#define libmesh_assert_msg(asserted, msg)                               \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      libmesh_error_msg(msg);                                           \
    } } while (0)

#define libmesh_exceptionless_assert_msg(asserted, msg)                 \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      libMesh::Threads::lock_singleton_spin_mutex();                    \
      libMesh::err << "Assertion `" #asserted "' failed." << std::endl; \
      libMesh::Threads::unlock_singleton_spin_mutex();                  \
      libmesh_exceptionless_error();                                    \
    } } while (0)

#define libmesh_assert_equal_to_msg(expr1,expr2, msg)                   \
  do {                                                                  \
    if (!((expr1) == (expr2))) {                                        \
      libmesh_error_msg(std::setprecision(17) << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl); \
    } } while (0)

#define libmesh_assert_not_equal_to_msg(expr1,expr2, msg)               \
  do {                                                                  \
    if (!((expr1) != (expr2))) {                                        \
      libmesh_error_msg(std::setprecision(17) << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl); \
    } } while (0)

template <template <class> class Comp>
struct casting_compare {

  template <typename T1, typename T2>
  bool operator()(const T1 & e1, const T2 & e2) const
  {
    typedef typename std::decay<T1>::type DT1;
    typedef typename std::decay<T2>::type DT2;
    return (Comp<DT2>()(static_cast<DT2>(e1), e2) &&
            Comp<DT1>()(e1, static_cast<DT1>(e2)));
  }

  template <typename T1>
  bool operator()(const T1 & e1, const T1 & e2) const
  {
    return Comp<T1>()(e1, e2);
  }
};

#define libmesh_assert_less_msg(expr1,expr2, msg)                       \
  do {                                                                  \
    if (!libMesh::casting_compare<std::less>()(expr1, expr2)) {         \
      libmesh_error_msg(std::setprecision(17) << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl); \
    } } while (0)

#define libmesh_assert_greater_msg(expr1,expr2, msg)                    \
  do {                                                                  \
    if (!libMesh::casting_compare<std::greater>()(expr1, expr2)) {      \
      libmesh_error_msg(std::setprecision(17) << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl); \
    } } while (0)

#define libmesh_assert_less_equal_msg(expr1,expr2, msg)                 \
  do {                                                                  \
    if (!libMesh::casting_compare<std::less_equal>()(expr1, expr2)) {   \
      libmesh_error_msg(std::setprecision(17) << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl); \
    } } while (0)

#define libmesh_assert_greater_equal_msg(expr1,expr2, msg)              \
  do {                                                                  \
    if (!libMesh::casting_compare<std::greater_equal>()(expr1, expr2)) { \
      libmesh_error_msg(std::setprecision(17) << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl); \
    } } while (0)

#endif


#define libmesh_assert(asserted) libmesh_assert_msg(asserted, "")
#define libmesh_exceptionless_assert(asserted) libmesh_exceptionless_assert_msg(asserted, "")
#define libmesh_assert_equal_to(expr1,expr2) libmesh_assert_equal_to_msg(expr1,expr2, "")
#define libmesh_assert_not_equal_to(expr1,expr2) libmesh_assert_not_equal_to_msg(expr1,expr2, "")
#define libmesh_assert_less(expr1,expr2) libmesh_assert_less_msg(expr1,expr2, "")
#define libmesh_assert_greater(expr1,expr2) libmesh_assert_greater_msg(expr1,expr2, "")
#define libmesh_assert_less_equal(expr1,expr2) libmesh_assert_less_equal_msg(expr1,expr2, "")
#define libmesh_assert_greater_equal(expr1,expr2) libmesh_assert_greater_equal_msg(expr1,expr2, "")

// The libmesh_error() macro prints a message and throws a LogicError
// exception
//
// The libmesh_not_implemented() macro prints a message and throws a
// NotImplemented exception
//
// The libmesh_file_error(const std::string & filename) macro prints a message
// and throws a FileError exception
//
// The libmesh_convergence_failure() macro
// throws a ConvergenceFailure exception
//
// The libmesh_degenerate_mapping() macro prints a message into and
// throws a DegenerateMap exception
//
// The libmesh_terminate() macro prints a message and throws a
// TerminationException exception
#define libmesh_error_msg(msg)                                          \
  do {                                                                  \
    std::stringstream message_stream;                                   \
    message_stream << msg << '\n';                                      \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME, message_stream); \
    LIBMESH_THROW(libMesh::LogicError(message_stream.str()));           \
  } while (0)

#define libmesh_error() libmesh_error_msg("")

#define libmesh_error_msg_if(cond, msg)         \
  do {                                          \
    if (cond)                                   \
      libmesh_error_msg(msg);                   \
  } while (0)

#define libmesh_exceptionless_error_msg(msg)                            \
  do {                                                                  \
    libMesh::Threads::lock_singleton_spin_mutex();                      \
    libMesh::err << msg << '\n';                                        \
    libMesh::Threads::unlock_singleton_spin_mutex();                    \
    libmesh_try { libMesh::MacroFunctions::report_error(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME); } \
    libmesh_catch (...) {}                                              \
    std::terminate();                                                   \
  } while (0)

#define libmesh_exceptionless_error() libmesh_exceptionless_error_msg("")

#define libmesh_not_implemented_msg(msg)                                \
  do {                                                                  \
    std::stringstream message_stream;                                   \
    message_stream << msg << '\n';                                      \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME, message_stream); \
    LIBMESH_THROW(libMesh::NotImplemented(message_stream.str()));       \
  } while (0)

#define libmesh_not_implemented() libmesh_not_implemented_msg("")

#define libmesh_file_error_msg(filename, msg)                           \
  do {                                                                  \
    std::stringstream message_stream;                                   \
    message_stream << msg << '\n';                                      \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, LIBMESH_DATE, LIBMESH_TIME, message_stream); \
    LIBMESH_THROW(libMesh::FileError(filename, message_stream.str()));  \
  } while (0)

#define libmesh_file_error(filename) libmesh_file_error_msg(filename,"")

#define libmesh_convergence_failure()                   \
  do {                                                  \
    LIBMESH_THROW(libMesh::ConvergenceFailure());       \
  } while (0)

#define libmesh_degenerate_mapping_msg(msg)                          \
  do {                                                               \
    std::stringstream message_stream;                                \
    message_stream << msg << '\n';                                   \
    LIBMESH_THROW(libMesh::DegenerateMap(message_stream.str())); \
  } while (0)

#define libmesh_degenerate_mapping(filename) libmesh_degenerate_mapping_msg("")

#define libmesh_terminate()                         \
  do {                                              \
    LIBMESH_THROW(libMesh::TerminationException()); \
  } while (0)

// The libmesh_example_requires() macro prints a message and calls
// "return 77;" if the condition specified by the macro is not true.  This
// macro is used in the example executables, which should run when the
// configure-time libMesh options support them but which should exit
// without failure otherwise.
//
// This macro only works in main(), because we have no better way than
// "return" from main to immediately exit successfully - std::exit()
// gets seen by at least some MPI stacks as failure.
//
// 77 is the automake code for a skipped test.

#define libmesh_example_requires(condition, option)                     \
  do {                                                                  \
    if (!(condition)) {                                                 \
      libMesh::out << "Configuring libMesh with " << option << " is required to run this example." << std::endl; \
      return 77;                                                        \
    } } while (0)

// The libmesh_do_once macro helps us avoid redundant repeated
// repetitions of the same warning messages
#undef libmesh_do_once
#define libmesh_do_once(do_this)                \
  do {                                          \
    static bool did_this_already = false;       \
    if (!did_this_already) {                    \
      did_this_already = true;                  \
      do_this;                                  \
    } } while (0)


// The libmesh_warning macro outputs a file/line/time stamped warning
// message, if warnings are enabled.
#ifdef LIBMESH_ENABLE_WARNINGS
#define libmesh_warning(message)                                        \
  libmesh_do_once(libMesh::out << message                               \
                  << __FILE__ << ", line " << __LINE__ << ", compiled " << LIBMESH_DATE << " at " << LIBMESH_TIME << " ***" << std::endl;)
#else
#define libmesh_warning(message)  ((void) 0)
#endif

// The libmesh_experimental macro warns that you are using
// bleeding-edge code
#undef libmesh_experimental
#define libmesh_experimental()                                          \
  libmesh_warning("*** Warning, This code is untested, experimental, or likely to see future API changes: ");


// The libmesh_deprecated macro warns that you are using obsoleted code
#undef libmesh_deprecated
#ifndef LIBMESH_ENABLE_DEPRECATED
#define libmesh_deprecated()                                            \
  libmesh_error_msg("*** Error, This code is deprecated, and likely to be removed in future library versions! ");
#else
#define libmesh_deprecated()                                            \
  libmesh_warning("*** Warning, This code is deprecated, and likely to be removed in future library versions! ");
#endif

// A function template for ignoring unused variables.  This is a way
// to shut up unused variable compiler warnings on a case by case
// basis.
template<class ...Args> inline void libmesh_ignore( const Args&... ) { }


// A workaround for the lack of C++17 merge() support in some
// compilers

#ifdef LIBMESH_HAVE_CXX17_SPLICING
template <typename T>
void libmesh_merge_move(T & target, T & source)
{
  target.merge(std::move(source));
}
#else
template <typename T>
void libmesh_merge_move(T & target, T & source)
{
  target.insert(source.begin(), source.end());
  source.clear(); // Avoid forwards-incompatibility
}
#endif // LIBMESH_HAVE_CXX17_SPLICING

/**
 * Mostly system independent demangler
 */
 std::string demangle(const char * name);

// cast_ref and cast_ptr do a dynamic cast and assert
// the result, if we have RTTI enabled and we're in debug or
// development modes, but they just do a faster static cast if we're
// in optimized mode.
//
// Use these casts when you're certain that a cast will succeed in
// correct code but you want to be able to double-check.
template <typename Tnew, typename Told>
inline Tnew cast_ref(Told & oldvar)
{
#if !defined(NDEBUG) && defined(LIBMESH_HAVE_RTTI) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  try
    {
      Tnew newvar = dynamic_cast<Tnew>(oldvar);
      return newvar;
    }
  catch (std::bad_cast &)
    {
      libMesh::err << "Failed to convert " << demangle(typeid(Told).name())
                   << " reference to " << demangle(typeid(Tnew).name())
                   << std::endl;
      libMesh::err << "The " << demangle(typeid(Told).name())
                   << " appears to be a "
                   << demangle(typeid(*(&oldvar)).name()) << std::endl;
      libmesh_error();
    }
#else
  return(static_cast<Tnew>(oldvar));
#endif
}

// We use two different function names to avoid an odd overloading
// ambiguity bug with icc 10.1.008
template <typename Tnew, typename Told>
inline Tnew cast_ptr (Told * oldvar)
{
#if !defined(NDEBUG) && defined(LIBMESH_HAVE_RTTI)
  Tnew newvar = dynamic_cast<Tnew>(oldvar);
  if (!newvar)
    {
      libMesh::err << "Failed to convert " << demangle(typeid(Told).name())
                   << " pointer to " << demangle(typeid(Tnew).name())
                   << std::endl;
      libMesh::err << "The " << demangle(typeid(Told).name())
                   << " appears to be a "
                   << demangle(typeid(*oldvar).name()) << std::endl;
      libmesh_error();
    }
  return newvar;
#else
  return(static_cast<Tnew>(oldvar));
#endif
}


#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename Tnew, typename Told>
inline Tnew libmesh_cast_ptr (Told * oldvar)
{
  libmesh_deprecated();

  // we use the less redundantly named libMesh::cast_ptr now
  return cast_ptr<Tnew>(oldvar);
}
#endif // LIBMESH_ENABLE_DEPRECATED


// cast_int asserts that the value of the castee is within the
// bounds which are exactly representable by the output type, if we're
// in debug or development modes, but it just does a faster static
// cast if we're in optimized mode.
//
// Use these casts when you're certain that a cast will succeed in
// correct code but you want to be able to double-check.
template <typename Tnew, typename Told>
inline Tnew cast_int (Told oldvar)
{
  libmesh_assert_equal_to
    (oldvar, static_cast<Told>(static_cast<Tnew>(oldvar)));

  return(static_cast<Tnew>(oldvar));
}


template <typename Tnew, typename Told>
inline Tnew libmesh_cast_int (Told oldvar)
{
  // we use the less redundantly named libMesh::cast_int now
  return cast_int<Tnew>(oldvar);
}

/**
 * This is a helper variable template for cases when we want to use a default compile-time
 * error with constexpr-based if conditions. The templating delays the triggering
 * of the static assertion until the template is instantiated.
 */
template <class T>
constexpr std::false_type always_false{};

static constexpr std::size_t libmesh_dim = LIBMESH_DIM;

// build a integer representation of version
#define LIBMESH_VERSION_ID(major,minor,patch) (((major) << 16) | ((minor) << 8) | ((patch) & 0xFF))


// libmesh_override is simply a synonym for override as we now require
// a C++11 compiler that supports this keyword.
#define libmesh_override override

// libmesh_delete is simply a synonym for '=delete' as we now require
// a C++11 compiler that supports this keyword.
#define libmesh_delete =delete

// libmesh_final is simply a synonym for 'final' as we now require
// a C++11 compiler that supports this keyword.
#define libmesh_final final

// Define backwards-compatible fallthrough attribute.  We could
// eventually also add support for other compiler-specific fallthrough
// attributes.
#ifdef LIBMESH_HAVE_CXX17_FALLTHROUGH_ATTRIBUTE
#define libmesh_fallthrough() [[fallthrough]]
#elif defined(LIBMESH_HAVE_DOUBLE_UNDERSCORE_ATTRIBUTE_FALLTHROUGH)
#define libmesh_fallthrough() __attribute__((fallthrough))
#else
#define libmesh_fallthrough() ((void) 0)
#endif

template <typename T>
class PassKey
{
  friend T;
  constexpr PassKey() = default;
};
} // namespace libMesh


// Backwards compatibility
namespace libMeshEnums
{
using namespace libMesh;
}

// Backwards compatibility with pre-TIMPI reference
namespace TIMPI {}

namespace libMesh {
  namespace Parallel {
    using namespace TIMPI;
  }
}


// Here we add missing types to the standard namespace.  For example,
// std::max(double, float) etc... are well behaved but not defined
// by the standard.  This also includes workarounds for super-strict
// implementations, for example Sun Studio and PGI C++.  However,
// this necessarily requires breaking the ISO-C++ standard, and is
// really just a hack.  As such, only do it if we are building the
// libmesh library itself.  Specifically, *DO NOT* export this to
// user code or install this header.
//
// We put this at the end of libmesh_common.h so we can make use of
// any exotic definitions of Real above.
#ifdef LIBMESH_IS_COMPILING_ITSELF
#  include "libmesh/libmesh_augment_std_namespace.h"
#endif


#ifdef _MSC_VER
#define LIBMESH_EXPORT __declspec(dllexport)
#else
#define LIBMESH_EXPORT
#endif


#endif // LIBMESH_LIBMESH_COMMON_H
