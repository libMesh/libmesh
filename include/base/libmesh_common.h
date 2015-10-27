// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// The library configuration options
#include "libmesh/libmesh_config.h"

// Use actual timestamps or constant dummies (to aid ccache)
#ifdef LIBMESH_ENABLE_TIMESTAMPS
#  define  __LIBMESH_TIME__ __TIME__
#  define  __LIBMESH_DATE__ __DATE__
#else
#  define  __LIBMESH_TIME__ "notime"
#  define  __LIBMESH_DATE__ "nodate"
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


// Include the MPI definition
#ifdef LIBMESH_HAVE_MPI
# include "libmesh/ignore_warnings.h"
# include <mpi.h>
# include "libmesh/restore_warnings.h"
#endif

// _basic_ library functionality
#include "libmesh/libmesh_base.h"
#include "libmesh/libmesh_exceptions.h"
extern "C" {
#include "libmesh/libmesh_C_isnan.h"
}

// Proxy class for libMesh::out/err output
#include "libmesh/ostream_proxy.h"

// Here we add missing types to the standard namespace.  For example,
// std::max(double, float) etc... are well behaved but not defined
// by the standard.  This also includes workarounds for super-strict
// implementations, for example Sun Studio and PGI C++.  However,
// this necessarily requires breaking the ISO-C++ standard, and is
// really just a hack.  As such, only do it if we are building the
// libmesh library itself.  Specifially, *DO NOT* export this to
// user code or install this header.
#ifdef LIBMESH_IS_COMPILING_ITSELF
#  include "libmesh/libmesh_augment_std_namespace.h"
#endif


namespace libMesh
{


// A namespace for functions used in the bodies of the macros below.
// The macros generally call these functions with __FILE__, __LINE__,
// __DATE__, and __TIME__ in the appropriate order.  These should not
// be called by users directly!  The implementations can be found in
// libmesh_common.C.
namespace MacroFunctions
{
void here(const char* file, int line, const char* date, const char* time);
void stop(const char* file, int line, const char* date, const char* time);
void report_error(const char* file, int line, const char* date, const char* time);
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

#ifdef MPI_REAL
#  undef MPI_REAL
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

#ifndef LIBMESH_DEFAULT_SINGLE_PRECISION
#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
static const Real TOLERANCE = 1.e-8;
# define MPI_REAL MPI_LONG_DOUBLE
#else
static const Real TOLERANCE = 1.e-6;
# define MPI_REAL MPI_DOUBLE
#endif
#else
static const Real TOLERANCE = 2.5e-3;
# define MPI_REAL MPI_FLOAT
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
template<typename T> inline T libmesh_conj(T a) { return a; }

template<typename T>
inline T libmesh_real(std::complex<T> a) { return std::real(a); }

template<typename T>
inline std::complex<T> libmesh_conj(std::complex<T> a) { return std::conj(a); }

// isnan isn't actually C++ standard yet; in contexts where it's not defined in
// cmath, libmesh_isnan will just end up returning false.
inline bool libmesh_isnan(float a) { return libmesh_C_isnan_float(a); }
inline bool libmesh_isnan(double a) { return libmesh_C_isnan_double(a); }
inline bool libmesh_isnan(long double a) { return libmesh_C_isnan_longdouble(a); }

template <typename T>
inline bool libmesh_isnan(std::complex<T> a) { return (libmesh_isnan(std::real(a)) || libmesh_isnan(std::imag(a))); }


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
#ifndef LIBMESH_DISABLE_COMMWORLD
/**
 * MPI Communicator to be used in the library.
 */
extern MPI_Comm COMM_WORLD;
#endif

/**
 * MPI Communicator used to initialize libMesh.
 */
extern MPI_Comm GLOBAL_COMM_WORLD;
#else
#ifndef LIBMESH_DISABLE_COMMWORLD
/**
 * Something to use with CHKERRABORT if we're just using PETSc's MPI
 * "uni" stub.
 */
extern int COMM_WORLD;
#endif

/**
 * Something to use with CHKERRABORT if we're just using PETSc's MPI
 * "uni" stub.
 */
extern int GLOBAL_COMM_WORLD;
#endif

// Let's define a couple output streams - these will default
// to cout/cerr, but LibMeshInit (or the user) can also set them to
// something more sophisticated.
//
// We use a proxy class rather than references so they can be
// reseated at runtime.

extern OStreamProxy out;
extern OStreamProxy err;

// This global variable is to help us deprecate AutoPtr.  We can't
// just use libmesh_deprecated() because then you get one print out
// per template instantiation, instead of one total print out.
extern bool warned_about_auto_ptr;

// These are useful macros that behave like functions in the code.
// If you want to make sure you are accessing a section of code just
// stick a libmesh_here(); in it, for example
#define libmesh_here()                                                  \
  do {                                                                  \
    libMesh::MacroFunctions::here(__FILE__, __LINE__, __LIBMESH_DATE__, __LIBMESH_TIME__); \
  } while(0)

// the libmesh_stop() macro will stop the code until a SIGCONT signal
// is recieved.  This is useful, for example, when determining the
// memory used by a given operation.  A libmesh_stop() could be
// instered before and after a questionable operation and the delta
// memory can be obtained from a ps or top.  This macro only works for
// serial cases.
#define libmesh_stop()                                                  \
  do {                                                                  \
    libMesh::MacroFunctions::stop(__FILE__, __LINE__, __LIBMESH_DATE__, __LIBMESH_TIME__); \
  } while(0)

// The libmesh_dbg_var() macro indicates that an argument to a function
// is used only in debug mode (i.e., when NDEBUG is not defined).
#ifndef NDEBUG
#define libmesh_dbg_var(var) var
#else
#define libmesh_dbg_var(var)
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

#define libmesh_assert_msg(asserted, msg)                               \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      libMesh::err << "Assertion `" #asserted "' failed." << std::endl; \
      libmesh_error_msg(msg);                                           \
    } } while(0)

#define libmesh_exceptionless_assert_msg(asserted, msg)                 \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      libMesh::err << "Assertion `" #asserted "' failed." << std::endl; \
      libmesh_exceptionless_error();                                    \
    } } while(0)

#define libmesh_assert_equal_to_msg(expr1,expr2, msg)                   \
  do {                                                                  \
    if (!(expr1 == expr2)) {                                            \
      libMesh::err << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      libmesh_error();                                                  \
    } } while(0)

#define libmesh_assert_not_equal_to_msg(expr1,expr2, msg)               \
  do {                                                                  \
    if (!(expr1 != expr2)) {                                            \
      libMesh::err << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      libmesh_error();                                                  \
    } } while(0)

#define libmesh_assert_less_msg(expr1,expr2, msg)                       \
  do {                                                                  \
    if (!(expr1 < expr2)) {                                             \
      libMesh::err << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      libmesh_error();                                                  \
    } } while(0)

#define libmesh_assert_greater_msg(expr1,expr2, msg)                    \
  do {                                                                  \
    if (!(expr1 > expr2)) {                                             \
      libMesh::err << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      libmesh_error();                                                  \
    } } while(0)

#define libmesh_assert_less_equal_msg(expr1,expr2, msg)                 \
  do {                                                                  \
    if (!(expr1 <= expr2)) {                                            \
      libMesh::err << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      libmesh_error();                                                  \
    } } while(0)

#define libmesh_assert_greater_equal_msg(expr1,expr2, msg)              \
  do {                                                                  \
    if (!(expr1 >= expr2)) {                                            \
      libMesh::err << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      libmesh_error();                                                  \
    } } while(0)
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
// The libmesh_file_error(const std::string& filename) macro prints a message
// and throws a FileError exception
//
// The libmesh_convergence_failure() macro
// throws a ConvergenceFailure exception
#define libmesh_error_msg(msg)                                          \
  do {                                                                  \
    libMesh::err << msg << std::endl;                                   \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, __LIBMESH_DATE__, __LIBMESH_TIME__); \
    LIBMESH_THROW(libMesh::LogicError());                               \
  } while(0)

#define libmesh_error() libmesh_error_msg("")

#define libmesh_exceptionless_error_msg(msg)                            \
  do {                                                                  \
    libMesh::err << msg << std::endl;                                   \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, __LIBMESH_DATE__, __LIBMESH_TIME__); \
    std::terminate();                                                   \
  } while(0)

#define libmesh_exceptionless_error() libmesh_exceptionless_error_msg("")

#define libmesh_not_implemented_msg(msg)                                \
  do {                                                                  \
    libMesh::err << msg << std::endl;                                   \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, __LIBMESH_DATE__, __LIBMESH_TIME__); \
    LIBMESH_THROW(libMesh::NotImplemented());                           \
  } while(0)

#define libmesh_not_implemented() libmesh_not_implemented_msg("")

#define libmesh_file_error_msg(filename, msg)                           \
  do {                                                                  \
    libMesh::MacroFunctions::report_error(__FILE__, __LINE__, __LIBMESH_DATE__, __LIBMESH_TIME__); \
    libMesh::err << msg << std::endl;                                   \
    LIBMESH_THROW(libMesh::FileError(filename));                        \
  } while(0)

#define libmesh_file_error(filename) libmesh_file_error_msg(filename,"")

#define libmesh_convergence_failure()                   \
  do {                                                  \
    LIBMESH_THROW(libMesh::ConvergenceFailure());       \
  } while(0)

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
      libMesh::out << "Configuring libMesh with " #option " is required to run this example." << std::endl; \
      return 77;                                                        \
    } } while(0)

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
                  << __FILE__ << ", line " << __LINE__ << ", compiled " << __LIBMESH_DATE__ << " at " << __LIBMESH_TIME__ << " ***" << std::endl;)
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
#define libmesh_deprecated()                                            \
  libmesh_warning("*** Warning, This code is deprecated, and likely to be removed in future library versions! ");

// A function template for ignoring unused variables.  This is a way
// to shut up unused variable compiler warnings on a case by case
// basis.
template<class T> inline void libmesh_ignore( const T& ) { }


// cast_ref and cast_ptr do a dynamic cast and assert
// the result, if we have RTTI enabled and we're in debug or
// development modes, but they just do a faster static cast if we're
// in optimized mode.
//
// Use these casts when you're certain that a cast will succeed in
// correct code but you want to be able to double-check.
template <typename Tnew, typename Told>
inline Tnew cast_ref(Told& oldvar)
{
#if !defined(NDEBUG) && defined(LIBMESH_HAVE_RTTI) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  try
    {
      Tnew newvar = dynamic_cast<Tnew>(oldvar);
      return newvar;
    }
  catch (std::bad_cast)
    {
      libMesh::err << "Failed to convert " << typeid(Told).name()
                   << " reference to " << typeid(Tnew).name()
                   << std::endl;
      libMesh::err << "The " << typeid(Told).name()
                   << " appears to be a "
                   << typeid(*(&oldvar)).name() << std::endl;
      libmesh_error();
    }
#else
  return(static_cast<Tnew>(oldvar));
#endif
}

template <typename Tnew, typename Told>
inline Tnew libmesh_cast_ref(Told& oldvar)
{
  // we use the less redundantly named libMesh::cast_ref now
  libmesh_deprecated();
  return cast_ref<Tnew>(oldvar);
}

// We use two different function names to avoid an odd overloading
// ambiguity bug with icc 10.1.008
template <typename Tnew, typename Told>
inline Tnew cast_ptr (Told* oldvar)
{
#if !defined(NDEBUG) && defined(LIBMESH_HAVE_RTTI)
  Tnew newvar = dynamic_cast<Tnew>(oldvar);
  if (!newvar)
    {
      libMesh::err << "Failed to convert " << typeid(Told).name()
                   << " pointer to " << typeid(Tnew).name()
                   << std::endl;
      libMesh::err << "The " << typeid(Told).name()
                   << " appears to be a "
                   << typeid(*oldvar).name() << std::endl;
      libmesh_error();
    }
  return newvar;
#else
  return(static_cast<Tnew>(oldvar));
#endif
}


template <typename Tnew, typename Told>
inline Tnew libmesh_cast_ptr (Told* oldvar)
{
  // we use the less redundantly named libMesh::cast_ptr now
  return cast_ptr<Tnew>(oldvar);
}


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


// build a integer representation of version
#define LIBMESH_VERSION_ID(major,minor,patch) (((major) << 16) | ((minor) << 8) | ((patch) & 0xFF))


// Allow for marking functions with "override" if the compiler supports it.
// Note: override ensures that the function is virtual and is
// overriding a virtual function from the base class.
#ifdef LIBMESH_HAVE_CXX11_OVERRIDE
#define libmesh_override override
#else
#define libmesh_override
#endif

} // namespace libMesh


// Backwards compatibility
namespace libMeshEnums
{
using namespace libMesh;
}

#endif // LIBMESH_LIBMESH_COMMON_H
