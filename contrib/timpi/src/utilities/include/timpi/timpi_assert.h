// The TIMPI Message-Passing Parallelism Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef TIMPI_TIMPI_ASSERT_H
#define TIMPI_TIMPI_ASSERT_H

// The library configuration options
#include "timpi/timpi_config.h"

// C++ includes
#include <functional> // std::less, etc
#include <iostream>
#include <sstream>
#include <type_traits> // std::decay


// Use actual timestamps or constant dummies (to aid ccache)
#ifdef TIMPI_ENABLE_TIMESTAMPS
#  define  TIMPI_TIME __TIME__
#  define  TIMPI_DATE __DATE__
#else
#  define  TIMPI_TIME "notime"
#  define  TIMPI_DATE "nodate"
#endif

namespace TIMPI
{

// Assert macros generally call report_error with __FILE__, __LINE__,
// __DATE__, and __TIME__ in the appropriate order.  This should not
// be called by users directly!
void report_here(const char * file, int line, const char * date, const char * time);
void report_error(const char * file, int line, const char * date, const char * time);

// A function template for ignoring unused variables.  This is a way
// to shut up conditionally unused variable compiler warnings (as are
// often created by assertions)
template<class ...Args> inline void ignore( const Args&... ) { }

// The timpi_dbg_var() macro indicates that an argument to a function
// is used only in debug mode (i.e., when NDEBUG is not defined).
#ifndef NDEBUG
#define timpi_dbg_var(var) var
#else
#define timpi_dbg_var(var)
#endif

// The timpi_assert() macro acts like C's assert(), but throws a
// timpi_error() (including stack trace, etc) instead of just exiting
#ifdef NDEBUG

#define timpi_assert_msg(asserted, msg)  ((void) 0)
#define timpi_assert_equal_to_msg(expr1,expr2, msg)  ((void) 0)
#define timpi_assert_not_equal_to_msg(expr1,expr2, msg)  ((void) 0)
#define timpi_assert_less_msg(expr1,expr2, msg)  ((void) 0)
#define timpi_assert_greater_msg(expr1,expr2, msg)  ((void) 0)
#define timpi_assert_less_equal_msg(expr1,expr2, msg)  ((void) 0)
#define timpi_assert_greater_equal_msg(expr1,expr2, msg)  ((void) 0)

#else

#define timpi_assert_msg(asserted, msg)                               \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      std::cerr << "Assertion `" #asserted "' failed." << std::endl; \
      timpi_error_msg(msg);                                           \
    } } while (0)

#define timpi_assert_equal_to_msg(expr1,expr2, msg)                   \
  do {                                                                  \
    if (!((expr1) == (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_not_equal_to_msg(expr1,expr2, msg)               \
  do {                                                                  \
    if (!((expr1) != (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

// Older Clang has a parsing bug that can't seem to handle the handle the decay statement when
// expanding these macros. We'll just fall back to the previous behavior with those older compilers
#if defined(__clang__) && (__clang_major__ <= 3)

#define timpi_assert_less_msg(expr1,expr2, msg)                       \
  do {                                                                  \
    if (!((expr1) < (expr2))) {                                         \
      std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_greater_msg(expr1,expr2, msg)                    \
  do {                                                                  \
    if (!((expr1) > (expr2))) {                                         \
      std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_less_equal_msg(expr1,expr2, msg)                 \
  do {                                                                  \
    if (!((expr1) <= (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_greater_equal_msg(expr1,expr2, msg)              \
  do {                                                                  \
    if (!((expr1) >= (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

// Newer clangs and all other supported compilers
#else

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

#define timpi_assert_less_msg(expr1,expr2, msg)                       \
  do {                                                                   \
    if (!TIMPI::casting_compare<std::less>()(expr1, expr2)) {         \
      std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_greater_msg(expr1,expr2, msg)                    \
  do {                                                                   \
    if (!TIMPI::casting_compare<std::greater>()(expr1, expr2)) {      \
      std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_less_equal_msg(expr1,expr2, msg)                 \
  do {                                                                   \
    if (!TIMPI::casting_compare<std::less_equal>()(expr1, expr2)) {   \
      std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#define timpi_assert_greater_equal_msg(expr1,expr2, msg)              \
  do {                                                                   \
    if (!TIMPI::casting_compare<std::greater_equal>()(expr1, expr2)) { \
      std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      timpi_error();                                                  \
    } } while (0)

#endif // clang <4.0

#endif


#define timpi_assert(asserted) timpi_assert_msg(asserted, "")
#define timpi_assert_equal_to(expr1,expr2) timpi_assert_equal_to_msg(expr1,expr2, "")
#define timpi_assert_not_equal_to(expr1,expr2) timpi_assert_not_equal_to_msg(expr1,expr2, "")
#define timpi_assert_less(expr1,expr2) timpi_assert_less_msg(expr1,expr2, "")
#define timpi_assert_greater(expr1,expr2) timpi_assert_greater_msg(expr1,expr2, "")
#define timpi_assert_less_equal(expr1,expr2) timpi_assert_less_equal_msg(expr1,expr2, "")
#define timpi_assert_greater_equal(expr1,expr2) timpi_assert_greater_equal_msg(expr1,expr2, "")


// Macro to notate and debug functions which should only be called in
// parallel on every processor at once

#ifndef NDEBUG
#define timpi_parallel_only(comm_obj) do {                            \
    timpi_assert((comm_obj).verify(std::string(__FILE__)));           \
    timpi_assert((comm_obj).verify(std::to_string(__LINE__))); } while (0)
#else
#define timpi_parallel_only(comm_obj)  ((void) 0)
#endif


#ifdef TIMPI_ENABLE_EXCEPTIONS
#define TIMPI_THROW(e) do { throw e; } while (0)
#else
#define TIMPI_THROW(e) do { std::abort(); } while (0)
#endif

// The timpi_error() macro prints a message and throws a LogicError
// exception
//
// The timpi_not_implemented() macro prints a message and throws a
// NotImplemented exception

#define timpi_error_msg(msg)                                               \
  do {                                                                        \
    std::cerr << msg << std::endl;                                            \
    std::stringstream msg_stream;                                             \
    msg_stream << msg;                                                        \
    TIMPI::report_error(__FILE__, __LINE__, TIMPI_DATE, TIMPI_TIME); \
    TIMPI_THROW(std::logic_error(msg_stream.str()));                       \
  } while (0)

#define timpi_error() timpi_error_msg("")

#define timpi_not_implemented_msg(msg)                                     \
  do {                                                                        \
    std::cerr << msg << std::endl;                                            \
    TIMPI::report_error(__FILE__, __LINE__, TIMPI_DATE, TIMPI_TIME); \
    TIMPI_THROW(std::logic_error("Error: not implemented!"));              \
  } while (0)

#define timpi_not_implemented() timpi_not_implemented_msg("")

// The timpi_do_once macro helps us avoid redundant repeated
// repetitions of the same warning messages
#undef timpi_do_once
#define timpi_do_once(do_this)                \
  do {                                          \
    static bool did_this_already = false;       \
    if (!did_this_already) {                    \
      did_this_already = true;                  \
      do_this;                                  \
    } } while (0)


// The timpi_warning macro outputs a file/line/time stamped warning
// message, if warnings are enabled.
#ifdef TIMPI_ENABLE_WARNINGS
#define timpi_warning(message)       \
  timpi_do_once(std::cout << message \
                   << __FILE__ << ", line " << __LINE__ << ", compiled " << TIMPI_DATE << " at " << TIMPI_TIME << " ***" << std::endl;)
#else
#define timpi_warning(message)  ((void) 0)
#endif

// The timpi_experimental macro warns that you are using
// bleeding-edge code
#undef timpi_experimental
#define timpi_experimental()                                          \
  timpi_warning("*** Warning, This code is untested, experimental, or likely to see future API changes: ");


// The timpi_deprecated macro warns that you are using obsoleted code
#undef timpi_deprecated
#ifndef TIMPI_ENABLE_DEPRECATED
#define timpi_deprecated()                                            \
  timpi_error_msg("*** Error, This code is deprecated, and likely to be removed in future library versions! ");
#else
#define timpi_deprecated()                                            \
  timpi_warning("*** Warning, This code is deprecated, and likely to be removed in future library versions! ");
#endif

// A function template for ignoring unused variables.  This is a way
// to shut up unused variable compiler warnings on a case by case
// basis.
template<class ...Args> inline void timpi_ignore( const Args&... ) { }


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
  timpi_assert_equal_to
    (oldvar, static_cast<Told>(static_cast<Tnew>(oldvar)));

  return(static_cast<Tnew>(oldvar));
}

} // namespace TIMPI


#endif // TIMPI_TIMPI_ASSERT_H
