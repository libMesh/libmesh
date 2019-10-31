// The PassMess Message-Passing Parallelism Library.
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



#ifndef PASSMESS_PASSMESS_ASSERT_H
#define PASSMESS_PASSMESS_ASSERT_H

// The library configuration options
#include "libmesh/libmesh_config.h"

// C++ includes
#include <functional> // std::less, etc
#include <iostream>
#include <sstream>
#include <type_traits> // std::decay


// Use actual timestamps or constant dummies (to aid ccache)
#ifdef LIBMESH_ENABLE_TIMESTAMPS
#  define  PASSMESS_TIME __TIME__
#  define  PASSMESS_DATE __DATE__
#else
#  define  PASSMESS_TIME "notime"
#  define  PASSMESS_DATE "nodate"
#endif

namespace PassMess
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

// The passmess_dbg_var() macro indicates that an argument to a function
// is used only in debug mode (i.e., when NDEBUG is not defined).
#ifndef NDEBUG
#define passmess_dbg_var(var) var
#else
#define passmess_dbg_var(var)
#endif

// The passmess_assert() macro acts like C's assert(), but throws a
// passmess_error() (including stack trace, etc) instead of just exiting
#ifdef NDEBUG

#define passmess_assert_msg(asserted, msg)  ((void) 0)
#define passmess_assert_equal_to_msg(expr1,expr2, msg)  ((void) 0)
#define passmess_assert_not_equal_to_msg(expr1,expr2, msg)  ((void) 0)
#define passmess_assert_less_msg(expr1,expr2, msg)  ((void) 0)
#define passmess_assert_greater_msg(expr1,expr2, msg)  ((void) 0)
#define passmess_assert_less_equal_msg(expr1,expr2, msg)  ((void) 0)
#define passmess_assert_greater_equal_msg(expr1,expr2, msg)  ((void) 0)

#else

#define passmess_assert_msg(asserted, msg)                               \
  do {                                                                  \
    if (!(asserted)) {                                                  \
      std::cerr << "Assertion `" #asserted "' failed." << std::endl; \
      passmess_error_msg(msg);                                           \
    } } while (0)

#define passmess_assert_equal_to_msg(expr1,expr2, msg)                   \
  do {                                                                  \
    if (!((expr1) == (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_not_equal_to_msg(expr1,expr2, msg)               \
  do {                                                                  \
    if (!((expr1) != (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

// Older Clang has a parsing bug that can't seem to handle the handle the decay statement when
// expanding these macros. We'll just fall back to the previous behavior with those older compilers
#if defined(__clang__) && (__clang_major__ <= 3)

#define passmess_assert_less_msg(expr1,expr2, msg)                       \
  do {                                                                  \
    if (!((expr1) < (expr2))) {                                         \
      std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_greater_msg(expr1,expr2, msg)                    \
  do {                                                                  \
    if (!((expr1) > (expr2))) {                                         \
      std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_less_equal_msg(expr1,expr2, msg)                 \
  do {                                                                  \
    if (!((expr1) <= (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_greater_equal_msg(expr1,expr2, msg)              \
  do {                                                                  \
    if (!((expr1) >= (expr2))) {                                        \
      std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
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

#define passmess_assert_less_msg(expr1,expr2, msg)                       \
  do {                                                                   \
    if (!PassMess::casting_compare<std::less>()(expr1, expr2)) {         \
      std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_greater_msg(expr1,expr2, msg)                    \
  do {                                                                   \
    if (!PassMess::casting_compare<std::greater>()(expr1, expr2)) {      \
      std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_less_equal_msg(expr1,expr2, msg)                 \
  do {                                                                   \
    if (!PassMess::casting_compare<std::less_equal>()(expr1, expr2)) {   \
      std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#define passmess_assert_greater_equal_msg(expr1,expr2, msg)              \
  do {                                                                   \
    if (!PassMess::casting_compare<std::greater_equal>()(expr1, expr2)) { \
      std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << '\n' << msg << std::endl; \
      passmess_error();                                                  \
    } } while (0)

#endif // clang <4.0

#endif


#define passmess_assert(asserted) passmess_assert_msg(asserted, "")
#define passmess_assert_equal_to(expr1,expr2) passmess_assert_equal_to_msg(expr1,expr2, "")
#define passmess_assert_not_equal_to(expr1,expr2) passmess_assert_not_equal_to_msg(expr1,expr2, "")
#define passmess_assert_less(expr1,expr2) passmess_assert_less_msg(expr1,expr2, "")
#define passmess_assert_greater(expr1,expr2) passmess_assert_greater_msg(expr1,expr2, "")
#define passmess_assert_less_equal(expr1,expr2) passmess_assert_less_equal_msg(expr1,expr2, "")
#define passmess_assert_greater_equal(expr1,expr2) passmess_assert_greater_equal_msg(expr1,expr2, "")

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#define PASSMESS_THROW(e) do { throw e; } while (0)
#else
#define PASSMESS_THROW(e) do { std::abort(); } while (0)
#endif

// The passmess_error() macro prints a message and throws a LogicError
// exception
//
// The passmess_not_implemented() macro prints a message and throws a
// NotImplemented exception

#define passmess_error_msg(msg)                                               \
  do {                                                                        \
    std::cerr << msg << std::endl;                                            \
    std::stringstream msg_stream;                                             \
    msg_stream << msg;                                                        \
    PassMess::report_error(__FILE__, __LINE__, PASSMESS_DATE, PASSMESS_TIME); \
    PASSMESS_THROW(std::logic_error(msg_stream.str()));                       \
  } while (0)

#define passmess_error() passmess_error_msg("")

#define passmess_not_implemented_msg(msg)                                     \
  do {                                                                        \
    std::cerr << msg << std::endl;                                            \
    PassMess::report_error(__FILE__, __LINE__, PASSMESS_DATE, PASSMESS_TIME); \
    PASSMESS_THROW(std::logic_error("Error: not implemented!");               \
  } while (0)

#define passmess_not_implemented() passmess_not_implemented_msg("")

// The passmess_do_once macro helps us avoid redundant repeated
// repetitions of the same warning messages
#undef passmess_do_once
#define passmess_do_once(do_this)                \
  do {                                          \
    static bool did_this_already = false;       \
    if (!did_this_already) {                    \
      did_this_already = true;                  \
      do_this;                                  \
    } } while (0)


// The passmess_warning macro outputs a file/line/time stamped warning
// message, if warnings are enabled.
#ifdef LIBMESH_ENABLE_WARNINGS
#define passmess_warning(message)       \
  passmess_do_once(std::cout << message \
                   << __FILE__ << ", line " << __LINE__ << ", compiled " << PASSMESS_DATE << " at " << PASSMESS_TIME << " ***" << std::endl;)
#else
#define passmess_warning(message)  ((void) 0)
#endif

// The passmess_experimental macro warns that you are using
// bleeding-edge code
#undef passmess_experimental
#define passmess_experimental()                                          \
  passmess_warning("*** Warning, This code is untested, experimental, or likely to see future API changes: ");


// The passmess_deprecated macro warns that you are using obsoleted code
#undef passmess_deprecated
#ifndef LIBMESH_ENABLE_DEPRECATED
#define passmess_deprecated()                                            \
  passmess_error_msg("*** Error, This code is deprecated, and likely to be removed in future library versions! ");
#else
#define passmess_deprecated()                                            \
  passmess_warning("*** Warning, This code is deprecated, and likely to be removed in future library versions! ");
#endif

// A function template for ignoring unused variables.  This is a way
// to shut up unused variable compiler warnings on a case by case
// basis.
template<class ...Args> inline void passmess_ignore( const Args&... ) { }


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
  passmess_assert_equal_to
    (oldvar, static_cast<Told>(static_cast<Tnew>(oldvar)));

  return(static_cast<Tnew>(oldvar));
}

} // namespace PassMess


#endif // PASSMESS_PASSMESS_ASSERT_H
