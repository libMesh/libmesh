// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_AUTO_PTR_H
#define LIBMESH_AUTO_PTR_H

#include "libmesh/libmesh_config.h"
#include "libmesh_common.h" // for libmesh_deprecated()
#include "libmesh/safe_bool.h"

// LibMesh's AutoPtr is now libmesh_deprecated(), just like the
// std::auto_ptr it is based on.
//
// New library code should use the UniquePtr typedef (which will
// eventually be a C++11 alias declaration).  LibMesh's UniquePtr is
// one of:
// 1.) std::unique_ptr
// 2.) Howard Hinnant's C++03 compatible boost::unique_ptr
// 3.) The deprecated libMesh AutoPtr
// in that order, depending on what your compiler supports.  If you
// are using a compiler that cannot compile Howard Hinnant's
// unique_ptr implementation, you should probably think about
// upgrading.
#ifdef LIBMESH_ENABLE_UNIQUE_PTR
#ifdef LIBMESH_HAVE_CXX11_UNIQUE_PTR
// If C++11 std::unique_ptr is available, alias declarations are also
// guaranteed to be available.
#  include <memory>
namespace libMesh
{
template<typename T>
using UniquePtr = std::unique_ptr<T>;
}
#elif LIBMESH_HAVE_HINNANT_UNIQUE_PTR
// As per Roy's suggestion, use a combination of a macro and a 'using'
// statement to make libMesh's UniquePtr type equivalent to the
// boost::unique_ptr implementation of Howard Hinnant.
#  include "libmesh/unique_ptr.hpp"
#  define UniquePtr unique_ptr
namespace libMesh
{
// Declare that we are using boost's unique_ptr type
using boost::unique_ptr;
}
#else
#  define UniquePtr AutoPtr
#endif
#else
// libMesh was configured with --disable-unique-ptr, so we'll use
// libMesh's AutoPtr class instead.
#define UniquePtr AutoPtr
#endif

namespace libMesh
{

/**
 *  A wrapper class to provide AutoPtr with reference semantics.  For
 *  example, an AutoPtr can be assigned (or constructed from) the result of
 *  a function which returns an AutoPtr by value.
 *
 *  All the AutoPtrRef stuff should happen behind the scenes.
 */
template<typename Tp1>
struct AutoPtrRef
{
  /**
   * The actual pointer.
   */
  Tp1* _ptr;

  /**
   * Constructor.
   */
  explicit
  AutoPtrRef(Tp1* p)
    : _ptr(p) {}
};


/**
 *  @brief  A simple smart pointer providing strict ownership semantics.
 *
 *  The Standard says:
 *  <pre>
 *  An @c AutoPtr owns the object it holds a pointer to.  Copying an
 *  @c AutoPtr copies the pointer and transfers ownership to the destination.
 *  If more than one @c AutoPtr owns the same object at the same time the
 *  behavior of the program is undefined.
 *
 *  The uses of @c AutoPtr include providing temporary exception-safety for
 *  dynamically allocated memory, passing ownership of dynamically allocated
 *  memory to a function, and returning dynamically allocated memory from a
 *  function.  @c AutoPtr does not meet the CopyConstructible and Assignable
 *  requirements for Standard Library <a href="tables.html#65">container</a>
 *  elements and thus instantiating a Standard Library container with an
 *  @c AutoPtr results in undefined behavior.
 *  </pre>
 *  Quoted from [20.4.5]/3.
 *
 * This class is adopted from the GCC 3.2.1 source tree and should
 * function as a replacement for \p std::auto_ptr<>.  Unfortunately
 * the \p std::auto_ptr<> is not particularly portable since various
 * compilers implement various revisions of the standard.  Using
 * \p AutoPtr<> instead of \p std::auto_ptr<> allows for easy
 * portability.
 *
 * The following are the original copyright declarations distributed with this class:
 *
 * Copyright (C) 2001, 2002 Free Software Foundation, Inc.
 *
 * This file is part of the GNU ISO C++ Library.  This library is free
 * software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this library; see the file COPYING.  If not, write to the Free
 * Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA.
 *
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU General Public License.  This exception does not however
 * invalidate any other reasons why the executable file might be covered by
 * the GNU General Public License.
 *
 * Copyright (c) 1997-1999
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */
template<typename Tp>
class AutoPtr : public safe_bool<AutoPtr<Tp> >
{
private:

  /**
   * The actual dumb pointer this class wraps.
   */
  Tp * _ptr;

public:
  /**
   * The pointed-to type.
   */
  typedef Tp element_type;

  /**
   *  @brief  An %AutoPtr is usually constructed from a raw pointer.
   *  @param  p  A pointer (defaults to NULL).
   *
   *  This object now @e owns the object pointed to by @a p.
   */
  explicit
  AutoPtr(element_type * p = 0)
    : _ptr(p)
  {
    // Note: we can't call libmesh_deprecated() here, since global
    // AutoPtr variables are sometimes created before the libMesh::out
    // stream is ready.
  }

  /**
   *  @brief  An %AutoPtr can be constructed from another %AutoPtr.
   *  @param  a  Another %AutoPtr of the same type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.
   */
  AutoPtr(AutoPtr & a)
    : _ptr(a.release())
  {
  }

  /**
   *  @brief  An %AutoPtr can be constructed from another %AutoPtr.
   *  @param  a  Another %AutoPtr of a different but related type.
   *
   *  A pointer-to-Tp1 must be convertible to a pointer-to-Tp/element_type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.
   */
  template<typename Tp1>
  AutoPtr(AutoPtr<Tp1> & a)
    : _ptr(a.release())
  {
  }

  /**
   *  @brief  %AutoPtr assignment operator.
   *  @param  a  Another %AutoPtr of the same type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.  The object that this one @e used to own and
   *  track has been deleted.
   */
  AutoPtr &
  operator=(AutoPtr & a)
  {
    reset(a.release());
    return *this;
  }

  /**
   *  @brief  %AutoPtr assignment operator.
   *  @param  a  Another %AutoPtr of a different but related type.
   *
   *  A pointer-to-Tp1 must be convertible to a pointer-to-Tp/element_type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.  The object that this one @e used to own and
   *  track has been deleted.
   */
  template <typename Tp1>
  AutoPtr &
  operator=(AutoPtr<Tp1> & a)
  {
    reset(a.release());
    return *this;
  }

  /**
   *  When the %AutoPtr goes out of scope, the object it owns is deleted.
   *  If it no longer owns anything (i.e., @c get() is @c NULL), then this
   *  has no effect.
   *
   *  @if maint
   *  The C++ standard says there is supposed to be an empty throw
   *  specification here, but omitting it is standard conforming.  Its
   *  presence can be detected only if _Tp::~_Tp() throws, but this is
   *  prohibited.  [17.4.3.6]/2
   *  @endif maint
   */
  ~AutoPtr()
  {
    if (!libMesh::warned_about_auto_ptr)
      {
        libMesh::warned_about_auto_ptr = true;
        libmesh_deprecated();
      }
    delete _ptr;
  }

  /**
   *  @brief  Smart pointer dereferencing.
   *
   *  If this %AutoPtr no longer owns anything, then this operation will
   *  crash.  (For a smart pointer, "no longer owns anything" is the same as
   *  being a null pointer, and you know what happens when you dereference
   *  one of those...)
   */
  element_type &
  operator*() const  { return *_ptr; }

  /**
   *  @brief  Smart pointer dereferencing.
   *
   *  This returns the pointer itself, which the language then will
   *  automatically cause to be dereferenced.
   */
  element_type *
  operator->() const  { return _ptr; }

  /**
   *  @brief  Bypassing the smart pointer.
   *  @return  The raw pointer being managed.
   *
   *  You can get a copy of the pointer that this object owns, for
   *  situations such as passing to a function which only accepts a raw
   *  pointer.
   *
   *  @note  This %AutoPtr still owns the memory.
   */
  element_type *
  get() const  { return _ptr; }

  /**
   *  @brief  Bypassing the smart pointer.
   *  @return  The raw pointer being managed.
   *
   *  You can get a copy of the pointer that this object owns, for
   *  situations such as passing to a function which only accepts a raw
   *  pointer.
   *
   *  @note  This %AutoPtr no longer owns the memory.  When this object
   *  goes out of scope, nothing will happen.
   */
  element_type *
  release()
  {
    element_type * tmp = _ptr;
    _ptr = 0;
    return tmp;
  }

  /**
   *  @brief  Forcibly deletes the managed object.
   *  @param  p  A pointer (defaults to NULL).
   *
   *  This object now @e owns the object pointed to by @a p.  The previous
   *  object has been deleted.
   */
  void
  reset(element_type * p = 0)
  {
    if (p != _ptr)
      {
        delete _ptr;
        _ptr = p;
      }
  }

  /** @{
   *  @brief  Automatic conversions
   *
   *  These operations convert an %AutoPtr into and from an AutoPtrRef
   *  automatically as needed.  This allows constructs such as
   *  @code
   *    AutoPtr<Derived>  func_returning_AutoPtr(.....);
   *    ...
   *    AutoPtr<Base> ptr = func_returning_AutoPtr(.....);
   *  @endcode
   */
  AutoPtr(AutoPtrRef<element_type> ref)
    : _ptr(ref._ptr) {}

  /**
   * op= for AutoPtr.  Allows you to write:
   * @code
   * AutoPtr<Base> ptr = func_returning_AutoPtr(.....);
   * @endcode
   */
  AutoPtr &
  operator=(AutoPtrRef<element_type> ref)
  {
    if (ref._ptr != this->get())
      {
        delete _ptr;
        _ptr = ref._ptr;
      }
    return *this;
  }

  /**
   * A "safe" replacement for operator bool () that behaves more like
   * an explicit conversion operator even in C++98. This allows code
   * like if (!foo) to work with AutoPtr.
   */
  bool boolean_test() const
  {
    return (this->get() != NULL);
  }

  /**
   * op() for AutoPtrRef<Tp1>.  Calls the release member.
   */
  template<typename Tp1>
  operator AutoPtrRef<Tp1>()
  { return AutoPtrRef<Tp1>(this->release()); }

  /**
   * op() for AutoPtr<Tp1>.  Calls the release member.
   */
  template<typename Tp1>
  operator AutoPtr<Tp1>()
  { return AutoPtr<Tp1>(this->release()); }
};



} // namespace libMesh

#endif // LIBMESH_AUTO_PTR_H
