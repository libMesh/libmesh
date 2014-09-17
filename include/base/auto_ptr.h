// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// LibMesh's AutoPtr was once equivalent to the (currently deprecated)
// std::auto_ptr, it is now either std::unique_ptr or Howard Hinnant's
// C++03 compatible boost::unique_ptr, depending on your compiler's
// capabilities.

#ifdef LIBMESH_HAVE_CXX11_UNIQUE_PTR
#  include <memory>
#  define AutoPtr std::unique_ptr
#elif LIBMESH_HAVE_HINNANT_UNIQUE_PTR
#  include "libmesh/unique_ptr.hpp"
#  define AutoPtr boost::unique_ptr
#else
#  define AutoPtr libMesh::DeprecatedAutoPtr
#endif


namespace libMesh
{

/**
 *  A wrapper class to provide DeprecatedAutoPtr with reference semantics.  For
 *  example, an DeprecatedAutoPtr can be assigned (or constructed from) the result of
 *  a function which returns an DeprecatedAutoPtr by value.
 *
 *  All the DeprecatedAutoPtrRef stuff should happen behind the scenes.
 */
template<typename Tp1>
struct DeprecatedAutoPtrRef
{
  /**
   * The actual pointer.
   */
  Tp1* _ptr;

  /**
   * Constructor.
   */
  explicit
  DeprecatedAutoPtrRef(Tp1* p)
    : _ptr(p) {}
};


/**
 *  @brief  A simple smart pointer providing strict ownership semantics.
 *
 *  The Standard says:
 *  <pre>
 *  An @c DeprecatedAutoPtr owns the object it holds a pointer to.  Copying an
 *  @c DeprecatedAutoPtr copies the pointer and transfers ownership to the destination.
 *  If more than one @c DeprecatedAutoPtr owns the same object at the same time the
 *  behavior of the program is undefined.
 *
 *  The uses of @c DeprecatedAutoPtr include providing temporary exception-safety for
 *  dynamically allocated memory, passing ownership of dynamically allocated
 *  memory to a function, and returning dynamically allocated memory from a
 *  function.  @c DeprecatedAutoPtr does not meet the CopyConstructible and Assignable
 *  requirements for Standard Library <a href="tables.html#65">container</a>
 *  elements and thus instantiating a Standard Library container with an
 *  @c DeprecatedAutoPtr results in undefined behavior.
 *  </pre>
 *  Quoted from [20.4.5]/3.
 *
 * This class is adopted from the GCC 3.2.1 source tree and should
 * function as a replacement for \p std::auto_ptr<>.  Unfortunately
 * the \p std::auto_ptr<> is not particularly portable since various
 * compilers implement various revisions of the standard.  Using
 * \p DeprecatedAutoPtr<> instead of \p std::auto_ptr<> allows for easy
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
class DeprecatedAutoPtr
{
private:

  /**
   * The actual dumb pointer this class wraps.
   */
  Tp* _ptr;

public:
  /**
   * The pointed-to type.
   */
  typedef Tp element_type;

  /**
   *  @brief  An %DeprecatedAutoPtr is usually constructed from a raw pointer.
   *  @param  p  A pointer (defaults to NULL).
   *
   *  This object now @e owns the object pointed to by @a p.
   */
  explicit
  DeprecatedAutoPtr(element_type* p = 0)
    : _ptr(p) {}

  /**
   *  @brief  An %DeprecatedAutoPtr can be constructed from another %DeprecatedAutoPtr.
   *  @param  a  Another %DeprecatedAutoPtr of the same type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.
   */
  DeprecatedAutoPtr(DeprecatedAutoPtr& a)
    : _ptr(a.release()) {}

  /**
   *  @brief  An %DeprecatedAutoPtr can be constructed from another %DeprecatedAutoPtr.
   *  @param  a  Another %DeprecatedAutoPtr of a different but related type.
   *
   *  A pointer-to-Tp1 must be convertible to a pointer-to-Tp/element_type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.
   */
  template<typename Tp1>
  DeprecatedAutoPtr(DeprecatedAutoPtr<Tp1>& a)
    : _ptr(a.release()) {}

  /**
   *  @brief  %DeprecatedAutoPtr assignment operator.
   *  @param  a  Another %DeprecatedAutoPtr of the same type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.  The object that this one @e used to own and
   *  track has been deleted.
   */
  DeprecatedAutoPtr&
  operator=(DeprecatedAutoPtr& a)
  {
    reset(a.release());
    return *this;
  }

  /**
   *  @brief  %DeprecatedAutoPtr assignment operator.
   *  @param  a  Another %DeprecatedAutoPtr of a different but related type.
   *
   *  A pointer-to-Tp1 must be convertible to a pointer-to-Tp/element_type.
   *
   *  This object now @e owns the object previously owned by @a a, which has
   *  given up ownsership.  The object that this one @e used to own and
   *  track has been deleted.
   */
  template <typename Tp1>
  DeprecatedAutoPtr&
  operator=(DeprecatedAutoPtr<Tp1>& a)
  {
    reset(a.release());
    return *this;
  }

  /**
   *  When the %DeprecatedAutoPtr goes out of scope, the object it owns is deleted.
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
  ~DeprecatedAutoPtr() { delete _ptr; }

  /**
   *  @brief  Smart pointer dereferencing.
   *
   *  If this %DeprecatedAutoPtr no longer owns anything, then this operation will
   *  crash.  (For a smart pointer, "no longer owns anything" is the same as
   *  being a null pointer, and you know what happens when you dereference
   *  one of those...)
   */
  element_type&
  operator*() const  { return *_ptr; }

  /**
   *  @brief  Smart pointer dereferencing.
   *
   *  This returns the pointer itself, which the language then will
   *  automatically cause to be dereferenced.
   */
  element_type*
  operator->() const  { return _ptr; }

  /**
   *  @brief  Bypassing the smart pointer.
   *  @return  The raw pointer being managed.
   *
   *  You can get a copy of the pointer that this object owns, for
   *  situations such as passing to a function which only accepts a raw
   *  pointer.
   *
   *  @note  This %DeprecatedAutoPtr still owns the memory.
   */
  element_type*
  get() const  { return _ptr; }

  /**
   *  @brief  Bypassing the smart pointer.
   *  @return  The raw pointer being managed.
   *
   *  You can get a copy of the pointer that this object owns, for
   *  situations such as passing to a function which only accepts a raw
   *  pointer.
   *
   *  @note  This %DeprecatedAutoPtr no longer owns the memory.  When this object
   *  goes out of scope, nothing will happen.
   */
  element_type*
  release()
  {
    element_type* tmp = _ptr;
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
  reset(element_type* p = 0)
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
   *  These operations convert an %DeprecatedAutoPtr into and from an DeprecatedAutoPtrRef
   *  automatically as needed.  This allows constructs such as
   *  @code
   *    DeprecatedAutoPtr<Derived>  func_returning_DeprecatedAutoPtr(.....);
   *    ...
   *    DeprecatedAutoPtr<Base> ptr = func_returning_DeprecatedAutoPtr(.....);
   *  @endcode
   */
  DeprecatedAutoPtr(DeprecatedAutoPtrRef<element_type> ref)
    : _ptr(ref._ptr) {}

  /**
   * op= for DeprecatedAutoPtr.  Allows you to write:
   * @code
   * DeprecatedAutoPtr<Base> ptr = func_returning_DeprecatedAutoPtr(.....);
   * @endcode
   */
  DeprecatedAutoPtr&
  operator=(DeprecatedAutoPtrRef<element_type> ref)
  {
    if (ref._ptr != this->get())
      {
        delete _ptr;
        _ptr = ref._ptr;
      }
    return *this;
  }

  /**
   * op() for DeprecatedAutoPtrRef<Tp1>.  Calls the release member.
   */
  template<typename Tp1>
  operator DeprecatedAutoPtrRef<Tp1>()
  { return DeprecatedAutoPtrRef<Tp1>(this->release()); }

  /**
   * op() for DeprecatedAutoPtr<Tp1>.  Calls the release member.
   */
  template<typename Tp1>
  operator DeprecatedAutoPtr<Tp1>()
  { return DeprecatedAutoPtr<Tp1>(this->release()); }
};



} // namespace libMesh

#endif // LIBMESH_AUTO_PTR_H
