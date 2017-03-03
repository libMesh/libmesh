// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_THREADS_ALLOCATORS_H
#define LIBMESH_THREADS_ALLOCATORS_H

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/threads.h"

// Threading building blocks includes
#ifdef LIBMESH_HAVE_TBB_API
#  include "libmesh/ignore_warnings.h"
#  include "tbb/scalable_allocator.h"
#  include "libmesh/restore_warnings.h"
#endif

// C++ includes
#include <memory> // for std::allocator
#include <cstddef>

namespace libMesh
{


/**
 * The Threads namespace is for wrapper functions
 * for common general multithreading algorithms and tasks.
 */
namespace Threads
{
#ifdef LIBMESH_HAVE_TBB_API

//-------------------------------------------------------------------
/**
 * Scalable allocator to be used in multithreaded code chunks which
 * allocate a lot of dynamic memory.  This allocator can be faster
 * than the std::allocator when there are multiple threads.
 */
template <typename T>
class scalable_allocator : public tbb::scalable_allocator<T>
{
public:
  typedef T * pointer;
  typedef const T * const_pointer;
  //     typedef T & reference;              // Intel 7.1 tries to instantiate an allocator<void>,
  //     typedef const T & const_reference;  // so we can't typedef a reference to void.
  typedef T value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  template<typename U>
  struct rebind
  {
    typedef scalable_allocator<U> other;
  };

  scalable_allocator () :
    tbb::scalable_allocator<T>() {}

  scalable_allocator (const scalable_allocator & a) :
    tbb::scalable_allocator<T>(a) {}

  template<typename U>
  scalable_allocator(const scalable_allocator<U> & a) :
    tbb::scalable_allocator<T>(a) {}
};



#else //LIBMESH_HAVE_TBB_API



//-------------------------------------------------------------------
/**
 * Just use std::allocator when the Threading Building Blocks is absent.
 */
template <typename T>
class scalable_allocator : public std::allocator<T>
{
public:
  typedef T * pointer;
  typedef const T * const_pointer;
  //     typedef T & reference;              // Intel 7.1 tries to instantiate an allocator<void>,
  //     typedef const T & const_reference;  // so we can't typedef a reference to void.
  typedef T value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  template<typename U>
  struct rebind
  {
    typedef scalable_allocator<U> other;
  };

  scalable_allocator () :
    std::allocator<T>() {}

  scalable_allocator (const scalable_allocator & a) :
    std::allocator<T>(a) {}

  template<typename U>
  scalable_allocator(const scalable_allocator<U> & a) :
    std::allocator<T>(a) {}
};


#endif // #ifdef LIBMESH_HAVE_TBB_API

} // namespace Threads

} // namespace libMesh

#endif // LIBMESH_THREADS_ALLOCATORS_H
