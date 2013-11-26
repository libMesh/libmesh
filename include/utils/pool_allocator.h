// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_POOL_ALLOCATOR_H
#define LIBMESH_POOL_ALLOCATOR_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_BOOST
// See: http://stackoverflow.com/questions/17000542/boost-pool-can-i-wean-it-from-boost-system
#define BOOST_POOL_NO_MT       // disable multi-threading
#define BOOST_THREAD_MUTEX_HPP // define the #include-guard to disable the header
#  include <boost/pool/pool.hpp>
#  include <boost/pool/object_pool.hpp>
#  include <boost/pool/pool_alloc.hpp>
#endif

#include <memory> // std::allocator

namespace libMesh
{
  // If Boost is enabled, wrappers to use their allocators.
#ifdef LIBMESH_HAVE_BOOST

  /**
   * An allocator which can be used in standard containers.  Uses
   * pool-based memory allocation to efficiently allocate many small
   * objects.  Note that object destruction returns memory to the pool
   * rather than deallocate it.  It must be explicitly deallocated
   * prior to program termination.
   */
  template <typename T>
  class PoolAllocator : public boost::pool_allocator<T>
  {
  public:

    /**
     * Methods required for copy construction of containers using this allocator.
     */
    template<typename U>
    struct rebind {
      typedef PoolAllocator<U> other;
    };


    PoolAllocator() :
      boost::pool_allocator<T>()
    {}

    explicit PoolAllocator(const PoolAllocator &o) :
      boost::pool_allocator<T>(o)
    {}

    /**
     * Frees every memory block that doesn't have any allocated chunks.
     * Returns true if at least one memory block was freed.
     */
    static bool release_memory ()
    {
      return boost::singleton_pool<boost::pool_allocator_tag, sizeof(T)>::release_memory();
    }

    /**
     * Frees every memory block. This function invalidates any pointers previously returned
     * by allocation functions. Returns true if at least one memory block was freed.
     */
    static bool purge_memory ()
    {
      return boost::singleton_pool<boost::pool_allocator_tag, sizeof(T)>::purge_memory();
    }
  };



  /**
   * An allocator which can be used in standard containers.  Uses
   * pool-based memory allocation to efficiently allocate many small
   * objects.  Note that object destruction returns memory to the pool
   * rather than deallocate it.  It must be explicitly deallocated
   * prior to program termination.
   */
  template <typename T>
  class FastPoolAllocator : public boost::fast_pool_allocator<T>
  {
  public:

    /**
     * Methods required for copy construction of containers using this allocator.
     */
    template<typename U>
    struct rebind {
      typedef FastPoolAllocator<U> other;
    };


    FastPoolAllocator() :
      boost::fast_pool_allocator<T>()
    {}

    explicit FastPoolAllocator(const FastPoolAllocator &o) :
      boost::fast_pool_allocator<T>(o)
    {}


    /**
     * Frees every memory block that doesn't have any allocated chunks.
     * Returns true if at least one memory block was freed.
     */
    static bool release_memory ()
    {
      return boost::singleton_pool<boost::fast_pool_allocator_tag, sizeof(T)>::release_memory();
    }

    /**
     * Frees every memory block. This function invalidates any pointers previously returned
     * by allocation functions. Returns true if at least one memory block was freed.
     */
    static bool purge_memory ()
    {
      return boost::singleton_pool<boost::fast_pool_allocator_tag, sizeof(T)>::purge_memory();
    }
  };

  // Otherwise fall back to std::allocator<>.
#else

  /**
   * An allocator which can be used in standard containers.
   * A wrapper for \p std::allocator<> when Boost is not available.
   */
  template <typename T>
  class PoolAllocator : public std::allocator<T>
  {
  public:

    /**
     * Methods required for copy construction of containers using this allocator.
     */
    template<typename U>
    struct rebind {
      typedef PoolAllocator<U> other;
    };

    PoolAllocator() :
      std::allocator<T>()
    {}

    explicit PoolAllocator(const PoolAllocator &o) :
      std::allocator<T>(o)
    {}

    /**
     * Frees every memory block that doesn't have any allocated chunks.
     * Returns true if at least one memory block was freed.
     */
    static bool release_memory () { /* no-op for std::allocator<> - already freed. */ return false; }

    /**
     * Frees every memory block. This function invalidates any pointers previously returned
     * by allocation functions. Returns true if at least one memory block was freed.
     */
    static bool purge_memory ()   { /* no-op for std::allocator<> - already freed. */ return false; }
  };



  /**
   * An allocator which can be used in standard containers.
   * A wrapper for \p std::allocator<> when Boost is not available.
   */
  template <typename T>
  class FastPoolAllocator : public std::allocator<T>
  {
  public:

    /**
     * Methods required for copy construction of containers using this allocator.
     */
    template<typename U>
    struct rebind {
      typedef FastPoolAllocator<U> other;
    };

    FastPoolAllocator() :
      std::allocator<T>()
    {}

    explicit FastPoolAllocator(const FastPoolAllocator &o) :
      std::allocator<T>(o)
    {}

    /**
     * Frees every memory block that doesn't have any allocated chunks.
     * Returns true if at least one memory block was freed.
     */
    static bool release_memory () { /* no-op for std::allocator<> - already freed. */ return false; }

    /**
     * Frees every memory block. This function invalidates any pointers previously returned
     * by allocation functions. Returns true if at least one memory block was freed.
     */
    static bool purge_memory ()   { /* no-op for std::allocator<> - already freed. */ return false; }
  };

#endif

} // end namespace libMesh


#endif // LIBMESH_POOL_ALLOCATOR_H
