// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

#ifndef __pool_allocator_h__
#define __pool_allocator_h__

#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_BOOST
#  include <boost/pool/pool.hpp>
#  include <boost/pool/object_pool.hpp>
#  include <boost/pool/pool_alloc.hpp>
#endif

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
  template PoolAllocator : public std::allocator<T>
  {
  public:
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
  template FastPoolAllocator : public std::allocator<T>
  {
  public:
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


#endif // #define __pool_allocator_h__
