// The libMesh Finite Element Library.
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


#ifndef LIBMESH_PACKING_H
#define LIBMESH_PACKING_H

// libMesh Includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <cstddef>
#include <iterator>
#include <vector>

namespace libMesh
{

namespace Parallel
{

/**
 * Define data types and (un)serialization functions for use when
 * encoding a potentially-variable-size object of type T.
 *
 * Users will need to specialize this class for their particular data
 * types.
 */
template <typename T>
class Packing {
public:
  // Should be an MPI sendable type in specializations, e.g.
  // typedef char buffer_type;
  // typedef unsigned int buffer_type;

  // Should copy an encoding of the provided object into the provided
  // output iterator (which is of type buffer_type)
  template <typename OutputIter, typename Context>
  static void pack(const T & object,
                   OutputIter data_out,
                   const Context * context);

  // Should return the number of array entries (of type buffer_type)
  // required to encode the provided object
  template <typename Context>
  static unsigned int packable_size(const T & object,
                                    const Context * context);

  // Should return the number of array entries which were used to
  // encode the provided serialization of an object which begins at
  // \p iter
  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter);

  // Decode a potentially-variable-size object from a subsequence of a
  // data array, returning a heap-allocated pointer to the result.
  template <typename BufferIter, typename Context>
  static T unpack(BufferIter in, Context * ctx);
};


/**
 * Decode a range of potentially-variable-size objects from a data
 * array.
 */
template <typename Context, typename buffertype,
          typename OutputIter, typename T>
inline void unpack_range (const typename std::vector<buffertype> & buffer,
                          Context * context,
                          OutputIter out,
                          const T * output_type /* used only to infer T */);

/**
 * Encode a range of potentially-variable-size objects to a data
 * array.
 *
 * The data will be buffered in vectors with lengths that do not
 * exceed the sum of \p approx_buffer_size and the size of an
 * individual packed object.
 */
template <typename Context, typename buffertype, typename Iter>
inline Iter pack_range (const Context * context,
                        Iter range_begin,
                        const Iter range_end,
                        typename std::vector<buffertype> & buffer,
                        std::size_t approx_buffer_size = 1000000);

/**
 * Return the total buffer size needed to encode a range of
 * potentially-variable-size objects to a data array.
 */
template <typename Context, typename Iter>
inline std::size_t packed_range_size (const Context * context,
                                      Iter range_begin,
                                      const Iter range_end);

// ------------------------------------------------------------
// Packing member functions, global functions

/**
 * Helper function for range packing
 */
template <typename Context, typename Iter>
inline std::size_t packed_range_size (const Context * context,
                                      Iter range_begin,
                                      const Iter range_end)
{
  typedef typename std::iterator_traits<Iter>::value_type T;

  std::size_t buffer_size = 0;
  for (Iter range_count = range_begin;
       range_count != range_end;
       ++range_count)
    {
      buffer_size += Parallel::Packing<T>::packable_size(*range_count, context);
    }
  return buffer_size;
}


/**
 * Helper function for range packing
 */
template <typename Context, typename buffertype, typename Iter>
inline Iter pack_range (const Context * context,
                        Iter range_begin,
                        const Iter range_end,
                        std::vector<buffertype> & buffer,
                        // When we serialize into buffers, we need to use large buffers to optimize MPI
                        // bandwidth, but not so large as to risk allocation failures.  max_buffer_size
                        // is measured in number of buffer type entries; number of bytes may be 4 or 8
                        // times larger depending on configuration.
                        std::size_t approx_buffer_size)
{
  typedef typename std::iterator_traits<Iter>::value_type T;

  // Count the total size of and preallocate buffer for efficiency.
  // Prepare to stop early if the buffer would be too large.
  std::size_t buffer_size = 0;
  Iter range_stop = range_begin;
  for (; range_stop != range_end && buffer_size < approx_buffer_size;
       ++range_stop)
    {
      std::size_t next_buffer_size =
        Parallel::Packing<T>::packable_size(*range_stop, context);
      buffer_size += next_buffer_size;
    }
  buffer.reserve(buffer.size() + buffer_size);

  // Pack the objects into the buffer
  for (; range_begin != range_stop; ++range_begin)
    {
#ifndef NDEBUG
      std::size_t old_size = buffer.size();
#endif

      Parallel::Packing<T>::pack
        (*range_begin, back_inserter(buffer), context);

#ifndef NDEBUG
      unsigned int my_packable_size =
        Parallel::Packing<T>::packable_size(*range_begin, context);
      unsigned int my_packed_size =
        Parallel::Packing<T>::packed_size (buffer.begin() + old_size);
      libmesh_assert_equal_to (my_packable_size, my_packed_size);
      libmesh_assert_equal_to (buffer.size(), old_size + my_packable_size);
#endif
    }

  return range_stop;
}



/**
 * Helper function for range unpacking
 */
template <typename Context, typename buffertype,
          typename OutputIter, typename T>
inline void unpack_range (const std::vector<buffertype> & buffer,
                          Context * context,
                          OutputIter out_iter,
                          const T * /* output_type */)
{
  // Loop through the buffer and unpack each object, returning the
  // object pointer via the output iterator
  typename std::vector<buffertype>::const_iterator
    next_object_start = buffer.begin();

  while (next_object_start < buffer.end())
    {
      *out_iter++ = Parallel::Packing<T>::unpack(next_object_start, context);
      next_object_start +=
        Parallel::Packing<T>::packed_size(next_object_start);
    }

  // We should have used up the exact amount of data in the buffer
  libmesh_assert (next_object_start == buffer.end());
}


} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_PACKING_H
