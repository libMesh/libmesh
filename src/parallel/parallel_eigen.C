// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/parallel_eigen.h"
#include "libmesh/int_range.h"

// C++ includes
#include <cstring> // memcpy

// Helper functions in anonymous namespace

namespace
{
  using namespace libMesh;

  static const unsigned int header_size = 2;

  // use "(a+b-1)/b" trick to get a/b to round up
  // static const unsigned int item_per_Real =
  //     (sizeof(Real) + sizeof(unsigned int) - 1) / sizeof(unsigned int);
}

namespace libMesh
{

  namespace Parallel
  {

    template <>
    unsigned int
    Packing<DynamicEigenMatrix>::packable_size(const DynamicEigenMatrix &mtx,
                                               const void *)
    {
      const auto rows = mtx.rows();
      const auto cols = mtx.cols();
      // use "(a+b-1)/b" trick to get a/b to round up
      std::size_t ints_per_data =
          (sizeof(Real) * rows * cols + sizeof(unsigned int) - 1) / sizeof(unsigned int);

      return header_size + ints_per_data;
    }

    template <>
    unsigned int
    Packing<DynamicEigenMatrix>::packed_size(const std::vector<unsigned int>::const_iterator in)
    {
      const auto rows = *in;
      const auto cols = *(in + 1);
      // use "(a+b-1)/b" trick to get a/b to round up
      const unsigned int ints_per_data =
          (sizeof(Real) * rows * cols + sizeof(unsigned int) - 1) / sizeof(unsigned int);

      return header_size + ints_per_data;
    }

    template <>
    unsigned int
    Packing<DynamicEigenMatrix>::packed_size (std::vector<unsigned int>::iterator in)
    {
      return packed_size(std::vector<unsigned int>::const_iterator(in));
    }

    template <>
    void
    Packing<DynamicEigenMatrix>::pack(const DynamicEigenMatrix &mtx,
                                      std::back_insert_iterator<std::vector<unsigned int>> data_out,
                                      const void *)
    {
      const auto rows = mtx.rows();
      const auto cols = mtx.cols();
      *data_out++ = rows;
      *data_out++ = cols;
      // use "(a+b-1)/b" trick to get a/b to round up
      const unsigned int ints_per_data =
          (sizeof(Real) * rows * cols + sizeof(unsigned int) - 1) / sizeof(unsigned int);

      // technically this is undefined behavior (only pointers to byte sized types are allowed in the standard)
      const auto *raw_data = reinterpret_cast<const unsigned int *>(mtx.data());
      for (const auto i : make_range(ints_per_data))
        *data_out++ = (raw_data[i]);
    }

    template <>
    DynamicEigenMatrix
    Packing<DynamicEigenMatrix>::unpack(std::vector<unsigned int>::const_iterator in,
                                        void *)
    {
      const auto rows = *in++;
      const auto cols = *in++;
      DynamicEigenMatrix mtx(rows, cols);

      // use "(a+b-1)/b" trick to get a/b to round up
      const unsigned int ints_per_data =
          (sizeof(Real) * rows * cols + sizeof(unsigned int) - 1) / sizeof(unsigned int);

      auto *raw_data = reinterpret_cast<unsigned int *>(mtx.data());
      for (const auto i : make_range(ints_per_data))
        raw_data[i] = *in++;

      return mtx;
    }

  } // namespace Parallel

} // namespace libMesh
