// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PARALLEL_EIGEN_H
#define LIBMESH_PARALLEL_EIGEN_H

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/int_range.h"

// TIMPI includes
#include "timpi/packing.h"

#ifdef LIBMESH_HAVE_EIGEN

// libEigen includes
#include <Eigen/Core>

namespace libMesh
{

  namespace Parallel
  {
    template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    class Packing<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>
    {
    public:
      typedef unsigned int buffer_type;

      template <typename Context>
      static unsigned int packable_size(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &object,
                                        const Context *context);

      template <typename BufferIter>
      static unsigned int packed_size(BufferIter iter);

      template <typename OutputIter, typename Context>
      static void pack(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &object,
                       OutputIter data_out,
                       const Context *context);

      template <typename BufferIter, typename Context>
      static Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> unpack(BufferIter in, Context *context);

    private:
      template <typename T2>
      static inline constexpr bool IsFixed = TIMPI::StandardType<T2>::is_fixed_type;
      static std::size_t BufferCount(std::size_t n) { return (sizeof(Scalar) * n + sizeof(buffer_type) - 1) / sizeof(buffer_type); }
    };



    template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    template <typename Context>
    unsigned int
    Packing<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>::packable_size(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &mtx,
                                               const Context * context)
    {
      const auto rows = mtx.rows();
      const auto cols = mtx.cols();

      // compute packable size of the underlying data
      std::size_t ints_per_data;
      if constexpr (IsFixed<Scalar>)
      {
        libmesh_ignore(context);
        ints_per_data = BufferCount(rows * cols);
      }
      else
      {
        ints_per_data = 0;
        for (const auto i : make_range(rows * cols))
          ints_per_data += Packing<Scalar>::packable_size(mtx.data()[i], context);
      }

      constexpr std::size_t header_size = (Rows == Eigen::Dynamic) + (Cols == Eigen::Dynamic);
      return header_size + ints_per_data;
    }

    template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    template <typename BufferIter>
    unsigned int
    Packing<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>::packed_size(BufferIter in)
    {
      constexpr std::size_t header_size = (Rows == Eigen::Dynamic) + (Cols == Eigen::Dynamic);
      const std::size_t rows = Rows == Eigen::Dynamic ? *in++ : Rows;
      const std::size_t cols = Cols == Eigen::Dynamic ? *in++ : Cols;

      // compute packable size of the underlying data
      std::size_t ints_per_data;
      if constexpr (IsFixed<Scalar>)
        ints_per_data = BufferCount(rows * cols);
      else
      {
        ints_per_data = 0;
        for (const auto i : make_range(rows * cols))
        {
          ints_per_data += Packing<Scalar>::packed_size(in + ints_per_data);
          libmesh_ignore(i);
        }
      }

      return header_size + ints_per_data;
    }

    template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    template <typename OutputIter, typename Context>
    void
    Packing<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>::pack(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &mtx,
                                      OutputIter data_out,
                                      const Context *context)
    {
      const auto rows = mtx.rows();
      const auto cols = mtx.cols();

      // store dynamic dimensions
      if constexpr (Rows == Eigen::Dynamic)
        *data_out++ = rows;
      if constexpr (Cols == Eigen::Dynamic)
        *data_out++ = cols;

      // pack underlying data
      if constexpr (IsFixed<Scalar>)
      {
        libmesh_ignore(context);
        const auto *raw_data = reinterpret_cast<const buffer_type *>(mtx.data());
        for (const auto i : make_range(BufferCount(rows * cols)))
          *data_out++ = (raw_data[i]);
        }
      else
        for (const auto i : make_range(rows * cols))
          Packing<Scalar>::pack(mtx.data()[i], data_out, context);
    }

    template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    template <typename BufferIter, typename Context>
    Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>
    Packing<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>::unpack(BufferIter in, Context *context)
    {
      const std::size_t rows = Rows == Eigen::Dynamic ? *in++ : Rows;
      const std::size_t cols = Cols == Eigen::Dynamic ? *in++ : Cols;
      auto mtx = Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>(rows, cols);

      // unpack underlying data
      if constexpr (IsFixed<Scalar>)
      {
        libmesh_ignore(context);
        auto *raw_data = reinterpret_cast<buffer_type *>(mtx.data());
        for (const auto i : make_range(BufferCount(rows * cols)))
         raw_data[i] = *in++;
      }
      else
        for (const auto i : make_range(rows * cols))
        {
          // unpack matrix entry and advance iterator
          mtx.data()[i] = Packing<Scalar>::unpack(in, context);
          in += Packing<Scalar>::packed_size(in);
        }

      return mtx;
    }

  } // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_HAVE_EIGEN

#endif // LIBMESH_PARALLEL_EIGEN_H
