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

#ifndef LIBMESH_PARALLEL_EIGEN_H
#define LIBMESH_PARALLEL_EIGEN_H

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"

// TIMPI includes
#include "timpi/packing.h"

// libEigen includes
#include <Eigen/Core>

namespace libMesh
{

  namespace Parallel
  {
    typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> DynamicEigenMatrix;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> DynamicEigenVector;

    template <>
    class Packing<DynamicEigenMatrix>
    {
    public:
      typedef unsigned int buffer_type;

      template <typename OutputIter, typename Context>
      static void pack(const DynamicEigenMatrix &object,
                       OutputIter data_out,
                       const Context *context);

      template <typename Context>
      static unsigned int packable_size(const DynamicEigenMatrix &object,
                                        const Context *context);

      template <typename BufferIter>
      static unsigned int packed_size(BufferIter iter);

      template <typename BufferIter, typename Context>
      static DynamicEigenMatrix unpack(BufferIter in, Context *ctx);
    };

    // template <typename _Scalar, int _Rows, int _Cols, int _Options>
    // class Packing<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options>>;

  } // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_PARALLEL_EIGEN_H
