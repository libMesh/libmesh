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

#ifndef LIBMESH_FE_OUTPUT_TYPE_H
#define LIBMESH_FE_OUTPUT_TYPE_H

#include "libmesh/libmesh_common.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/tensor_tools.h"

namespace libMesh
{
template <typename> class VectorValue;

/**
 * Most finite element types in libMesh are scalar-valued
 */
template <FEFamily T, typename RealType = Real>
struct FEOutputType
{
  typedef RealType type;
};


/**
 * Specialize for non-scalar-valued elements
 */
template<typename RealType>
struct FEOutputType<LAGRANGE_VEC,RealType>
{
  typedef VectorValue<RealType> type;
};

template<typename RealType>
struct FEOutputType<NEDELEC_ONE,RealType>
{
  typedef VectorValue<RealType> type;
};

template<typename RealType>
struct FEOutputType<MONOMIAL_VEC,RealType>
{
  typedef VectorValue<RealType> type;
};

template <typename MathType, typename RealType = Real>
struct FEOutputTwoTypeBase;

template <typename RealType>
struct FEOutputTwoTypeBase<Real, RealType>
{
  typedef RealType OutputShape;
};

template <typename RealType>
struct FEOutputTwoTypeBase<VectorValue<Real>, RealType>
{
  typedef VectorValue<RealType> OutputShape;
};

template <typename MathType, typename RealType = Real>
struct FEOutputTwoType : public FEOutputTwoTypeBase<MathType, RealType>
{
  using typename FEOutputTwoTypeBase<MathType,RealType>::OutputShape;
  typedef typename TensorTools::IncrementRank<OutputShape>::type          OutputGradient;
  typedef typename TensorTools::IncrementRank<OutputGradient>::type       OutputTensor;
  typedef typename TensorTools::DecrementRank<OutputShape>::type          OutputDivergence;
};

} // namespace libMesh

#endif // LIBMESH_FE_OUTPUT_TYPE_H
