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

#ifndef LIBMESH_FE_CONSTRAINTS_IMPL_H
#define LIBMESH_FE_CONSTRAINTS_IMPL_H

#include "fe.h"

namespace libMesh
{
// We put this here such that we can carefully control instantiations for different dimensions.
// Specifically we only want to instantiate this for 2 and 3 dimensions. Otherwise we want
// compile time errors
template <unsigned int Dim, FEFamily T, typename RealType>
void
FE<Dim, T, RealType>::compute_constraints(DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                          const ElemTempl<Real> * elem)
{
  FEShim<Dim, T, RealType>::compute_constraints(constraints, dof_map, variable_number, elem);
}
}

#define INSTANTIATE_FE_CONSTRAINTS0(Dim, RealType)                                                 \
  template void FE<Dim, CLOUGH, RealType>::compute_constraints(                                    \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, HERMITE, RealType>::compute_constraints(                                   \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, HIERARCHIC, RealType>::compute_constraints(                                \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, L2_HIERARCHIC, RealType>::compute_constraints(                             \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, LAGRANGE, RealType>::compute_constraints(                                  \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, LAGRANGE_VEC, RealType>::compute_constraints(                              \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, L2_LAGRANGE, RealType>::compute_constraints(                               \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, MONOMIAL, RealType>::compute_constraints(                                  \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, SCALAR, RealType>::compute_constraints(                                    \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, XYZ, RealType>::compute_constraints(                                       \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, NEDELEC_ONE, RealType>::compute_constraints(                               \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, MONOMIAL_VEC, RealType>::compute_constraints(                              \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *)

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
#define INSTANTIATE_FE_CONSTRAINTS(Dim, RealType)                                                  \
  INSTANTIATE_FE_CONSTRAINTS0(Dim, RealType);                                                      \
  template void FE<Dim, BERNSTEIN, RealType>::compute_constraints(                                 \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, SZABAB, RealType>::compute_constraints(                                    \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *);                               \
  template void FE<Dim, RATIONAL_BERNSTEIN, RealType>::compute_constraints(                        \
      DofConstraints &, DofMap &, const unsigned int, const ElemTempl<Real> *)
#else
#define INSTANTIATE_FE_CONSTRAINTS(Dim, RealType) INSTANTIATE_FE_CONSTRAINTS0(Dim, RealType)
#endif

#endif // LIBMESH_FE_CONSTRAINTS_IMPL_H
