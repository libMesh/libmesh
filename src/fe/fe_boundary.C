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

#include "libmesh/fe_boundary_impl.h"
#include "libmesh/fe_xyz_boundary_impl.h"

namespace libMesh
{

// Explicit FEMap Instantiations
FACE_EDGE_SHAPE_ERROR_INSTANTIATE(init_face_shape_functions,Real);
FACE_EDGE_SHAPE_ERROR_INSTANTIATE(init_edge_shape_functions,Real);
INIT_FACE_EDGE_SHAPE_FUNCTIONS_INSTANTIATE(Real);

REINIT_AND_SIDE_MAPS(LAGRANGE,Real);
REINIT_AND_SIDE_MAPS(LAGRANGE_VEC,Real);
REINIT_AND_SIDE_MAPS(L2_LAGRANGE,Real);
REINIT_AND_SIDE_MAPS(HIERARCHIC,Real);
REINIT_AND_SIDE_MAPS(L2_HIERARCHIC,Real);
REINIT_AND_SIDE_MAPS(CLOUGH,Real);
REINIT_AND_SIDE_MAPS(HERMITE,Real);
REINIT_AND_SIDE_MAPS(MONOMIAL,Real);
REINIT_AND_SIDE_MAPS(MONOMIAL_VEC,Real);
REINIT_AND_SIDE_MAPS(SCALAR,Real);
REINIT_AND_SIDE_MAPS(XYZ,Real);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
REINIT_AND_SIDE_MAPS(BERNSTEIN,Real);
REINIT_AND_SIDE_MAPS(SZABAB,Real);
REINIT_AND_SIDE_MAPS(RATIONAL_BERNSTEIN,Real);
#endif

ERRORS_IN_0D_INSTANTIATE(CLOUGH,Real);
ERRORS_IN_0D_INSTANTIATE(HERMITE,Real);
ERRORS_IN_0D_INSTANTIATE(HIERARCHIC,Real);
ERRORS_IN_0D_INSTANTIATE(L2_HIERARCHIC,Real);
ERRORS_IN_0D_INSTANTIATE(LAGRANGE,Real);
ERRORS_IN_0D_INSTANTIATE(L2_LAGRANGE,Real);
ERRORS_IN_0D_INSTANTIATE(LAGRANGE_VEC,Real);
ERRORS_IN_0D_INSTANTIATE(MONOMIAL,Real);
ERRORS_IN_0D_INSTANTIATE(MONOMIAL_VEC,Real);
ERRORS_IN_0D_INSTANTIATE(NEDELEC_ONE,Real);
ERRORS_IN_0D_INSTANTIATE(SCALAR,Real);
SUB_ERRORS_IN_0D_INSTANTIATE(XYZ,Real);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
ERRORS_IN_0D_INSTANTIATE(BERNSTEIN,Real);
ERRORS_IN_0D_INSTANTIATE(SZABAB,Real);
ERRORS_IN_0D_INSTANTIATE(RATIONAL_BERNSTEIN,Real);
#endif

ALL_FAMILY_1D_ERRORS_INSTANTIATE(Real);

template void FE<2,SUBDIVISION,Real>::reinit(ElemTempl<Real> const *, unsigned int, Real, const std::vector<PointTempl<Real>> * const, const std::vector<Real> * const);
template void FE<2,NEDELEC_ONE,Real>::reinit(ElemTempl<Real> const *, unsigned int, Real, const std::vector<PointTempl<Real>> * const, const std::vector<Real> * const);
template void FE<2,NEDELEC_ONE,Real>::side_map(ElemTempl<Real> const *, ElemTempl<Real> const *, const unsigned int, const std::vector<PointTempl<Real>> &, std::vector<PointTempl<Real>> &);
template void FE<2,NEDELEC_ONE,Real>::edge_reinit(ElemTempl<Real> const *, unsigned int, Real, const std::vector<PointTempl<Real>> * const, const std::vector<Real> * const);

template void FE<3,NEDELEC_ONE,Real>::reinit(ElemTempl<Real> const *, unsigned int, Real, const std::vector<PointTempl<Real>> * const, const std::vector<Real> * const);
template void FE<3,NEDELEC_ONE,Real>::side_map(ElemTempl<Real> const *, ElemTempl<Real> const *, const unsigned int, const std::vector<PointTempl<Real>> &, std::vector<PointTempl<Real>> &);
template void FE<3,NEDELEC_ONE,Real>::edge_reinit(ElemTempl<Real> const *, unsigned int, Real, const std::vector<PointTempl<Real>> * const, const std::vector<Real> * const);

template void FEMapTempl<Real>::compute_face_map(int, const std::vector<Real> &, const ElemTempl<Real> *);
template void FEMapTempl<Real>::compute_edge_map(int, const std::vector<Real> &, const ElemTempl<Real> *);

// Intel 9.1 complained it needed this in devel mode.
//template void FE<2,XYZ>::init_face_shape_functions(const std::vector<Point> &, const Elem *);
//template void FE<3,XYZ>::init_face_shape_functions(const std::vector<Point> &, const Elem *);

} // namespace libMesh
