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

#include "libmesh/fe_map_impl.h"

namespace libMesh
{
template class FEMapTempl<Real>;

// Explicit instantiation of FEMap member functions
template void FEMapTempl<Real>::init_reference_to_physical_map<0>( const std::vector<PointTempl<Real>> &, const ElemTempl<Real> *);
template void FEMapTempl<Real>::init_reference_to_physical_map<1>( const std::vector<PointTempl<Real>> &, const ElemTempl<Real> *);
template void FEMapTempl<Real>::init_reference_to_physical_map<2>( const std::vector<PointTempl<Real>> &, const ElemTempl<Real> *);
template void FEMapTempl<Real>::init_reference_to_physical_map<3>( const std::vector<PointTempl<Real>> &, const ElemTempl<Real> *);

// subdivision elements are implemented only for 2D meshes & reimplement
// the inverse_maps method separately
INSTANTIATE_SUBDIVISION_MAPS(Real);

} // namespace libMesh
