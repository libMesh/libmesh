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

#ifndef LIBMESH_FE_TRANSFORMATION_BASE_IMPL_H
#define LIBMESH_FE_TRANSFORMATION_BASE_IMPL_H

#include "libmesh/fe_transformation_base.h"
#include "libmesh/h1_fe_transformation.h"
#include "libmesh/hcurl_fe_transformation.h"
#include "libmesh/fe_type.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/fe_base.h"

namespace libMesh
{

template <typename MathType, typename RealType>
std::unique_ptr<FETransformationBase<MathType,RealType>> FETransformationBase<MathType,RealType>::build( const FEType & fe_type )
{
  switch (fe_type.family)
    {
      // H1 Conforming Elements
    case LAGRANGE:
    case HIERARCHIC:
    case BERNSTEIN:
    case SZABAB:
    case CLOUGH: // PB: Really H2
    case HERMITE: // PB: Really H2
    case SUBDIVISION:
    case LAGRANGE_VEC:
    case MONOMIAL: // PB: Shouldn't this be L2 conforming?
    case MONOMIAL_VEC: // PB: Shouldn't this be L2 conforming?
    case XYZ: // PB: Shouldn't this be L2 conforming?
    case RATIONAL_BERNSTEIN:
    case L2_HIERARCHIC: // PB: Shouldn't this be L2 conforming?
    case L2_LAGRANGE: // PB: Shouldn't this be L2 conforming?
    case JACOBI_20_00: // PB: For infinite elements...
    case JACOBI_30_00: // PB: For infinite elements...
      return libmesh_make_unique<H1FETransformation<MathType,RealType>>();

      // HCurl Conforming Elements
    case NEDELEC_ONE:
      return libmesh_make_unique<HCurlFETransformation<MathType,RealType>>();

      // HDiv Conforming Elements
      // L2 Conforming Elements

      // Other...
    case SCALAR:
      // Should never need this for SCALARs
      return libmesh_make_unique<H1FETransformation<MathType,RealType>>();

    default:
      libmesh_error_msg("Unknown family = " << fe_type.family);
    }
}

} // namespace libMesh

#endif // LIBMESH_FE_TRANSFORMATION_BASE_IMPL_H
