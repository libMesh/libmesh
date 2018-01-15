// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/quadrature_jacobi.h"

namespace libMesh
{

// See the file: quadrature_jacobi_1D.C for implementation of specific
// element types.

QuadratureType QJacobi::type() const
{
  if ((_alpha == 1) && (_beta == 0))
    return QJACOBI_1_0;

  else if ((_alpha == 2) && (_beta == 0))
    return QJACOBI_2_0;

  else
    libmesh_error_msg("Invalid Jacobi quadrature rule: alpha = " << _alpha << ", beta = " << _beta);
}

}
