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



// Local includes
#include "libmesh/fe_scalar_impl.h"
#include "libmesh/fe_scalar_shape_0D_impl.h"
#include "libmesh/fe_scalar_shape_1D_impl.h"
#include "libmesh/fe_scalar_shape_2D_impl.h"
#include "libmesh/fe_scalar_shape_3D_impl.h"

namespace libMesh
{

template struct FEShim<0,SCALAR,Real>;
template struct FEShim<1,SCALAR,Real>;
template struct FEShim<2,SCALAR,Real>;
template struct FEShim<3,SCALAR,Real>;

} // namespace libMesh
