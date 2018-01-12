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

#include "libmesh/tensor_tools.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_n_tensor.h"

namespace libMesh
{
namespace TensorTools
{
// Needed for ExactSolution to compile
Number curl_from_grad(const VectorValue<Number> &)
{
  libmesh_error_msg("Operation not defined for scalar quantities.");
}

VectorValue<Number> curl_from_grad(const TensorValue<Number> & grad)
{
  const Number duz_dy = grad(2,1);
  const Number duy_dz = grad(1,2);
  const Number dux_dz = grad(0,2);
  const Number duz_dx = grad(2,0);
  const Number duy_dx = grad(1,0);
  const Number dux_dy = grad(0,1);

  return VectorValue<Number>(duz_dy - duy_dz,
                             dux_dz - duz_dx,
                             duy_dx - dux_dy);
}

// Needed for ExactSolution to compile. Will implement when needed.
TensorValue<Number> curl_from_grad( const TypeNTensor<3,Number> & /* grad */ )
{
  libmesh_not_implemented();
}

// Needed for ExactSolution to compile
Number div_from_grad( const VectorValue<Number> & /* grad */ )
{
  libmesh_error_msg("Operation not defined for scalar quantities.");
}

Number div_from_grad( const TensorValue<Number> & grad )
{
  const Number dux_dx = grad(0,0);
  const Number duy_dy = grad(1,1);
  const Number duz_dz = grad(2,2);

  return dux_dx + duy_dy + duz_dz;
}

// Needed for ExactSolution to compile. Will implement when needed.
VectorValue<Number> div_from_grad( const TypeNTensor<3,Number> & /* grad */ )
{
  libmesh_not_implemented();
}

}
}
