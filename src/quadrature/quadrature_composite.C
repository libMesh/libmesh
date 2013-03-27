// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_composite.h"
#include "libmesh/elem.h"



namespace libMesh
{


template <class QSubCell>
void QComposite<QSubCell>::init (const Elem &elem,		     
				 const std::vector<Real> &vertex_distance_func,
				 unsigned int /* p_level */)
{
  libmesh_assert_equal_to (vertex_distance_func.size(), elem.n_vertices());
  
  libmesh_error();
}
  

//--------------------------------------------------------------
// Explicit instantiations
template class QComposite<QGauss>;
  
} // namespace libMesh


