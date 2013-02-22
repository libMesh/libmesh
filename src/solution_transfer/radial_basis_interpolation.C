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



// C++ includes
#include <iomanip>

// Local includes
#include "libmesh/point.h"
#include "libmesh/radial_basis_interpolation.h"



namespace libMesh
{
  //--------------------------------------------------------------------------------
  // RadialBasisInterpolation methods
  template <unsigned int KDDim>
  void RadialBasisInterpolation<KDDim>::clear()
  {
    // Call base class clear method
    InverseDistanceInterpolation<KDDim>::clear();
  }



  template <unsigned int KDDim>
  void RadialBasisInterpolation<KDDim>::interpolate_field_data (const std::vector<std::string> &field_names,
								const std::vector<Point>  &tgt_pts,
								std::vector<Number> &tgt_vals) const
  {
    libmesh_experimental();
    libmesh_not_implemented();
  }


  
// ------------------------------------------------------------
// Explicit Instantiations
  template class RadialBasisInterpolation<1>;
  template class RadialBasisInterpolation<2>;
  template class RadialBasisInterpolation<3>;

} // namespace libMesh

