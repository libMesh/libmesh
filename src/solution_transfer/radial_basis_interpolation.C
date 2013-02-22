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
#include "libmesh/radial_basis_interpolation.h"
#include "libmesh/radial_basis_functions.h"
#include "libmesh/libmesh_logging.h"



namespace libMesh
{
  //--------------------------------------------------------------------------------
  // RadialBasisInterpolation methods
  template <unsigned int KDDim, class RBF>
  void RadialBasisInterpolation<KDDim,RBF>::clear()
  {
    // Call base class clear method
    InverseDistanceInterpolation<KDDim>::clear();
  }



  template <unsigned int KDDim, class RBF>
  void RadialBasisInterpolation<KDDim,RBF>::prepare_for_use()
  {
    // Call base class methods for prep
    InverseDistanceInterpolation<KDDim>::prepare_for_use();
    InverseDistanceInterpolation<KDDim>::construct_kd_tree();

#ifndef LIBMESH_HAVE_EIGEN
    
    std::err << "ERROR: this functionality presently requires Eigen!\n";
    libmesh_error();
    
#else
    START_LOG ("prepare_for_use()", "RadialBasisInterpolation<>");

      
    STOP_LOG  ("prepare_for_use()", "RadialBasisInterpolation<>");
 #endif
    
 }



  template <unsigned int KDDim, class RBF>
  void RadialBasisInterpolation<KDDim,RBF>::interpolate_field_data (const std::vector<std::string> & /* field_names */,
								    const std::vector<Point>  & /* tgt_pts */,
								    std::vector<Number> & /* tgt_vals */) const
  {
    libmesh_experimental();
    libmesh_not_implemented();

    RBF rbf;
  }


  
// ------------------------------------------------------------
// Explicit Instantiations
  template class RadialBasisInterpolation<3, WendlandRBF<3,0> >;
  template class RadialBasisInterpolation<3, WendlandRBF<3,2> >;
  template class RadialBasisInterpolation<3, WendlandRBF<3,4> >;
  template class RadialBasisInterpolation<3, WendlandRBF<3,8> >;

} // namespace libMesh

