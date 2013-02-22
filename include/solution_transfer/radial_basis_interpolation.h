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



#ifndef LIBMESH_RADIAL_BASIS_INTERPOLATION_H
#define LIBMESH_RADIAL_BASIS_INTERPOLATION_H

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/meshfree_interpolation.h"



namespace libMesh
{
/**
 * Radial Basis Function interplation.
 */
template <unsigned int KDDim>
class RadialBasisInterpolation : public InverseDistanceInterpolation<KDDim>
{
  
public:

  /**
   * Constructor. 
   */
  RadialBasisInterpolation () :
    InverseDistanceInterpolation<KDDim> ()
  {}

  /**
   * Clears all internal data structures and restores to a
   * pristine state.
   */
  virtual void clear();

  /**
   * Interpolate source data at target points.
   * Pure virtual, must be overriden in derived classes.
   */
  virtual void interpolate_field_data (const std::vector<std::string> &field_names,
				       const std::vector<Point>  &tgt_pts,
				       std::vector<Number> &tgt_vals) const;

};

} // namespace libMesh


#endif // #define LIBMESH_RADIAL_BASIS_INTERPOLATION_H
