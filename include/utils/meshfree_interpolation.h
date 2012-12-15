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



#ifndef LIBMESH_MESHFREE_INTERPOLATION_H
#define LIBMESH_MESHFREE_INTERPOLATION_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"

// C++ includes
#include <string>
#include <vector>



namespace libMesh
{

/**
 * Base class to support various mesh-free interpolation methods.
 * Such methods can be useful for example to pass data between 
 * two different domains which share a physical boundary, where
 * that boundary may be discretized differently in each domain.
 * This is the case for conjugate heat transfer applications where
 * the common interface has overlapping but distinct boundary 
 * discretizations.
 */
class MeshfreeInterpolation
{
public:
  MeshfreeInterpolation () {}

  /**
   * Clears all internal data structures and restores to a
   * pristine state.
   */
  virtual void clear();

  /**
   * The numer of field variables.
   */
  unsigned int n_field_variables () const 
  { return _names.size(); }
  
  /**
   * Sets source data at specified points.
   */     
  virtual void add_field_data (const std::vector<std::string> &field_names,
			       const std::vector<Point>  &pts,
			       const std::vector<Number> &vals);
  /**
   * Interpolate source data at target points.
   * Pure virtual, must be overriden in derived classes.
   */
  virtual void interpolate_field_data (const std::vector<std::string> &field_names,
				       const std::vector<Point>  &tgt_pts,
				       const std::vector<Number> &tgt_vals) const = 0;

protected:
    
  std::vector<std::string> _names;
  std::vector<Point>       _src_pts;
  std::vector<Number>      _src_vals;
};

} // namespace libMesh


#endif // #define LIBMESH_MESHFREE_INTERPOLATION_H
