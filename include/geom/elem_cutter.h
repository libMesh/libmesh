// The libMesh Finite Element Library.
// Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ELEM_CUTTER_H
#define LIBMESH_ELEM_CUTTER_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <vector>



namespace libMesh
{

// Forward declarations
class Elem;
class SerialMesh;
class TriangleInterface;
class TetGenMeshInterface;
  
/**
 * This class implements cutting a single element into a collection
 * of subelements.
 *
 * \author Benjamin S. Kirk, 2013
 */

// ------------------------------------------------------------
// ElemCutter class definition
class ElemCutter
{
public:

  /**
   * Constructor.  Initializes pointer data
   * without requiring a full \p SerialMesh in this header file.
   */
  ElemCutter();

  /**
   * Destructor.
   */
  ~ElemCutter();

  /**
   * This function implements cutting an element by a signed distance
   * function. The input array \p vertex_distance_func contains the
   * vertex values of a signed distance function, from which the cutting
   * interface is inferred from the 0 level set.  If all vertex values
   * are positive, the element is outside the cutting surface and is not cut.
   * Likewise if all vertex values are negative, the element is inside the
   * cutting surface and is not cut.
   */
  void operator()(const Elem &elem_in,
		  const std::vector<Real> &vertex_distance_func);


  /**
   * Returns a list of in general element pieces considered inside the
   * cutting surface.  These are subelements whose geometric union
   * defines the spatial domain of the inside portion of the cut element.   
   */
  const std::vector<Elem const*> & inside_elements () const
  { return _inside_elem; }
  
  /**
   * Returns a list of in general element pieces considered outside the
   * cutting surface.  These are subelements whose geometric union
   * defines the spatial domain of the outside portion of the cut element.   
   */
  const std::vector<Elem const*> & outside_elements() const
  { return _outside_elem; }
  
protected:

  /**
   * Finds the points where the cutting surface
   * intersects the element edges.
   */
  void find_intersection_points(const Elem &elem,
				const std::vector<Real> &vertex_distance_func);
  
  /**
   * cutting algoritm in 1D.
   */
  void cut_1D();
  
  /**
   * cutting algoritm in 2D.
   */
  void cut_2D(const Elem &elem,
	      const std::vector<Real> &vertex_distance_func);
  
  /**
   * cutting algoritm in 3D.
   */
  void cut_3D();
  
  std::vector<Elem const*> _inside_elem;
  std::vector<Elem const*> _outside_elem;

  AutoPtr<SerialMesh> _inside_mesh_2D;
  AutoPtr<SerialMesh> _outside_mesh_2D;
  AutoPtr<SerialMesh> _inside_mesh_3D;
  AutoPtr<SerialMesh> _outside_mesh_3D;

  AutoPtr<TriangleInterface>   _triangle_inside;
  AutoPtr<TriangleInterface>   _triangle_outside;
  AutoPtr<TetGenMeshInterface> _tetgen_inside;
  AutoPtr<TetGenMeshInterface> _tetgen_outside;

  std::vector<Point> _intersection_pts;
};


} // namespace libMesh


#endif // LIBMESH_ELEM_CUTTER_H

