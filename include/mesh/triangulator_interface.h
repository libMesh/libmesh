// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MESH_TRIANGULATOR_INTERFACE_H
#define LIBMESH_MESH_TRIANGULATOR_INTERFACE_H


#include "libmesh/libmesh_config.h"

// Local Includes
#include "libmesh/libmesh.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum ElemType : int;
}
#else
#include "libmesh/enum_elem_type.h"
#endif

// C++ includes
#include <vector>

namespace libMesh
{

// Forward Declarations
class UnstructuredMesh;
template <typename Output> class FunctionBase;

class TriangulatorInterface
{
public:
  /**
   * The constructor.  A reference to the mesh containing the points
   * which are to be triangulated must be provided.  Unless otherwise
   * specified, a convex hull will be computed for the set of input points
   * and the convex hull will be meshed.
   */
  explicit
  TriangulatorInterface(UnstructuredMesh & mesh);

  /**
   * Empty destructor.
   */
  virtual ~TriangulatorInterface() = default;

  /**
   * The TriangulationType is used with the general triangulate function
   * defined below.
   */
  enum TriangulationType
    {
      /**
       * First generate a convex hull from the set of points passed
       * in, and then triangulate this set of points.  This is
       * probably the most common type of usage.
       */
      GENERATE_CONVEX_HULL = 0,

      /**
       * Triangulate the interior of a Planar Straight Line Graph,
       * which is defined implicitly by the order of the "points"
       * vector: a straight line is assumed to lie between each
       * successive pair of points, with an additional line joining
       * the final and first points.
       *
       * Explicitly telling the triangulator to add additional points
       * may be important for this option.
       */
      PSLG = 1,

      /**
       * Does nothing, used as a "null" value.
       */
      INVALID_TRIANGULATION_TYPE
    };

  /**
   * The hole class and its several subclasses define the interface
   * and functionality of a "hole" which appears in a 2D mesh.
   * See mesh_triangle_holes.C/h for definitions.
   */
  class Hole;
  class AffineHole;
  class PolygonHole;
  class ArbitraryHole;

  /**
   * The region class defines the interface
   * and functionality of a "region" which appears in a 2D mesh.
   * See mesh_triangle_holes.C/h for definitions.
   */
  class Region;

  /**
   * This is the main public interface for this function.
   */
  virtual void triangulate() = 0;

  /**
   * Sets and/or gets the desired element type.
   */
  ElemType & elem_type() {return _elem_type;}

  /**
   * Sets and/or gets the desired triangle area. Set to zero to disable
   * area constraint.
   *
   * If a \p desired_area_function is set, then \p desired_area()
   * should be used to set a *minimum* desired area; this will reduce
   * "false negatives" by suggesting how finely to sample \p
   * desired_area_function inside large triangles, where ideally the
   * \p desired_area_function will be satisfied in the triangle
   * interior and not just at the triangle vertices.
   */
  Real & desired_area() {return _desired_area;}

  /**
   * Set a function giving desired triangle area as a function of
   * position.  Set this to nullptr to disable position-dependent area
   * constraint (falling back on desired_area()).
   *
   * This may not be implemented in all subclasses.
   */
  virtual void set_desired_area_function (FunctionBase<Real> *)
  { libmesh_not_implemented(); }

  /**
   * Get the function giving desired triangle area as a function of
   * position, or \p nullptr if no such function has been set.
   */
  virtual FunctionBase<Real> * get_desired_area_function ()
  { return nullptr; }

  /**
   * Sets and/or gets the minimum desired angle. Set to zero to
   * disable angle constraint.
   */
  Real & minimum_angle() {return _minimum_angle;}

  /**
   * Sets and/or gets the desired triangulation type.
   */
  TriangulationType & triangulation_type() {return _triangulation_type;}

  /**
   * Sets and/or gets the flag for inserting add'l points.
   */
  bool & insert_extra_points() {return _insert_extra_points;}

  /**
   * Complicated setter, for compatibility with insert_extra_points()
   */
  void set_interpolate_boundary_points (int n_points);

  /**
   * Complicated getter, for compatibility with insert_extra_points()
   */
  int get_interpolate_boundary_points () const;

  /**
   * Set whether or not the triangulation is allowed to refine the
   * mesh boundary when refining the interior.  This is true by
   * default, but may be set to false to make the mesh boundary more
   * predictable (and so easier to stitch to other meshes) later.
   *
   * This may not be implemented in all subclasses.
   */
  virtual void set_refine_boundary_allowed (bool)
  { libmesh_not_implemented(); }

  /**
   * Get whether or not the triangulation is allowed to refine the
   * mesh boundary when refining the interior.  True by default.
   */
  virtual bool refine_boundary_allowed ()
  { return true; }

  /**
   * Sets/gets flag which tells whether to do two steps of Laplace
   * mesh smoothing after generating the grid.
   */
  bool & smooth_after_generating() {return _smooth_after_generating;}

  /**
   * Whether not to silence internal messages to stdout
   */
  bool & quiet() {return _quiet;}

  /**
   * Attaches a vector of Hole* pointers which will be
   * meshed around.
   */
  void attach_hole_list(const std::vector<Hole*> * holes) {_holes = holes;}

  /**
   * When constructing a PSLG, if the node numbers do not define the
   * desired boundary segments implicitly through the ordering of the
   * points, you can use the segments vector to specify the segments
   * explicitly, Ex: unit square numbered counter-clockwise starting
   * from origin
   * segments[0] = (0,1)
   * segments[1] = (1,2)
   * segments[2] = (2,3)
   * segments[3] = (3,0)
   * (For the above case you could actually use the implicit
   * ordering!)
   */
  std::vector<std::pair<unsigned int, unsigned int>> segments;

  /**
   * Attaches boundary markers.
   * If segments is set, the number of markers must be equal to the size of segments,
   * otherwise, it is equal to the number of points.
   */
  void attach_boundary_marker(const std::vector<int> * markers) { _markers = markers; }

  /**
   * Attaches regions for using attribute to set subdomain IDs and better
   * controlling the triangle sizes within the regions.
   */
  void attach_region_list(const std::vector<Region*> * regions) { _regions = regions; }

protected:
  /**
   * Helper function to add extra points (midpoints of initial
   * segments) to a PSLG triangulation
   */
  void insert_any_extra_boundary_points();

  /**
   * Helper function to count points in and verify holes
   */
  unsigned int total_hole_points();

  /**
   * Reference to the mesh which is to be created by triangle.
   */
  UnstructuredMesh & _mesh;

  /**
   * A pointer to a vector of Hole*s.  If this is nullptr, there
   * are no holes!
   */
  const std::vector<Hole*> * _holes;

  /**
   * Boundary markers
   */
  const std::vector<int> * _markers;

  /**
   * A pointer to a vector of Regions*s.  If this is nullptr, there
   * are no regions!
   */
  const std::vector<Region*> * _regions;

  /**
   * The type of elements to generate.  (Defaults to
   * TRI3).
   */
  ElemType _elem_type;

  /**
   * The desired area for the elements in the resulting mesh.
   */
  Real _desired_area;

  /**
   * Minimum angle in triangles
   */
  Real _minimum_angle;

  /**
   * The type of triangulation to perform: choices are:
   * convex hull
   * PSLG
   */
  TriangulationType _triangulation_type;

  /**
   * Flag which tells whether or not to insert additional nodes
   * before triangulation.  This can sometimes be used to "de-regularize"
   * the resulting triangulation.
   *
   * This flag is supported for backwards compatibility; setting
   * _interpolate_boundary_points = 1 is equivalent.
   */
  bool _insert_extra_points;

  /**
   * Flag which tells how many additional nodes should be inserted
   * between each pair of original mesh points.
   */
  int _interpolate_boundary_points;

  /**
   * Flag which tells whether we should smooth the mesh after
   * it is generated.  True by default.
   */
  bool _smooth_after_generating;

  /**
   * Flag which tells if we want to suppress stdout outputs
   */
  bool _quiet;
};

} // namespace libMesh

#endif // ifndef LIBMESH_MESH_TRIANGULATOR_INTERFACE_H
