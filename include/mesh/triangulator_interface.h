// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/meshfree_interpolation.h"

// C++ includes
#include <set>
#include <vector>

namespace libMesh
{

// Forward Declarations
class UnstructuredMesh;
template <typename Output> class FunctionBase;
enum ElemType : int;

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
  class MeshedHole;

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
  virtual bool refine_boundary_allowed () const
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
   * Verifying that hole boundaries don't cross the outer boundary or
   * each other is something like O(N_bdys^2*N_points_per_bdy^2), so
   * we only do it if requested
   */
  void set_verify_hole_boundaries(bool v) {_verify_hole_boundaries = v;}

  bool get_verify_hole_boundaries() const {return _verify_hole_boundaries;}

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

  /**
  *  Generate an auto area function based on spacing of boundary points.
  */
  void generate_auto_area_function(const Parallel::Communicator &comm,
                                   const unsigned int num_nearest_pts,
                                   const unsigned int power,
                                   const Number background_value,
                                   const Real  background_eff_dist);
  
  /**
  *  Whether or not an auto area function has been generated.
  */
  bool has_auto_area_function() {return _auto_area_function != nullptr;}
  
  /**
  *  Caluclate the local desired area based on the auto area function.
  */
  Real get_auto_desired_area(const Point &p);

  /**
   * A set of ids to allow on the outer boundary loop: interpreted as
   * boundary ids of 2D elements and/or subdomain ids of 1D edges.  If
   * this is empty, then the outer boundary may be constructed from
   * boundary edges of any id!
   */
  void set_outer_boundary_ids(std::set<std::size_t> bdy_ids) { _bdy_ids = std::move(bdy_ids); }
  const std::set<std::size_t> & get_outer_boundary_ids() const { return _bdy_ids; }

protected:
  /**
   * Helper function to create PSLG segments from our other
   * boundary-defining options (1D mesh edges, 2D mesh
   * boundary sides), if no segments already exist.
   */
  void elems_to_segments();

  /**
   * Helper function to create PSLG segments from our node ordering,
   * up to the maximum node id, if no segments already exist.
   */
  void nodes_to_segments(dof_id_type max_node_id);

  /**
   * Helper function to add extra points (midpoints of initial
   * segments) to a PSLG triangulation
   */
  void insert_any_extra_boundary_points();

  /**
   * Helper function to check holes for intersections if requested.
   */
  void verify_holes(const Hole & outer_bdy);

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
   * A set of ids to allow on the outer boundary loop.
   */
  std::set<std::size_t> _bdy_ids;

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

  /**
   * Flag which tells if we want to check hole geometry
   */
  bool _verify_hole_boundaries;

  /**
  * The auto area function based on the spacing of boundary points
  */
  std::unique_ptr<InverseDistanceInterpolation<3>> _auto_area_function;
};

} // namespace libMesh

#endif // ifndef LIBMESH_MESH_TRIANGULATOR_INTERFACE_H
