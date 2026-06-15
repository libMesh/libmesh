// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/libmesh_config.h"

// libmesh includes
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_triangle_holes.h"
#include "libmesh/mesh_triangle_wrapper.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/utility.h"

#include "libmesh/meshfree_interpolation.h"

// C/C++ includes
#include <limits>
#include <sstream>


namespace libMesh
{
//
// Function definitions for the AutoAreaFunction class
//

// Constructor
AutoAreaFunction::AutoAreaFunction (const Parallel::Communicator &comm,
                                    const unsigned int num_nearest_pts,
                                    const unsigned int power,
                                    const Real background_value,
                                    const Real  background_eff_dist):
  _comm(comm),
  _num_nearest_pts(num_nearest_pts),
  _power(power),
  _background_value(background_value),
  _background_eff_dist(background_eff_dist),
  _auto_area_mfi(std::make_unique<InverseDistanceInterpolation<3>>(_comm, _num_nearest_pts, _power, _background_value, _background_eff_dist))
{
  this->_initialized = false;
  this->_is_time_dependent = false;
}

// Destructor
AutoAreaFunction::~AutoAreaFunction () = default;

void AutoAreaFunction::init_mfi (const std::vector<Point> & input_pts,
                                 const std::vector<Real> & input_vals)
{
  std::vector<std::string> field_vars{"f"};
  _auto_area_mfi->set_field_variables(field_vars);
  _auto_area_mfi->get_source_points() = input_pts;
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  std::vector<Number> input_complex_vals;
  for (const auto & input_val : input_vals)
    input_complex_vals.push_back(Complex (input_val, 0.0));
  _auto_area_mfi->get_source_vals() = input_complex_vals;
#else
  _auto_area_mfi->get_source_vals() = input_vals;
#endif
  _auto_area_mfi->prepare_for_use();
  this->_initialized = true;
}

Real AutoAreaFunction::operator() (const Point & p,
                                   const Real /*time*/)
{
  libmesh_assert(this->_initialized);

  std::vector<Point> target_pts;
  std::vector<Number> target_vals;

  target_pts.push_back(p);
  target_vals.resize(1);

  _auto_area_mfi->interpolate_field_data(_auto_area_mfi->field_variables(), target_pts, target_vals);

  return libmesh_real(target_vals.front());
}

//
// Function definitions for the TriangulatorInterface class
//

// Constructor
TriangulatorInterface::TriangulatorInterface(UnstructuredMesh & mesh)
  : _mesh(mesh),
    _holes(nullptr),
    _markers(nullptr),
    _regions(nullptr),
    _elem_type(TRI3),
    _desired_area(0.1),
    _minimum_angle(20.0),
    _triangulation_type(GENERATE_CONVEX_HULL),
    _insert_extra_points(false),
    _smooth_after_generating(true),
    _quiet(true),
    _fixup_tri7_center_nodes(false),
    _auto_area_function(nullptr)
{}


void TriangulatorInterface::set_interpolate_boundary_points (int n_points)
{
  // Maybe we'll reserve a meaning for negatives later?
  libmesh_assert(n_points >= 0);

  _interpolate_boundary_points = n_points;

  // backwards compatibility - someone (including us) might want to
  // query this via the old API.
  _insert_extra_points = n_points;
}



int TriangulatorInterface::get_interpolate_boundary_points () const
{
  // backwards compatibility - someone might have turned this off via
  // the old API
  if (!_insert_extra_points)
    return 0;

  return _interpolate_boundary_points;
}



void TriangulatorInterface::elems_to_segments()
{
  // Don't try to override manually specified segments
  if (!this->segments.empty())
    return;

  // If we have edges, they should form the polyline with the ordering
  // we want.  Let's turn them into segments for later use, because
  // we're going to delete the original elements to replace with our
  // triangulation.
  if (_mesh.n_elem())
    {
      // Mapping from points to node ids, to back those out from
      // MeshedHole results later
      std::map<Point, dof_id_type> point_id_map;

      for (Node * node : _mesh.node_ptr_range())
        {
          // We're not going to support overlapping nodes on the boundary
          libmesh_error_msg_if
            (point_id_map.count(*node),
             "TriangulatorInterface does not support overlapping nodes found at "
             << static_cast<Point&>(*node));

          point_id_map.emplace(*node, node->id());
        }

      // We don't support directly generating Tri6, so for
      // compatibility with future stitching we need to be working
      // with first-order elements.  Let's get rid of any non-vertex
      // nodes we just added.
      for (Elem * elem : _mesh.element_ptr_range())
        for (auto n : make_range(elem->n_vertices(), elem->n_nodes()))
          point_id_map.erase(elem->point(n));

      // We'll steal the ordering calculation from
      // the MeshedHole code
      const TriangulatorInterface::MeshedHole mh { _mesh, this->_bdy_ids };

      // If we've specified only a subset of the mesh as our outer
      // boundary, then we may have nodes that don't actually fall
      // inside that boundary.  Triangulator code doesn't like Steiner
      // points that aren't inside the triangulation domain, so we
      // need to get rid of them.
      //
      // Also, if we're using Edge3 elements to define our outer
      // boundary, we're only dealing with their 2 end nodes and we'll
      // need to get rid of their central nodes.
      std::unordered_set<Node *> nodes_to_delete;

      for (Elem * elem : _mesh.element_ptr_range())
        for (auto n : make_range(elem->n_vertices(), elem->n_nodes()))
          nodes_to_delete.insert(elem->node_ptr(n));

      if (!this->_bdy_ids.empty())
        {
          for (auto & node : _mesh.node_ptr_range())
            if (!mh.contains(*node))
              nodes_to_delete.insert(node);
        }

      // And now we're done with elements.  Delete them lest they have
      // dangling pointers to nodes we'll be deleting.
      _mesh.clear_elems();

      // Make segments from boundary nodes; also make sure we don't
      // delete them.
      const std::size_t np = mh.n_points();
      for (auto i : make_range(np))
        {
          const Point pt = mh.point(i);
          const dof_id_type id0 = libmesh_map_find(point_id_map, pt);
          nodes_to_delete.erase(_mesh.node_ptr(id0));
          const Point next_pt = mh.point((np+i+1)%np);
          const dof_id_type id1 = libmesh_map_find(point_id_map, next_pt);
          this->segments.emplace_back(id0, id1);
          for (auto m : make_range(mh.n_midpoints()))
          {
            this->segment_midpoints.emplace_back(mh.midpoint(m, i));
            this->segment_midpoints_keys.emplace_back(pt);
          }
        }

      for (Node * node : nodes_to_delete)
        _mesh.delete_node(node);

      if (this->_verify_hole_boundaries && _holes)
        this->verify_holes(mh);
    }
}



void TriangulatorInterface::nodes_to_segments(dof_id_type max_node_id)
{
  // Don't try to override manually specified segments, or try to add
  // segments if we're doing a convex hull
  if (!this->segments.empty() || _triangulation_type != PSLG)
    return;

  for (auto node_it = _mesh.nodes_begin(),
       node_end = _mesh.nodes_end();
       node_it != node_end;)
    {
      Node * node = *node_it;

      // If we're out of boundary nodes, the rest are going to be
      // Steiner points or hole points
      if (node->id() >= max_node_id)
        break;

      ++node_it;

      Node * next_node = (node_it == node_end) ?
        *_mesh.nodes_begin() : *node_it;

      this->segments.emplace_back(node->id(), next_node->id());
    }

  if (this->_verify_hole_boundaries && _holes)
    {
      std::vector<Point> outer_pts;
      for (auto segment : this->segments)
        outer_pts.push_back(_mesh.point(segment.first));

      ArbitraryHole ah(outer_pts);
      this->verify_holes(ah);
    }
}



void TriangulatorInterface::insert_any_extra_boundary_points()
{
  // If the initial PSLG is really simple, e.g. an L-shaped domain or
  // a square/rectangle, the resulting triangulation may be very
  // "structured" looking.  Sometimes this is a problem if your
  // intention is to work with an "unstructured" looking grid.  We can
  // attempt to work around this limitation by inserting midpoints
  // into the original PSLG.  Inserting additional points into a
  // set of points meant to be a convex hull usually makes less sense.

  const int n_interpolated = this->get_interpolate_boundary_points();
  if ((_triangulation_type==PSLG) && n_interpolated)
    {
      // If we were lucky enough to start with contiguous node ids,
      // let's keep them that way.
      dof_id_type nn = _mesh.max_node_id();

      std::vector<std::pair<unsigned int, unsigned int>> old_segments =
        std::move(this->segments);

      // We expect to have converted any elems and/or nodes into
      // segments by now.
      libmesh_assert(!old_segments.empty());

      this->segments.clear();

      // Insert a new point on each segment at evenly spaced locations
      // between existing boundary points.
      // np=index into new points vector
      // n =index into original points vector
      for (auto old_segment : old_segments)
        {
          Node * begin_node = _mesh.node_ptr(old_segment.first);
          Node * end_node = _mesh.node_ptr(old_segment.second);
          dof_id_type current_id = begin_node->id();
          for (auto i : make_range(n_interpolated))
            {
              // new points are equispaced along the original segments
              const Point new_point =
                ((n_interpolated-i) * *(Point *)(begin_node) +
                 (i+1) * *(Point *)(end_node)) /
                (n_interpolated + 1);
              Node * next_node = _mesh.add_point(new_point, nn++);
              this->segments.emplace_back(current_id,
                                          next_node->id());
              current_id = next_node->id();
            }
          this->segments.emplace_back(current_id,
                                      end_node->id());
        }
    }
}


void TriangulatorInterface::increase_triangle_order()
{
  switch (_elem_type)
    {
    case TRI3:
    // Nothing to do if we're not requested to increase order
      return;
    case TRI6:
      _mesh.all_second_order();
      break;
    case TRI7:
      _mesh.all_complete_order();
      break;
    default:
      libmesh_not_implemented();
    }

  // If we have any midpoint location data, we'll want to look it up
  // by point.  all_midpoints[{p, m}] will be the mth midpoint
  // location following after point p (when traversing a triangle
  // counter-clockwise)
  std::map<std::pair<Point, unsigned int>, Point> all_midpoints;
  unsigned int n_midpoints =
    this->segment_midpoints.size() / this->segments.size();
  libmesh_assert_equal_to(this->segments.size() * n_midpoints,
                          this->segment_midpoints.size());
  for (auto m : make_range(n_midpoints))
    for (auto i : make_range(this->segments.size()))
      {
        const Point & p = segment_midpoints_keys[i*n_midpoints+m];
        all_midpoints[{p,m}] =
          this->segment_midpoints[i*n_midpoints+m];
      }

  if (_holes)
    for (const Hole * hole : *_holes)
      {
        if (!hole->n_midpoints())
          continue;
        if (!n_midpoints)
          n_midpoints = hole->n_midpoints();
        else if (hole->n_midpoints() != n_midpoints)
          libmesh_not_implemented_msg
            ("Differing boundary midpoint counts " <<
             hole->n_midpoints() << " and " << n_midpoints);

        // Our inner holes are expected to have points in
        // counter-clockwise order, which is backwards from how we
        // want to traverse them when iterating in counter-clockwise
        // order over a triangle, so we'll need to reverse our maps
        // carefully here.
        const auto n_hole_points = hole->n_points();
        libmesh_assert(n_hole_points);
        for (auto m : make_range(n_midpoints))
          {
            for (auto i : make_range(n_hole_points-1))
              {
                const Point & p = hole->point(i+1);
                all_midpoints[{p,m}] = hole->midpoint(n_midpoints-m-1, i);
              }
            const Point & p = hole->point(0);
            all_midpoints[{p,m}] = hole->midpoint(n_midpoints-m-1, n_hole_points-1);
          }
      }

  // The n_midpoints > 1 case is for future proofing, but in the
  // present we have EDGE4 and no TRI10 yet.
  if (n_midpoints > 1)
    libmesh_not_implemented_msg
      ("Cannot construct triangles with more than 1 midpoint per edge");

  if (!n_midpoints)
    return;

  for (Elem * elem : _mesh.element_ptr_range())
    {
      // This should only be called right after we've finished
      // converting a triangulation to higher order
      libmesh_assert_equal_to(elem->n_vertices(), 3);
      libmesh_assert_not_equal_to(elem->default_order(), FIRST);

      for (auto n : make_range(3))
        {
          // Only hole/outer boundary segments need adjusted midpoints
          if (elem->neighbor_ptr(n))
            continue;

          const Point & p = elem->point(n);

          if (const auto it = all_midpoints.find({p,0});
              it != all_midpoints.end())
            elem->point(n+3) = it->second;
        }
    }

  // Moving boundary mid-edge nodes can displace the TRI7 interior node
  // and tangle the element map.  Repositioning the interior node is
  // opt-in (off by default); the validity check always runs.
  if (_elem_type == TRI7 && _fixup_tri7_center_nodes)
    this->fixup_tri7_center_nodes();

  this->verify_quadratic_elements();
}


void TriangulatorInterface::fixup_tri7_center_nodes()
{
  libmesh_assert_equal_to(_elem_type, TRI7);

  // Place the interior node at the image of the reference centroid
  // (xi, eta) = (1/3, 1/3) under the curved Tri6 map, using the Tri6
  // shape function values there as weights: -1/9 on the vertices and
  // 4/9 on the mid-edges.  This reduces to the straight-edge centroid
  // when no boundary midpoint has moved.
  static const Real wv = -Real(1)/9;
  static const Real wm =  Real(4)/9;

  for (Elem * elem : _mesh.element_ptr_range())
    {
      libmesh_assert_equal_to(elem->n_vertices(), 3);
      libmesh_assert_equal_to(elem->n_nodes(), 7u);

      elem->point(6) = wv * (elem->point(0) +
                             elem->point(1) +
                             elem->point(2)) +
                       wm * (elem->point(3) +
                             elem->point(4) +
                             elem->point(5));
    }
}


void TriangulatorInterface::verify_quadratic_elements()
{
  if (_elem_type != TRI6 && _elem_type != TRI7)
    return;

  // Once fixup_tri7_center_nodes() has placed node 6, the TRI6 and TRI7
  // mappings coincide and this Tri6 formula serves both.
  static const Real xi_samples[7]  = {Real(0),   Real(1),   Real(0),
                                      Real(1)/2, Real(1)/2, Real(0),
                                      Real(1)/3};
  static const Real eta_samples[7] = {Real(0),   Real(0),   Real(1),
                                      Real(0),   Real(1)/2, Real(1)/2,
                                      Real(1)/3};

  for (Elem * elem : _mesh.element_ptr_range())
    {
      libmesh_assert_equal_to(elem->n_vertices(), 3);
      libmesh_assert_greater_equal(elem->n_nodes(), 6u);

      const Point & x0 = elem->point(0);
      const Point & x1 = elem->point(1);
      const Point & x2 = elem->point(2);
      const Point & x3 = elem->point(3);
      const Point & x4 = elem->point(4);
      const Point & x5 = elem->point(5);

      // Tri6 mapping derivative coefficients (see Tri6::volume()):
      // dx/dxi = xi*a1 + eta*b1 + c1, dx/deta = xi*b1 + eta*b2 + c2.
      const Point a1 =  4*x0 + 4*x1 - 8*x3;
      const Point b1 =  4*x0 - 4*x3 + 4*x4 - 4*x5;
      const Point c1 = -3*x0 - 1*x1 + 4*x3;
      const Point b2 =  4*x0 + 4*x2 - 8*x5;
      const Point c2 = -3*x0 - 1*x2 + 4*x5;

      // Scale the tolerance by the straight-edge triangle area, which
      // is strictly positive for the valid TRI3 poly2tri input.
      const Real ref_area = 0.5 * cross_norm(x1 - x0, x2 - x0);
      const Real jac_tol = TOLERANCE * ref_area;

      Real min_jac = std::numeric_limits<Real>::max();
      unsigned int worst_sample = 0;
      for (unsigned int s = 0; s != 7; ++s)
        {
          const Real xi  = xi_samples[s];
          const Real eta = eta_samples[s];
          const Point dxi  = xi*a1 + eta*b1 + c1;
          const Point deta = xi*b1 + eta*b2 + c2;
          // z-component of the cross product; the elements are planar.
          const Real jac = dxi(0)*deta(1) - dxi(1)*deta(0);
          if (jac < min_jac)
            {
              min_jac = jac;
              worst_sample = s;
            }
        }

      if (min_jac > jac_tol)
        continue;

      // Build a diagnostic naming every snapped boundary side on this
      // element so the user can immediately see which curved-boundary
      // input caused the tangle.
      std::ostringstream sides;
      for (unsigned int n = 0; n != 3; ++n)
        if (!elem->neighbor_ptr(n))
          {
            const Point straight =
              0.5 * (elem->point(n) + elem->point((n+1) % 3));
            sides << " (boundary side " << n
                  << ": straight midpoint " << straight
                  << ", snapped midpoint " << elem->point(n+3) << ")";
          }

      libmesh_error_msg(
        "TriangulatorInterface: snapping a boundary midpoint produced a "
        "tangled quadratic triangle (element " << elem->id()
        << ", non-positive Jacobian " << min_jac
        << " at reference sample (" << xi_samples[worst_sample] << ", "
        << eta_samples[worst_sample] << "); reference triangle area "
        << ref_area << ")." << sides.str()
        << " Refine the boundary discretization so that recorded "
        "midpoints lie closer to their straight-line midpoints, "
        "then retry.");
    }
}


void TriangulatorInterface::verify_holes(const Hole & outer_bdy)
{
  for (const Hole * hole : *_holes)
    {
      for (const Hole * hole2 : *_holes)
        {
          if (hole == hole2)
            continue;

          for (auto i : make_range(hole2->n_points()))
            if (hole->contains(hole2->point(i)))
              libmesh_error_msg
                ("Found point " << hole2->point(i) <<
                 " on one hole boundary and another's interior");
        }

      for (auto i : make_range(hole->n_points()))
        if (!outer_bdy.contains(hole->point(i)))
          libmesh_error_msg
            ("Found point " << hole->point(i) <<
             " on hole boundary but outside outer boundary");
    }
}


unsigned int TriangulatorInterface::total_hole_points()
{
  // If the holes vector is non-nullptr (and non-empty) we need to determine
  // the number of additional points which the holes will add to the
  // triangulation.
  // Note that the number of points is always equal to the number of segments
  // that form the holes.
  unsigned int n_hole_points = 0;

  if (_holes)
    for (const auto & hole : *_holes)
    {
      n_hole_points += hole->n_points();
      // A hole at least has one enclosure.
      // Points on enclosures are ordered so that we can add segments implicitly.
      // Elements in segment_indices() indicates the starting points of all enclosures.
      // The last element in segment_indices() is the number of total points.
      libmesh_assert_greater(hole->segment_indices().size(), 1);
      libmesh_assert_equal_to(hole->segment_indices().back(), hole->n_points());
    }

  return n_hole_points;
}

void TriangulatorInterface::set_auto_area_function(const Parallel::Communicator &comm,
                                                   const unsigned int num_nearest_pts,
                                                   const unsigned int power,
                                                   const Real background_value,
                                                   const Real  background_eff_dist)
{
   _auto_area_function = std::make_unique<AutoAreaFunction>(comm, num_nearest_pts, power, background_value, background_eff_dist);
}

FunctionBase<Real> * TriangulatorInterface::get_auto_area_function()
{
  if (!_auto_area_function->initialized())
  {
    // Points and target element sizes for the interpolation
    std::vector<Point> function_points;
    std::vector<Real> function_sizes;
    calculate_auto_desired_area_samples(function_points, function_sizes);
    _auto_area_function->init_mfi(function_points, function_sizes);
  }
  return _auto_area_function.get();
}

void TriangulatorInterface::calculate_auto_desired_area_samples(std::vector<Point> & function_points,
                                                                std::vector<Real> & function_sizes,
                                                                const Real & area_factor)
{
  // Get the hole mesh of the outer boundary
  // Holes should already be attached if applicable when this function is called
  const TriangulatorInterface::MeshedHole bdry_mh { _mesh, this->_bdy_ids };
  // Collect all the centroid points of the outer boundary segments
  // and the corresponding element sizes
  for (unsigned int i = 0; i < bdry_mh.n_points(); i++)
  {
    function_points.push_back((bdry_mh.point(i) + bdry_mh.point((i + 1) % bdry_mh.n_points())) /
                              Real(2.0));
    function_sizes.push_back(
        (bdry_mh.point(i) - bdry_mh.point((i + 1) % bdry_mh.n_points())).norm());
  }
  // If holes are present, do the same for the hole boundaries
  if(_holes)
    for (const Hole * hole : *_holes)
    {
      for (unsigned int i = 0; i < hole->n_points(); i++)
      {
        function_points.push_back(
            (hole->point(i) + hole->point((i + 1) % hole->n_points())) / Real(2.0));
        function_sizes.push_back(
            (hole->point(i) - hole->point((i + 1) % hole->n_points())).norm());
      }
    }

  std::for_each(
      function_sizes.begin(), function_sizes.end(), [&area_factor](Real & a) { a = a * a * area_factor * std::sqrt(3.0) / 4.0; });

}
} // namespace libMesh

