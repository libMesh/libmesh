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


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_POLY2TRI

// libmesh includes
#include "libmesh/poly2tri_triangulator.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/function_base.h"
#include "libmesh/hashing.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/mesh_triangle_holes.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h"

// poly2tri includes
#include "libmesh/ignore_warnings.h" // utf-8 comments should be fine...
#include "poly2tri/poly2tri.h"
#include "libmesh/restore_warnings.h"

// Anonymous namespace - poly2tri doesn't define operator<(Point,Point)
namespace
{
using namespace libMesh;

struct P2TPointCompare
{
  bool operator()(const p2t::Point & a, const p2t::Point & b) const
  {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  }
};

p2t::Point to_p2t(const libMesh::Point & p)
{
#if LIBMESH_DIM > 2
  libmesh_error_msg_if
    (p(2) != 0,
     "Poly2TriTriangulator only supports point sets in the XY plane");
#endif

  return {p(0), p(1)};
}

Real distance_from_circumcircle(const Elem & elem,
                                const Point & p)
{
  libmesh_assert_equal_to(elem.n_vertices(), 3);

  const Point circumcenter = elem.quasicircumcenter();
  const Real radius = (elem.point(0) - circumcenter).norm();
  const Real p_dist = (p - circumcenter).norm();

  return p_dist - radius;
}


bool in_circumcircle(const Elem & elem,
                     const Point & p,
                     const Real tol = 0)
{
  return (distance_from_circumcircle(elem, p) < tol);

  /*
  libmesh_assert_equal_to(elem.n_vertices(), 3);

  const Point pv0 = elem.point(0) - p;
  const Point pv1 = elem.point(1) - p;
  const Point pv2 = elem.point(2) - p;

  return ((pv0.norm_sq() * (pv1(0)*pv2(1)-pv2(0)*pv1(1))) -
          (pv1.norm_sq() * (pv0(0)*pv2(1)-pv2(0)*pv0(1))) +
          (pv2.norm_sq() * (pv0(0)*pv1(1)-pv1(0)*pv0(1)))) > 0;
  */
}


std::pair<bool, unsigned short>
can_delaunay_swap(const Elem & elem,
                  unsigned short side,
                  Real tol)
{
  const Elem * neigh = elem.neighbor_ptr(side);
  if (!neigh)
    return {false, 0};

  unsigned short nn = 0;

  // What neighbor node does elem not share?
  for (; nn < 3; ++nn)
    {
      const Node * neigh_node = neigh->node_ptr(nn);
      if (neigh_node == elem.node_ptr(0) ||
          neigh_node == elem.node_ptr(1) ||
          neigh_node == elem.node_ptr(2))
        continue;

      // Might we need to do a diagonal swap here?  Avoid
      // undoing a borderline swap.
      if (in_circumcircle(elem, *neigh_node, tol))
        break;
    }

  if (nn == 3)
    return {false, 0};

  const unsigned short n = (side+2)%3;
  const RealVectorValue right =
    (elem.point((n+1)%3)-elem.point(n)).unit();
  const RealVectorValue mid =
    (neigh->point(nn)-elem.point(n)).unit();
  const RealVectorValue left =
    (elem.point((n+2)%3)-elem.point(n)).unit();

  // If the "middle" vector isn't really in the middle, we can't do a
  // swap without involving other triangles (or we can't at all if
  // there's a domain boundary in the way)
  if (mid*right < left*right ||
      left*mid < left*right)
    return {false, 0};

  return {true, nn};
}


[[maybe_unused]] void libmesh_assert_locally_delaunay(const Elem & elem)
{
  libmesh_ignore(elem);

#ifndef NDEBUG
  // -TOLERANCE, because we're fine with something a little inside the
  // circumcircle
  for (auto s : make_range(elem.n_sides()))
    libmesh_assert(!can_delaunay_swap(elem, s, -TOLERANCE).first);
#endif
}

template <typename Container>
inline
void libmesh_assert_delaunay(MeshBase & libmesh_dbg_var(mesh),
                             Container & new_elems)
{
  libmesh_ignore(new_elems);
#ifndef NDEBUG
  LOG_SCOPE("libmesh_assert_delaunay()", "Poly2TriTriangulator");

  for (auto & elem : mesh.element_ptr_range())
    libmesh_assert_locally_delaunay(*elem);

  for (auto & [raw_elem, unique_elem] : new_elems)
    {
      libmesh_ignore(unique_elem); // avoid warnings on old gcc
      libmesh_assert_locally_delaunay(*raw_elem);
    }
#endif
}


// Restore a triangulation's Delaunay property, starting with a set of
// all triangles that might initially not be locally Delaunay with
// their neighbors.
template <typename Container>
inline
void restore_delaunay(Container & check_delaunay_on,
                      BoundaryInfo & boundary_info)
{
  LOG_SCOPE("restore_delaunay()", "Poly2TriTriangulator");

  while (!check_delaunay_on.empty())
    {
      Elem & elem = **check_delaunay_on.begin();
      check_delaunay_on.erase(&elem);
      for (auto s : make_range(elem.n_sides()))
        {
          // Can we make a swap here?  With what neighbor, with what
          // far node?  Use a negative tolerance to avoid swapping
          // back and forth.
          auto [can_swap, nn] =
            can_delaunay_swap(elem, s, -TOLERANCE*TOLERANCE);
          if (!can_swap)
            continue;

          Elem * neigh = elem.neighbor_ptr(s);

          // If we made it here it's time to diagonal swap
          const unsigned short n = (s+2)%3;

          const std::array<Node *,4> nodes {elem.node_ptr(n),
            elem.node_ptr((n+1)%3), neigh->node_ptr(nn),
            elem.node_ptr((n+2)%3)};

          // If we have to swap then either we or any of our neighbors
          // might no longer be Delaunay
          for (auto ds : make_range(3))
            {
              if (elem.neighbor_ptr(ds))
                check_delaunay_on.insert(elem.neighbor_ptr(ds));
              if (neigh->neighbor_ptr(ds))
                check_delaunay_on.insert(neigh->neighbor_ptr(ds));
            }

          // An interior boundary between two newly triangulated
          // triangles shouldn't have any bcids
          libmesh_assert(!boundary_info.n_boundary_ids(neigh, (nn+1)%3));
          libmesh_assert(!boundary_info.n_boundary_ids(&elem, (n+1)%3));

          // The two changing boundary sides might have bcids
          std::vector<boundary_id_type> bcids;
          boundary_info.boundary_ids(&elem, (n+2)%3, bcids);
          if (!bcids.empty())
            {
              boundary_info.remove_side(&elem, (n+2)%3);
              boundary_info.add_side(neigh, (nn+1)%3, bcids);
            }

          boundary_info.boundary_ids(neigh, (nn+2)%3, bcids);
          if (!bcids.empty())
            {
              boundary_info.remove_side(neigh, (nn+2)%3);
              boundary_info.add_side(&elem, (n+1)%3, bcids);
            }

          elem.set_node((n+2)%3) = nodes[2];
          neigh->set_node((nn+2)%3) = nodes[0];

          // No need for a temporary array to swap these around, if we
          // do it in the right order.
          //
          // Watch me neigh->neigh
          Elem * neighneigh = neigh->neighbor_ptr((nn+2)%3);
          if (neighneigh)
            {
              unsigned int snn = neighneigh->which_neighbor_am_i(neigh);
              neighneigh->set_neighbor(snn, &elem);
            }

          Elem * elemoldneigh = elem.neighbor_ptr((n+2)%3);
          if (elemoldneigh)
            {
              unsigned int seon = elemoldneigh->which_neighbor_am_i(&elem);
              elemoldneigh->set_neighbor(seon, neigh);
            }

          elem.set_neighbor((n+1)%3, neigh->neighbor_ptr((nn+2)%3));
          neigh->set_neighbor((nn+1)%3, elem.neighbor_ptr((n+2)%3));
          elem.set_neighbor((n+2)%3, neigh);
          neigh->set_neighbor((nn+2)%3, &elem);

          // Start over after this much change, don't just loop to the
          // next neighbor
          break;
        }
    }
}


unsigned int segment_intersection(const Elem & elem,
                                  Point & source,
                                  const Point & target,
                                  unsigned int source_side)
{
  libmesh_assert_equal_to(elem.dim(), 2);

  const auto ns = elem.n_sides();

  for (auto s : make_range(ns))
    {
      // Don't go backwards just because some FP roundoff said to
      if (s == source_side)
        continue;

      const Point v0 = elem.point(s);
      const Point v1 = elem.point((s+1)%ns);

      // Calculate intersection parameters (fractions of the distance
      // along each segment)
      const Real raydx = target(0)-source(0),
                 raydy = target(1)-source(1),
                 edgedx = v1(0)-v0(0),
                 edgedy = v1(1)-v0(1);
      const Real denom = edgedx * raydy - edgedy * raydx;

      // divide-by-zero means the segments are parallel
      if (denom == 0)
        continue;

      const Real one_over_denom = 1 / denom;

      const Real targetsdx = v1(0)-target(0),
                 targetsdy = v1(1)-target(1);

      const Real t_num = targetsdx * raydy -
                         targetsdy * raydx;
      const Real t = t_num * one_over_denom;

      if (t < -TOLERANCE*TOLERANCE || t > 1 + TOLERANCE*TOLERANCE)
        continue;

      const Real u_num = targetsdx * edgedy - targetsdy * edgedx;
      const Real u = u_num * one_over_denom;

      if (u < -TOLERANCE*TOLERANCE || u > 1 + TOLERANCE*TOLERANCE)
        continue;

/*
      // Partial workaround for an old poly2tri bug (issue #39): if we
      // end up with boundary points that are nearly-collinear but
      // infinitesimally concave, p2t::CDT::Triangulate throws a "null
      // triangle" exception.  So let's try to be infinitesimally
      // convex instead.
      const Real ray_fraction = (1-u) * (1+TOLERANCE*TOLERANCE);
*/
      const Real ray_fraction = (1-u);

      source(0) += raydx * ray_fraction;
      source(1) += raydy * ray_fraction;
      return s;
    }

  return libMesh::invalid_uint;
}

}

namespace libMesh
{

//
// Function definitions for the Poly2TriTriangulator class
//

// Constructor
Poly2TriTriangulator::Poly2TriTriangulator(UnstructuredMesh & mesh,
                                           dof_id_type n_boundary_nodes)
  : TriangulatorInterface(mesh),
    _n_boundary_nodes(n_boundary_nodes),
    _refine_bdy_allowed(true)
{
}


Poly2TriTriangulator::~Poly2TriTriangulator() = default;


// Primary function responsible for performing the triangulation
void Poly2TriTriangulator::triangulate()
{
  LOG_SCOPE("triangulate()", "Poly2TriTriangulator");

  // We only operate on serialized meshes.  And it's not safe to
  // serialize earlier, because it would then be possible for the user
  // to re-parallelize the mesh in between there and here.
  MeshSerializer serializer(_mesh);

  // We don't yet support every set of Triangulator options in the
  // poly2tri implementation

  // We don't support convex hull triangulation, only triangulation of
  // (implicitly defined, by node ordering) polygons (with holes if
  // requested)
  if (_triangulation_type != PSLG)
    libmesh_not_implemented();

  // We currently don't handle region specifications
  if (_regions)
    libmesh_not_implemented();

  // We won't support quads any time soon, or 1D/3D in this interface
  // ever.
  if (_elem_type != TRI3 &&
      _elem_type != TRI6 &&
      _elem_type != TRI7)
    libmesh_not_implemented();

  // If we have no explicit segments defined, we may get them from
  // mesh elements
  this->elems_to_segments();

  // If we *still* have no explicit segments defined, we get them from
  // the order of nodes.
  this->nodes_to_segments(_n_boundary_nodes);

  // Insert additional new points in between existing boundary points,
  // if that is requested and reasonable
  this->insert_any_extra_boundary_points();

  // Triangulate the points we have, then see if we need to add more;
  // repeat until we don't need to add more.
  //
  // This is currently done redundantly in parallel; make sure no
  // processor quits before the others.
  do
    {
      libmesh_parallel_only(_mesh.comm());
      this->triangulate_current_points();
    }
  while (this->insert_refinement_points());

  libmesh_parallel_only(_mesh.comm());

  // Okay, we really do need to support boundary ids soon, but we
  // don't yet
  if (_markers)
    libmesh_not_implemented();

  _mesh.set_mesh_dimension(2);

  // To the naked eye, a few smoothing iterations usually looks better,
  // so we do this by default unless the user says not to.
  if (this->_smooth_after_generating)
    LaplaceMeshSmoother(_mesh).smooth(2);

  // Prepare the mesh for use before returning.  This ensures (among
  // other things) that it is partitioned and therefore users can
  // iterate over local elements, etc.
  _mesh.prepare_for_use();
}


void Poly2TriTriangulator::set_desired_area_function
  (FunctionBase<Real> * desired)
{
  if (desired)
    _desired_area_func = desired->clone();
  else
    _desired_area_func.reset();
}


FunctionBase<Real> * Poly2TriTriangulator::get_desired_area_function ()
{
  return _desired_area_func.get();
}


bool Poly2TriTriangulator::is_refine_boundary_allowed
  (const BoundaryInfo & boundary_info,
   const Elem & elem,
   unsigned int side)
{
  // We should only be calling this on a boundary side
  libmesh_assert(!elem.neighbor_ptr(side));

  std::vector<boundary_id_type> bcids;
  boundary_info.boundary_ids(&elem, side, bcids);

  // We should have one bcid on every boundary side.
  libmesh_assert_equal_to(bcids.size(), 1);

  if (bcids[0] == 0)
    return this->refine_boundary_allowed();

  // If we're not on an outer boundary side we'd better be on a hole
  // side
  libmesh_assert(this->_holes);

  const boundary_id_type hole_num = bcids[0]-1;
  libmesh_assert_less(hole_num, this->_holes->size());
  const Hole * hole = (*this->_holes)[hole_num];
  return hole->refine_boundary_allowed();
}


void Poly2TriTriangulator::triangulate_current_points()
{
  LOG_SCOPE("triangulate_current_points()", "Poly2TriTriangulator");

  // Will the triangulation have holes?
  const std::size_t n_holes = _holes != nullptr ? _holes->size() : 0;

  // Mapping from Poly2Tri points to libMesh nodes, so we can get the
  // connectivity translated back later.
  std::map<const p2t::Point, Node *, P2TPointCompare> point_node_map;

  // Poly2Tri data structures
  // Poly2Tri takes vectors of pointers-to-Point for some reason, but
  // we'll just make those shims to vectors of Point rather than
  // individually/manually heap allocating everything.
  std::vector<p2t::Point> outer_boundary_points;
  std::vector<std::vector<p2t::Point>> inner_hole_points(n_holes);

  dof_id_type nn = _mesh.max_node_id();
  libmesh_error_msg_if
    (!nn, "Poly2TriTriangulator cannot triangulate an empty mesh!");

  // Unless we're using an explicit segments list, we assume node ids
  // are contiguous here.
  if (this->segments.empty())
    libmesh_error_msg_if
      (_mesh.n_nodes() != nn,
       "Poly2TriTriangulator needs contiguous node ids or explicit segments!");

  // And if we have more nodes than outer boundary points, the rest
  // may be interior "Steiner points".  We use a set here so we can
  // cheaply search and erase from it later, when we're identifying
  // hole points.
  std::set<p2t::Point, P2TPointCompare> steiner_points;

  // If we were asked to use all mesh nodes as boundary nodes, now's
  // the time to see how many that is.
  if (_n_boundary_nodes == DofObject::invalid_id)
    {
      _n_boundary_nodes = _mesh.n_nodes();
      libmesh_assert_equal_to(std::ptrdiff_t(_n_boundary_nodes),
                              std::distance(_mesh.nodes_begin(),
                                            _mesh.nodes_end()));

    }
  else
    libmesh_assert_less_equal(_n_boundary_nodes,
                              _mesh.n_nodes());

  // Prepare poly2tri points for our nodes, sorted into outer boundary
  // points and interior Steiner points.

  if (this->segments.empty())
    {
      // If we have no segments even after taking elems into account,
      // the nodal id ordering defines our outer polyline ordering
      for (auto & node : _mesh.node_ptr_range())
        {
          const p2t::Point pt = to_p2t(*node);

          // If we're out of boundary nodes, the rest are going to be
          // Steiner points or hole points
          if (node->id() < _n_boundary_nodes)
            outer_boundary_points.push_back(pt);
          else
            steiner_points.insert(pt);

          // We're not going to support overlapping nodes on the boundary
          if (point_node_map.count(pt))
            libmesh_not_implemented();

          point_node_map.emplace(pt, node);
        }
    }
  // If we have explicit segments defined, that's our outer polyline
  // ordering:
  else
    {
      // Let's make sure our segments are in order
      dof_id_type last_id = DofObject::invalid_id;

      // Add nodes from every segment, in order, to the outer polyline
      for (auto [segment_start, segment_end] : this->segments)
        {
          if (last_id != DofObject::invalid_id)
            libmesh_error_msg_if(segment_start != last_id,
                                 "Disconnected triangulator segments");
          last_id = segment_end;

          Node * node = _mesh.query_node_ptr(segment_start);

          libmesh_error_msg_if(!node,
                               "Triangulator segments reference nonexistent node id " <<
                               segment_start);

          outer_boundary_points.emplace_back((*node)(0), (*node)(1));
          p2t::Point * pt = &outer_boundary_points.back();

          // We're not going to support overlapping nodes on the boundary
          if (point_node_map.count(*pt))
            libmesh_not_implemented_msg
              ("Triangulating overlapping boundary nodes is unsupported");

          point_node_map.emplace(*pt, node);
        }

        libmesh_error_msg_if(last_id != this->segments[0].first,
                             "Non-closed-loop triangulator segments");

      // If we have points that aren't in any segments, those will be
      // Steiner points
      for (auto & node : _mesh.node_ptr_range())
        {
          const p2t::Point pt = to_p2t(*node);

          auto it = point_node_map.find(pt);

          if (it == point_node_map.end())
            {
              steiner_points.insert(pt);
              point_node_map.emplace(pt, node);
            }
          else
            libmesh_assert_equal_to(it->second, node);
        }
    }

  // If we have any elements from a previous triangulation, we're
  // going to replace them with our own.  If we have any elements that
  // were used to create our segments, we're done creating and we no
  // longer need them.
  _mesh.clear_elems();

  // Keep track of what boundary ids we want to assign to each new
  // triangle.  We'll give the outer boundary BC 0, and give holes ids
  // starting from 1.  We've already got the point_node_map to find
  // nodes, so we can just key on pairs of node ids to identify a side.
  std::unordered_map<std::pair<dof_id_type,dof_id_type>,
                     boundary_id_type, libMesh::hash> side_boundary_id;

  const boundary_id_type outer_bcid = 0;
  const std::size_t n_outer = outer_boundary_points.size();

  for (auto i : make_range(n_outer))
    {
      const Node * node1 =
        libmesh_map_find(point_node_map, outer_boundary_points[i]),
                 * node2 =
        libmesh_map_find(point_node_map, outer_boundary_points[(i+1)%n_outer]);

      side_boundary_id.emplace(std::make_pair(node1->id(),
                                              node2->id()),
                               outer_bcid);
    }

  // Create poly2tri triangulator with our mesh points
  std::vector<p2t::Point *> outer_boundary_pointers(n_outer);
  std::transform(outer_boundary_points.begin(),
                 outer_boundary_points.end(),
                 outer_boundary_pointers.begin(),
                 [](p2t::Point & p) { return &p; });


  // Make sure shims for holes last as long as the CDT does; the
  // poly2tri headers don't make clear whether or not they're hanging
  // on to references to these vectors, and it would be reasonable for
  // them to do so.
  std::vector<std::vector<p2t::Point *>> inner_hole_pointers(n_holes);

  p2t::CDT cdt{outer_boundary_pointers};

  // Add any holes
  for (auto h : make_range(n_holes))
    {
      const Hole * initial_hole = (*_holes)[h];
      auto it = replaced_holes.find(initial_hole);
      const Hole & our_hole =
        (it == replaced_holes.end()) ?
        *initial_hole : *it->second;
      auto & poly2tri_hole = inner_hole_points[h];

      for (auto i : make_range(our_hole.n_points()))
        {
          Point p = our_hole.point(i);
          poly2tri_hole.emplace_back(to_p2t(p));

          const auto & pt = poly2tri_hole.back();

          // This won't be a steiner point.
          steiner_points.erase(pt);

          // If we see a hole point already in the mesh, we'll share
          // that node.  This might be a problem if it's a boundary
          // node, but it might just be the same hole point already
          // added during a previous triangulation refinement step.
          if (point_node_map.count(pt))
            {
              libmesh_assert_equal_to
                (point_node_map[pt],
                 _mesh.query_node_ptr(point_node_map[pt]->id()));
            }
          else
            {
              Node * node = _mesh.add_point(p, nn++);
              point_node_map[pt] = node;
            }
        }

      const boundary_id_type inner_bcid = h+1;
      const std::size_t n_inner = poly2tri_hole.size();

      for (auto i : make_range(n_inner))
        {
          const Node * node1 =
            libmesh_map_find(point_node_map, poly2tri_hole[i]),
                     * node2 =
            libmesh_map_find(point_node_map, poly2tri_hole[(i+1)%n_inner]);

          side_boundary_id.emplace(std::make_pair(node1->id(),
                                                  node2->id()),
                                   inner_bcid);
        }

      auto & poly2tri_ptrs = inner_hole_pointers[h];
      poly2tri_ptrs.resize(n_inner);

      std::transform(poly2tri_hole.begin(),
                     poly2tri_hole.end(),
                     poly2tri_ptrs.begin(),
                     [](p2t::Point & p) { return &p; });

      cdt.AddHole(poly2tri_ptrs);
    }

  // Add any steiner points.  We had them in a set, but post-C++11
  // that won't give us non-const element access (even if we
  // pinky-promise not to change the elements in any way that affects
  // our Comparator), and Poly2Tri wants non-const elements (to store
  // edge data?), so we need to move them here.
  std::vector<p2t::Point> steiner_vector(steiner_points.begin(), steiner_points.end());
  steiner_points.clear();
  for (auto & p : steiner_vector)
    cdt.AddPoint(&p);

  // Triangulate!
  cdt.Triangulate();

  // Get poly2tri triangles, turn them into libMesh triangles
  std::vector<p2t::Triangle *> triangles = cdt.GetTriangles();

  // Do our own numbering, even on DistributedMesh
  dof_id_type next_id = 0;

  BoundaryInfo & boundary_info = _mesh.get_boundary_info();
  boundary_info.clear();

  // Add the triangles to our Mesh data structure.
  for (auto ptri_ptr : triangles)
    {
      p2t::Triangle & ptri = *ptri_ptr;

      // We always use TRI3 here, since that's what we have nodes for;
      // if we need a higher order we can convert at the end.
      auto elem = Elem::build_with_id(TRI3, next_id++);
      for (auto v : make_range(3))
        {
          const p2t::Point & vertex = *ptri.GetPoint(v);

          Node * node = libmesh_map_find(point_node_map, vertex);
          libmesh_assert(node);
          elem->set_node(v) = node;
        }

      // We expect a consistent triangle orientation
      libmesh_assert(!elem->is_flipped());

      Elem * added_elem = _mesh.add_elem(std::move(elem));

      for (auto v : make_range(3))
        {
          const Node & node1 = added_elem->node_ref(v),
                     & node2 = added_elem->node_ref((v+1)%3);

          auto it = side_boundary_id.find(std::make_pair(node1.id(), node2.id()));
          if (it == side_boundary_id.end())
            it = side_boundary_id.find(std::make_pair(node2.id(), node1.id()));
          if (it != side_boundary_id.end())
            boundary_info.add_side(added_elem, v, it->second);
        }
    }
}



bool Poly2TriTriangulator::insert_refinement_points()
{
  LOG_SCOPE("insert_refinement_points()", "Poly2TriTriangulator");

  if (this->minimum_angle() != 0)
    libmesh_not_implemented();

  // We need neighbor pointers for ray casting and cavity finding
  UnstructuredMesh & mesh = dynamic_cast<UnstructuredMesh &>(this->_mesh);
  mesh.find_neighbors();

  if (this->desired_area() == 0 &&
      this->get_desired_area_function() == nullptr &&
      !this->has_auto_area_function())
    return false;

  BoundaryInfo & boundary_info = _mesh.get_boundary_info();

  // We won't immediately add these, lest we invalidate iterators on a
  // ReplicatedMesh.  They'll still be in the mesh neighbor topology
  // for the purpose of doing Delaunay cavity stuff, so we need to
  // manage memory here, but there's no point in adding them to the
  // Mesh just to remove them again afterward when we hit up poly2tri.

  // We'll need to be able to remove new elems from new_elems, in
  // cases where a later refinement insertion has a not-yet-added
  // element in its cavity, so we'll use a map here to make searching
  // possible.
  //
  // For parallel consistency, we can't order a container we plan to
  // iterate through based on Elem * or a hash of it.  We'll be doing
  // Delaunay swaps so we can't iterate based on geometry.  These are
  // not-yet-added elements so we can't iterate based on proper
  // element ids ... but we can set a temporary element id to use for
  // the purpose.
  struct cmp {
    bool operator()(Elem * a, Elem * b) const {
      libmesh_assert(a == b || a->id() != b->id());
      return (a->id() < b->id());
    }
  } comp;

  std::map<Elem *, std::unique_ptr<Elem>, decltype(comp)> new_elems(comp);

  // We should already be Delaunay when we get here, otherwise we
  // won't be able to stay Delaunay later.  But we're *not* always
  // Delaunay when we get here?  What the hell, poly2tri?  Fixing this
  // is expensive!
  {
    // restore_delaunay should get to the same Delaunay triangulation up to
    // isomorphism regardless of ordering ... but we actually care
    // about the isomorphisms!  If a triangle's nodes are permuted on
    // one processor vs another that's an issue.  So sort our input
    // carefully.
    std::set<Elem *, decltype(comp)> all_elems
      { mesh.elements_begin(), mesh.elements_end(), comp };

    restore_delaunay(all_elems, boundary_info);

    libmesh_assert_delaunay(mesh, new_elems);
  }

  // Map of which points follow which in the boundary polylines.  If
  // we have to add new boundary points, we'll use this to construct
  // an updated this->segments to retriangulate with.  If we have to
  // add new hole points, we'll use this to insert points into an
  // ArbitraryHole.
  std::unordered_map<Point, Node *> next_boundary_node;

  // In cases where we've been working with contiguous node id ranges;
  // let's keep it that way.
  dof_id_type nn = _mesh.max_node_id();
  dof_id_type ne = _mesh.max_elem_id();

  // We can't handle duplicated nodes.  We shouldn't ever create one,
  // but let's make sure of that.
#ifdef DEBUG
  std::unordered_set<Point> mesh_points;
  for (const Node * node : mesh.node_ptr_range())
    {
      libmesh_assert(!mesh_points.count(*node));
      mesh_points.insert(*node);
    }
#endif

  auto add_point = [&mesh,
#ifdef DEBUG
                    &mesh_points,
#endif
                    &nn](const Point & p)
    {
#ifdef DEBUG
      libmesh_assert(!mesh_points.count(p));
      mesh_points.insert(p);
#endif
      return mesh.add_point(p, nn++);
    };

  for (auto & elem : mesh.element_ptr_range())
    {
      // element_ptr_range skips deleted elements ... right?
      libmesh_assert(elem);
      libmesh_assert(elem->valid_id());

      // We only handle triangles in our triangulation
      libmesh_assert_equal_to(elem->level(), 0u);
      libmesh_assert_equal_to(elem->type(), TRI3);

      // If this triangle is as small as we desire, move along
      if (!should_refine_elem(*elem))
        continue;

      // Otherwise add a Steiner point.  We'd like to add the
      // circumcenter ...
      Point new_pt = elem->quasicircumcenter();

      // And to give it a node;
      Node * new_node = nullptr;

      // But that might be outside our triangle, or even outside the
      // boundary.  We'll find a triangle that should contain our new
      // point
      Elem * cavity_elem = elem; // Start looking at elem anyway

      // We'll refine a boundary later if necessary.
      auto boundary_refine = [this, &next_boundary_node,
                              &cavity_elem, &new_node]
        (unsigned int side)
      {
        libmesh_ignore(this); // Only used in dbg/devel
        libmesh_assert(new_node);
        libmesh_assert(new_node->valid_id());

        Node * old_segment_start = cavity_elem->node_ptr(side),
             * old_segment_end   = cavity_elem->node_ptr((side+1)%3);
        libmesh_assert(old_segment_start);
        libmesh_assert(old_segment_start->valid_id());
        libmesh_assert(old_segment_end);
        libmesh_assert(old_segment_end->valid_id());

        auto it = next_boundary_node.find(*old_segment_start);
        if (it != next_boundary_node.end())
          {
            libmesh_assert(it->second == old_segment_end);
            it->second = new_node;
          }
        else
          {
            // This would be an O(N) sanity check if we already
            // have a segments vector or any holes.  :-P
            libmesh_assert(!this->segments.empty() ||
                           (_holes && !_holes->empty()) ||
                           (old_segment_end->id() ==
                            old_segment_start->id() + 1));
            next_boundary_node[*old_segment_start] = new_node;
          }

        next_boundary_node[*new_node] = old_segment_end;
      };

      // Let's find a triangle containing our new point, or at least
      // containing the end of a ray leading from our current triangle
      // to the new point.
      Point ray_start = elem->vertex_average();

      // What side are we coming from, and what side are we going to?
      unsigned int source_side = invalid_uint;
      unsigned int side = invalid_uint;

      while (!cavity_elem->contains_point(new_pt))
        {
          side = segment_intersection(*cavity_elem, ray_start, new_pt, source_side);

          libmesh_assert_not_equal_to (side, invalid_uint);

          Elem * neigh = cavity_elem->neighbor_ptr(side);
          // If we're on a boundary, stop there.  Refine the boundary
          // if we're allowed, the boundary element otherwise.
          if (!neigh)
            {
              if (this->is_refine_boundary_allowed(boundary_info,
                                                   *cavity_elem,
                                                   side))
                {
                  new_pt = ray_start;
                  new_node = add_point(new_pt);
                  boundary_refine(side);
                }
              else
                {
                  // Should we just add the vertex average of the
                  // boundary element, to minimize the number of
                  // slivers created?
                  //
                  // new_pt = cavity_elem->vertex_average();
                  //
                  // That works for a while, but it
                  // seems to be able to "run away" and leave us with
                  // crazy slivers on boundaries if we push interior
                  // refinement too far while disabling boundary
                  // refinement.
                  //
                  // Let's go back to refining our original problem
                  // element.
                  cavity_elem = elem;
                  new_pt = cavity_elem->vertex_average();
                  new_node = add_point(new_pt);

                  // This was going to be a side refinement but it's
                  // now an internal refinement
                  side = invalid_uint;
                }

              break;
            }

          source_side = neigh->which_neighbor_am_i(cavity_elem);
          cavity_elem = neigh;
          side = invalid_uint;
        }

      // If we're ready to create a new node and we're not on a
      // boundary ... should we be?  We don't want to create any
      // sliver elements or confuse poly2tri or anything.
      if (side == invalid_uint && !new_node)
        {
          unsigned int worst_side = libMesh::invalid_uint;
          Real worst_cos = 1;
          for (auto s : make_range(3u))
            {
              // We never snap to a non-domain-boundary
              if (cavity_elem->neighbor_ptr(s))
                continue;

              Real ax = cavity_elem->point(s)(0) - new_pt(0),
                   ay = cavity_elem->point(s)(1) - new_pt(1),
                   bx = cavity_elem->point((s+1)%3)(0) - new_pt(0),
                   by = cavity_elem->point((s+1)%3)(1) - new_pt(1);
              const Real my_cos = (ax*bx+ay*by) /
                                  std::sqrt(ax*ax+ay*ay) /
                                  std::sqrt(bx*bx+by*by);

              if (my_cos < worst_cos)
                {
                  worst_side = s;
                  worst_cos = my_cos;
                }
            }

          // If we'd create a sliver element on the side, let's just
          // refine the side instead, if we're allowed.
          if (worst_cos < -0.6) // -0.5 is the best we could enforce?
            {
              side = worst_side;

              if (this->is_refine_boundary_allowed(boundary_info,
                                                   *cavity_elem,
                                                   side))
                {
                  // Let's just try bisecting for now
                  new_pt = (cavity_elem->point(side) +
                            cavity_elem->point((side+1)%3)) / 2;
                  new_node = add_point(new_pt);
                  boundary_refine(side);
                }
              else // Do the best we can under these restrictions
                {
                  new_pt = cavity_elem->vertex_average();
                  new_node = add_point(new_pt);

                  // This was going to be a side refinement but it's
                  // now an internal refinement
                  side = invalid_uint;
                }
            }
          else
            new_node = add_point(new_pt);
        }
      else
        libmesh_assert(new_node);

      // Find the Delaunay cavity around the new point.
      std::set<Elem *, decltype(comp)> cavity(comp);

      std::set<Elem *, decltype(comp)> unchecked_cavity ({cavity_elem}, comp);
      while (!unchecked_cavity.empty())
        {
          std::set<Elem *, decltype(comp)> checking_cavity(comp);
          checking_cavity.swap(unchecked_cavity);
          for (Elem * checking_elem : checking_cavity)
            {
              for (auto s : make_range(3u))
                {
                  Elem * neigh = checking_elem->neighbor_ptr(s);
                  if (!neigh)
                    continue;
                  if (checking_cavity.find(neigh) !=
                      checking_cavity.end())
                    continue;
                  if (cavity.find(neigh) !=
                      cavity.end())
                    continue;

                  if (in_circumcircle(*neigh, new_pt, TOLERANCE*TOLERANCE))
                    unchecked_cavity.insert(neigh);
                }
            }

          libmesh_merge_move(cavity, checking_cavity);
        }

      // Retriangulate the Delaunay cavity.
      // Each of our cavity triangle edges that are exterior to the
      // cavity will be a source of one new triangle.

      // Set of elements that might need Delaunay swaps
      std::set<Elem *, decltype(comp)> check_delaunay_on(comp);

      // Keep maps for doing neighbor pointer assignment.  Not going
      // to iterate through these so hashing pointers is fine.
      std::unordered_map<Node *, std::pair<Elem *, boundary_id_type>>
        neighbors_CCW, neighbors_CW;

      for (Elem * old_elem : cavity)
        {
          for (auto s : make_range(3u))
            {
              Elem * neigh = old_elem->neighbor_ptr(s);
              if (!neigh || !cavity.count(neigh))
                {
                  Node * node_CW = old_elem->node_ptr(s),
                       * node_CCW = old_elem->node_ptr((s+1)%3);

                  auto set_neighbors =
                    [&neighbors_CW, &neighbors_CCW, &node_CW,
                      &node_CCW, &boundary_info]
                    (Elem * new_neigh, boundary_id_type bcid)
                  {
                    // Set clockwise neighbor and vice-versa if possible
                    auto CW_it = neighbors_CW.find(node_CW);
                    if (CW_it == neighbors_CW.end())
                      {
                        libmesh_assert(!neighbors_CCW.count(node_CW));
                        neighbors_CCW[node_CW] = std::make_pair(new_neigh, bcid);
                      }
                    else
                      {
                        Elem * neigh_CW = CW_it->second.first;
                        if (new_neigh)
                          {
                            new_neigh->set_neighbor(0, neigh_CW);
                            boundary_id_type bcid_CW = CW_it->second.second;
                            if (bcid_CW != BoundaryInfo::invalid_id)
                              boundary_info.add_side(new_neigh, 0, bcid_CW);

                          }
                        if (neigh_CW)
                          {
                            neigh_CW->set_neighbor(2, new_neigh);
                            if (bcid != BoundaryInfo::invalid_id)
                              boundary_info.add_side(neigh_CW, 2, bcid);
                          }
                        neighbors_CW.erase(CW_it);
                      }

                    // Set counter-CW neighbor and vice-versa if possible
                    auto CCW_it = neighbors_CCW.find(node_CCW);
                    if (CCW_it == neighbors_CCW.end())
                      {
                        libmesh_assert(!neighbors_CW.count(node_CCW));
                        neighbors_CW[node_CCW] = std::make_pair(new_neigh, bcid);
                      }
                    else
                      {
                        Elem * neigh_CCW = CCW_it->second.first;
                        if (new_neigh)
                          {
                            boundary_id_type bcid_CCW = CCW_it->second.second;
                            new_neigh->set_neighbor(2, neigh_CCW);
                            if (bcid_CCW != BoundaryInfo::invalid_id)
                              boundary_info.add_side(new_neigh, 2, bcid_CCW);
                          }
                        if (neigh_CCW)
                          {
                            neigh_CCW->set_neighbor(0, new_neigh);
                            if (bcid != BoundaryInfo::invalid_id)
                              boundary_info.add_side(neigh_CCW, 0, bcid);
                          }
                        neighbors_CCW.erase(CCW_it);
                      }
                  };

                  // We aren't going to try to add a sliver element if we
                  // have a new boundary node here.  We do need to
                  // keep track of other elements' neighbors, though.
                  if (old_elem == cavity_elem &&
                      s == side)
                    {
                      std::vector<boundary_id_type> bcids;
                      boundary_info.boundary_ids(old_elem, s, bcids);
                      libmesh_assert_equal_to(bcids.size(), 1);
                      set_neighbors(nullptr, bcids[0]);
                      continue;
                    }

                  auto new_elem = Elem::build_with_id(TRI3, ne++);
                  new_elem->set_node(0) = new_node;
                  new_elem->set_node(1) = node_CW;
                  new_elem->set_node(2) = node_CCW;
                  libmesh_assert(!new_elem->is_flipped());

                  // Set in-and-out-of-cavity neighbor pointers
                  new_elem->set_neighbor(1, neigh);
                  if (neigh)
                    {
                      const unsigned int neigh_s =
                        neigh->which_neighbor_am_i(old_elem);
                      neigh->set_neighbor(neigh_s, new_elem.get());
                    }
                  else
                    {
                      std::vector<boundary_id_type> bcids;
                      boundary_info.boundary_ids(old_elem, s, bcids);
                      boundary_info.add_side(new_elem.get(), 1, bcids);
                    }

                  // Set in-cavity neighbors' neighbor pointers
                  set_neighbors(new_elem.get(), BoundaryInfo::invalid_id);

                  // C++ allows function argument evaluation in any
                  // order, but we need get() to precede move
                  Elem * new_elem_ptr = new_elem.get();
                  new_elems.emplace(new_elem_ptr, std::move(new_elem));

                  check_delaunay_on.insert(new_elem_ptr);
                }
            }

          boundary_info.remove(old_elem);
        }

      // Now that we're done using our cavity elems (including with a
      // cavity.find() that used a comparator that dereferences the
      // elmeents!) it's safe to delete them.
      for (Elem * old_elem : cavity)
        {
          auto it = new_elems.find(old_elem);
          if (it == new_elems.end())
            mesh.delete_elem(old_elem);
          else
            new_elems.erase(it);
        }

      // Everybody found their match?
      libmesh_assert(neighbors_CW.empty());
      libmesh_assert(neighbors_CCW.empty());

      // Because we're preserving boundaries here, our naive cavity
      // triangulation might not be a Delaunay triangulation.  Let's
      // check and if necessary fix that; we depend on it when doing
      // future point insertions.
      restore_delaunay(check_delaunay_on, boundary_info);

      // This is too expensive to do on every cavity in devel mode
#ifdef DEBUG
      libmesh_assert_delaunay(mesh, new_elems);
#endif
    }

  // If we added any new boundary nodes, we're going to need to keep
  // track of the changes they made to the outer polyline and/or to
  // any holes.
  if (!next_boundary_node.empty())
    {
      auto checked_emplace = [this](dof_id_type new_first,
                                    dof_id_type new_second)
      {
#ifdef DEBUG
        for (auto [first, second] : this->segments)
          {
            libmesh_assert_not_equal_to(first, new_first);
            libmesh_assert_not_equal_to(second, new_second);
          }
        if (!this->segments.empty())
          libmesh_assert_equal_to(this->segments.back().second, new_first);
#endif
        libmesh_assert_not_equal_to(new_first, new_second);

        this->segments.emplace_back (new_first, new_second);
      };

      // Keep track of the outer polyline
      if (this->segments.empty())
        {
          dof_id_type last_id = DofObject::invalid_id;

          // Custom loop because we increment node_it 1+ times inside
          for (auto node_it = _mesh.nodes_begin(),
               node_end = _mesh.nodes_end();
               node_it != node_end;)
            {
              Node & node = **node_it;
              ++node_it;

              const dof_id_type node_id = node.id();

              // Don't add Steiner points
              if (node_id >= _n_boundary_nodes)
                break;

              // Connect up the previous node, if we didn't already
              // connect it after some newly inserted nodes
              if (!this->segments.empty())
                last_id = this->segments.back().second;

              if (last_id != DofObject::invalid_id &&
                  last_id != node_id)
                checked_emplace(last_id, node_id);

              last_id = node_id;

              // Connect to any newly inserted nodes
              Node * this_node = &node;
              auto it = next_boundary_node.find(*this_node);
              while (it != next_boundary_node.end())
                {
                  libmesh_assert(this_node->valid_id());
                  Node * next_node = it->second;
                  libmesh_assert(next_node->valid_id());

                  if (node_it != node_end &&
                      next_node == *node_it)
                    ++node_it;

                  checked_emplace(this_node->id(), next_node->id());

                  this_node = next_node;
                  if (this_node->id() == this->segments.front().first)
                    break;

                  it = next_boundary_node.find(*this_node);
                }
            }

          // We expect a closed loop here
          if (this->segments.back().second != this->segments.front().first)
            checked_emplace(this->segments.back().second,
                            this->segments.front().first);
        }
      else
        {
          std::vector<std::pair<unsigned int, unsigned int>> old_segments;
          old_segments.swap(this->segments);

          auto old_it  = old_segments.begin();

          Node * node = _mesh.node_ptr(old_it->first);
          Node * const first_node = node;

          do
            {
              const dof_id_type node_id = node->id();
              auto it = next_boundary_node.find(*node);
              if (it == next_boundary_node.end())
                {
                  while (node_id != old_it->first)
                    {
                      ++old_it;
                      libmesh_assert(old_it != old_segments.end());
                    }
                  node = mesh.node_ptr(old_it->second);
                }
              else
                {
                  node = it->second;
                }

              checked_emplace(node_id, node->id());
            }
          while (node != first_node);
        }

      // Keep track of any holes
      if (this->_holes)
        {
          // Do we have any holes that need to be newly replaced?
          for (const Hole * hole : *this->_holes)
            {
              if (this->replaced_holes.count(hole))
                continue;

              bool hole_point_insertion = false;
                for (auto p : make_range(hole->n_points()))
                  if (next_boundary_node.count(hole->point(p)))
                    {
                      hole_point_insertion = true;
                      break;
                    }
              if (hole_point_insertion)
                this->replaced_holes.emplace
                  (hole, std::make_unique<ArbitraryHole>(*hole));
            }

          // If we have any holes that are being replaced, make sure
          // their replacements are up to date.
          for (const Hole * hole : *this->_holes)
            {
              auto hole_it = replaced_holes.find(hole);
              if (hole_it == replaced_holes.end())
                continue;

              ArbitraryHole & arb = *hole_it->second;

              // We only need to update a replacement that's just had
              // new points inserted
              bool point_inserted = false;
              for (const Point & point : arb.get_points())
                if (next_boundary_node.count(point))
                  {
                    point_inserted = true;
                    break;
                  }

              if (!point_inserted)
                continue;

              // Find all points in the replacement hole
              std::vector<Point> new_points;

              // Our outer polyline is expected to have points in
              // counter-clockwise order, so it proceeds "to the left"
              // from the point of view of rays inside the domain
              // pointing outward, and our next_boundary_node ordering
              // was filled accordingly.
              //
              // Our inner holes are expected to have points in
              // counter-clockwise order, but for holes "to the left
              // as viewed from the hole interior is the *opposite* of
              // "to the left as viewed from the domain interior".  We
              // need to build the updated hole ordering "backwards".

              // We should never see duplicate points when we add one
              // to a hole; if we do then we did something wrong.
              auto push_back_new_point = [&new_points](const Point & p) {
                // O(1) assert in devel
                libmesh_assert(new_points.empty() ||
                               new_points.back() != p);
#ifdef DEBUG
                // O(N) asserts in dbg
                for (auto old_p : new_points)
                  libmesh_assert_not_equal_to(old_p, p);
#endif
                new_points.push_back(p);
              };

              for (auto point_it = arb.get_points().rbegin(),
                   point_end = arb.get_points().rend();
                   point_it != point_end;)
                {
                  Point point = *point_it;
                  ++point_it;

                  if (new_points.empty() ||
                      (point != new_points.back() &&
                       point != new_points.front()))
                    push_back_new_point(point);

                  auto it = next_boundary_node.find(point);
                  while (it != next_boundary_node.end())
                    {
                      point = *it->second;
                      if (point == new_points.front())
                        break;
                      if (point_it != point_end &&
                          point == *point_it)
                        ++point_it;
                      push_back_new_point(point);
                      it = next_boundary_node.find(point);
                    }
                }

              std::reverse(new_points.begin(), new_points.end());

              arb.set_points(std::move(new_points));
            }
        }
    }

  // Okay, *now* we can add the new elements.
  for (auto & [raw_elem, unique_elem] : new_elems)
    {
      libmesh_assert_equal_to(raw_elem, unique_elem.get());
      libmesh_assert(!raw_elem->is_flipped());
      libmesh_ignore(raw_elem); // Old gcc warns "unused variable"
      mesh.add_elem(std::move(unique_elem));
    }

  // Did we add anything?
  return !new_elems.empty();
}


bool Poly2TriTriangulator::should_refine_elem(Elem & elem)
{
  const Real min_area_target = this->desired_area();
  FunctionBase<Real> * area_func = this->get_desired_area_function();

  // If this isn't a question, why are we here?
  libmesh_assert(min_area_target > 0 ||
                 area_func != nullptr ||
                 this->has_auto_area_function());

  const Real area = elem.volume();

  // If we don't have position-dependent area targets we can make a
  // decision quickly
  if (!area_func && !this->has_auto_area_function())
    return (area > min_area_target);

  // If we do?
  //
  // See if we're meeting the local area target at all the elem
  // vertices first
  for (auto v : make_range(elem.n_vertices()))
    {
      // If we have an auto area function, we'll use it and override other area options
      const Real local_area_target = this->has_auto_area_function() ? this->get_auto_desired_area(elem.point(v)) : (*area_func)(elem.point(v));
      libmesh_error_msg_if
        (local_area_target <= 0,
         "Non-positive desired element areas are unachievable");
      if (area > local_area_target)
        return true;
    }

  // If our vertices are happy, it's still possible that our interior
  // isn't.  Are we allowed not to bother checking it?
  if (!min_area_target)
    return false;

  libmesh_not_implemented_msg
    ("Combining a minimum desired_area with an area function isn't yet supported.");
}


} // namespace libMesh







#endif // LIBMESH_HAVE_POLY2TRI
