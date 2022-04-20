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


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_POLY2TRI

// libmesh includes
#include "libmesh/poly2tri_triangulator.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/function_base.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/mesh_triangle_holes.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h"

// poly2tri includes
#include "poly2tri/poly2tri.h"

// Anonymous namespace - poly2tri doesn't define operator<(Point,Point)
namespace
{
struct P2TPointCompare
{
  bool operator()(const p2t::Point & a, const p2t::Point & b) const
  {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  }
};

bool in_circumcircle(const libMesh::Elem & elem,
                     const libMesh::Point & p)
{
  using namespace libMesh;

  libmesh_assert_equal_to(elem.n_vertices(), 3);

  const Point pv0 = elem.point(0) - p;
  const Point pv1 = elem.point(1) - p;
  const Point pv2 = elem.point(2) - p;

  return ((pv0.norm_sq() * (pv1(0)*pv2(1)-pv2(0)*pv1(1))) -
          (pv1.norm_sq() * (pv0(0)*pv2(1)-pv2(0)*pv0(1))) +
          (pv2.norm_sq() * (pv0(0)*pv1(1)-pv1(0)*pv0(1)))) > 0;
}

unsigned int segment_intersection(const libMesh::Elem & elem,
                                  libMesh::Point & source,
                                  const libMesh::Point & target,
                                  unsigned int source_side)
{
  using namespace libMesh;

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

      if (t < -TOLERANCE || t > 1 + TOLERANCE)
        continue;

      const Real u_num = targetsdx * edgedy - targetsdy * edgedx;
      const Real u = u_num * one_over_denom;

      if (u < -TOLERANCE || u > 1 + TOLERANCE)
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
    _serializer(_mesh),
    _n_boundary_nodes(n_boundary_nodes),
    _refine_bdy_allowed(true)
{
}


Poly2TriTriangulator::~Poly2TriTriangulator() = default;


// Primary function responsible for performing the triangulation
void Poly2TriTriangulator::triangulate()
{
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

  // Insert additional new points in between existing boundary points,
  // if that is requested and reasonable
  this->insert_any_extra_boundary_points();

  do
    {
      this->triangulate_current_points();
    }
  while (!this->insert_refinement_points());

  // Okay, we really do need to support boundary ids soon, but we
  // don't yet
  if (_markers)
    libmesh_not_implemented();

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


void Poly2TriTriangulator::triangulate_current_points()
{
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

  // Unless we're using an explicit segments list, we assume node ids
  // are contiguous here.
  if (this->segments.empty())
    libmesh_error_msg_if
      (_mesh.n_nodes() != _mesh.max_node_id(),
       "Poly2TriTriangulator needs contiguous node ids or explicit segments!");

  // And if we have more nodes than outer boundary points, the rest
  // may be interior "Steiner points".  We use a set here so we can
  // cheaply search and erase from it later, when we're identifying
  // hole points.
  std::set<p2t::Point, P2TPointCompare> steiner_points;

  // If we have any elements, we assume they come from a preceding
  // triangulation, and we clear them.
  _mesh.clear_elems();

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

  // If we have no explicit segments defined, our nodal id ordering
  // defines our outer polyline ordering:
  if (this->segments.empty())
    {
      for (auto & node : _mesh.node_ptr_range())
        {
          p2t::Point pt((*node)(0), (*node)(1));

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
            libmesh_not_implemented();

          point_node_map.emplace(*pt, node);
        }

        libmesh_error_msg_if(last_id != this->segments[0].first,
                             "Non-closed-loop triangulator segments");

      // If we have points that aren't in any segments, those will be
      // Steiner points
      for (auto & node : _mesh.node_ptr_range())
        {
          p2t::Point pt {(*node)(0), (*node)(1)};

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

  // Create poly2tri triangulator with our mesh points
  std::vector<p2t::Point *>
    outer_boundary_pointers(outer_boundary_points.size());
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
          poly2tri_hole.emplace_back(p(0), p(1));

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
              Node * node = _mesh.add_point(p);
              point_node_map[pt] = node;
            }
        }

      auto & poly2tri_ptrs = inner_hole_pointers[h];
      poly2tri_ptrs.resize(poly2tri_hole.size());

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

  for (auto ptri_ptr : triangles)
    {
      p2t::Triangle & ptri = *ptri_ptr;

      // We always use TRI3 here, since that's what we have nodes for;
      // if we need a higher order we'll convert at the end.
      auto elem = Elem::build_with_id(_elem_type, next_id++);
      for (auto v : make_range(3))
        {
          const p2t::Point & vertex = *ptri.GetPoint(v);

          Node * node = point_node_map[vertex];
          libmesh_assert(node);
          elem->set_node(v) = node;
        }

      _mesh.add_elem(std::move(elem));
    }
}



bool Poly2TriTriangulator::insert_refinement_points()
{
  if (this->minimum_angle() != 0)
    libmesh_not_implemented();

  // We need neighbor pointers for ray casting and cavity finding
  UnstructuredMesh & mesh = dynamic_cast<UnstructuredMesh &>(this->_mesh);
  mesh.find_neighbors();

  if (this->desired_area() == 0 &&
      this->get_desired_area_function() == nullptr)
    return false;

  // We won't immediately add these, lest we invalidate iterators on a
  // ReplicatedMesh.  They'll still be in the mesh neighbor topology
  // for the purpose of doing Delaunay cavity stuff, so we need to
  // manage memory here, but there's no point in adding them to the
  // Mesh just to remove them again afterward when we hit up poly2tri.

  // We'll need to be able to remove new elems from new_elems, in
  // cases where a later refinement insertion has a not-yet-added
  // element in its cavity, so we'll use a map here to make searching
  // possible.
  std::unordered_map<Elem *, std::unique_ptr<Elem>> new_elems;

  // Map of which points follow which in the boundary polylines.  If
  // we have to add new boundary points, we'll use this to construct
  // an updated this->segments to retriangulate with.  If we have to
  // add new hole points, we'll use this to insert points into an
  // ArbitraryHole.
  std::unordered_map<Point, Node *> next_boundary_node;

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
              if (this->refine_boundary_allowed())
                {
                  new_pt = ray_start;
                  new_node = mesh.add_point(new_pt);
                  boundary_refine(side);
                }
              else
                {
                  new_pt = cavity_elem->vertex_average();
                  new_node = mesh.add_point(new_pt);
                }

              break;
            }

          source_side = neigh->which_neighbor_am_i(cavity_elem);
          cavity_elem = neigh;
          side = invalid_uint;
        }

      // If we're not on a boundary ... should we be?  We don't want
      // to create any sliver elements or confuse poly2tri or
      // anything.
      if (side == invalid_uint)
        {
          libmesh_assert(!new_node);

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

              if (this->refine_boundary_allowed())
                {
                  // Let's just try bisecting for now
                  new_pt = (cavity_elem->point(side) +
                            cavity_elem->point((side+1)%3)) / 2;
                  new_node = mesh.add_point(new_pt);
                  boundary_refine(side);
                }
              else // Do the best we can under these restrictions
                {
                  new_pt = cavity_elem->vertex_average();
                  new_node = mesh.add_point(new_pt);
                }
            }
          else
            new_node = mesh.add_point(new_pt);
        }
      else
        libmesh_assert(new_node);


      // Find the Delaunay cavity around the new point.
      std::set<Elem *> unchecked_cavity {cavity_elem};
      std::set<Elem *> cavity;
      std::set<Node *> cavity_nodes {cavity_elem->node_ptr(0),
                                     cavity_elem->node_ptr(1),
                                     cavity_elem->node_ptr(2)};
      while (!unchecked_cavity.empty())
        {
          std::set<Elem *> checking_cavity;
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

                  if (in_circumcircle(*neigh, new_pt))
                    {
                      unchecked_cavity.insert(neigh);
                      for (auto v : make_range(3u))
                        cavity_nodes.insert(neigh->node_ptr(v));
                    }
                }
            }

          libmesh_merge_move(cavity, checking_cavity);
        }

      // Retriangulate the Delaunay cavity.
      // Each of our cavity triangle edges that are exterior to the
      // cavity will be a source of one new triangle.

      // Keep maps for doing neighbor pointer assignment
      std::map<Node *, Elem *> neighbors_CCW, neighbors_CW;

      for (Elem * old_elem : cavity)
        {
          for (auto s : make_range(3u))
            {
              Elem * neigh = old_elem->neighbor_ptr(s);
              if (!neigh || !cavity.count(neigh))
                {
                  Node * node_CW = old_elem->node_ptr(s),
                       * node_CCW = old_elem->node_ptr((s+1)%3);

                  auto set_neighbors = [&neighbors_CW, &neighbors_CCW,
                                        &node_CW, &node_CCW](Elem * new_neigh)
                  {
                    // Set clockwise neighbor and vice-versa if possible
                    auto CW_it = neighbors_CW.find(node_CW);
                    if (CW_it == neighbors_CW.end())
                      {
                        libmesh_assert(!neighbors_CCW.count(node_CW));
                        neighbors_CCW[node_CW] = new_neigh;
                      }
                    else
                      {
                        Elem * neigh_CW = CW_it->second;
                        if (new_neigh)
                          new_neigh->set_neighbor(0, neigh_CW);
                        if (neigh_CW)
                          neigh_CW->set_neighbor(2, new_neigh);
                        neighbors_CW.erase(CW_it);
                      }

                    // Set counter-CW neighbor and vice-versa if possible
                    auto CCW_it = neighbors_CCW.find(node_CCW);
                    if (CCW_it == neighbors_CCW.end())
                      {
                        libmesh_assert(!neighbors_CW.count(node_CCW));
                        neighbors_CW[node_CCW] = new_neigh;
                      }
                    else
                      {
                        Elem * neigh_CCW = CCW_it->second;
                        if (new_neigh)
                          new_neigh->set_neighbor(2, neigh_CCW);
                        if (neigh_CCW)
                          neigh_CCW->set_neighbor(0, new_neigh);
                        neighbors_CCW.erase(CCW_it);
                      }
                  };

                  // We aren't going to try to add a sliver element if we
                  // have a new boundary node here.  We do need to
                  // keep track of other elements' neighbors, though.
                  if (old_elem == cavity_elem &&
                      s == side)
                    {
                      set_neighbors(nullptr);
                      continue;
                    }

                  auto new_elem = Elem::build(TRI3);
                  new_elem->set_node(0) = new_node;
                  new_elem->set_node(1) = node_CW;
                  new_elem->set_node(2) = node_CCW;

                  // Set in-and-out-of-cavity neighbor pointers
                  new_elem->set_neighbor(1, neigh);
                  if (neigh)
                    {
                      const unsigned int neigh_s =
                        neigh->which_neighbor_am_i(old_elem);
                      neigh->set_neighbor(neigh_s, new_elem.get());
                    }

                  // Set in-cavity neighbors' neighbor pointers
                  set_neighbors(new_elem.get());

                  // C++ allows function argument evaluation in any
                  // order, but we need get() to precede move
                  Elem * new_elem_ptr = new_elem.get();
                  new_elems.emplace(new_elem_ptr, std::move(new_elem));
                }
            }

          auto it = new_elems.find(old_elem);
          if (it == new_elems.end())
            mesh.delete_elem(old_elem);
          else
            new_elems.erase(it);
        }

      // Everybody found their match?
      libmesh_assert(neighbors_CW.empty());
      libmesh_assert(neighbors_CCW.empty());
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

              for (auto point_it = arb.get_points().rbegin(),
                   point_end = arb.get_points().rend();
                   point_it != point_end;)
                {
                  Point point = *point_it;
                  ++point_it;

                  if (new_points.empty() ||
                      (point != new_points.back() &&
                       point != new_points.front()))
                    new_points.push_back(point);

                  auto it = next_boundary_node.find(point);
                  while (it != next_boundary_node.end())
                    {
                      point = *it->second;
                      if (point == new_points.front())
                        break;
                      if (point_it != point_end &&
                          point == *point_it)
                        ++point_it;
                      new_points.push_back(point);
                      it = next_boundary_node.find(point);
                    }
                }

              std::reverse(new_points.begin(), new_points.end());

              arb.set_points(std::move(new_points));
            }
        }
    }

  // Okay, *now* we can add the new elements.
  for (auto & new_elem_pair : new_elems)
    mesh.add_elem(std::move(new_elem_pair.second));

  // Did we add anything?
  return new_elems.empty();
}


bool Poly2TriTriangulator::should_refine_elem(Elem & elem)
{
  const Real min_area_target = this->desired_area();
  FunctionBase<Real> * area_func = this->get_desired_area_function();

  // If this isn't a question, why are we here?
  libmesh_assert(min_area_target > 0 ||
                 area_func != nullptr);

  const Real area = elem.volume();

  // If we don't have position-dependent area targets we can make a
  // decision quickly
  if (!area_func)
    return (area > min_area_target);

  // If we do?
  //
  // See if we're meeting the local area target at all the elem
  // vertices first
  Real local_area_target = (*area_func)(elem.point(0));
  for (auto v : make_range(1u, elem.n_vertices()))
    {
      if (area > local_area_target)
        return true;
      local_area_target =
        std::min(local_area_target,
                 (*area_func)(elem.point(v)));
    }

  if (area > local_area_target)
    return true;

  // If our vertices are happy, it's still possible that our interior
  // isn't.  Are we allowed not to bother checking it?
  if (!min_area_target)
    return false;

  libmesh_not_implemented();
}


} // namespace libMesh







#endif // LIBMESH_HAVE_POLY2TRI
