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
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/mesh_triangle_holes.h"
#include "libmesh/poly2tri_triangulator.h"
#include "libmesh/unstructured_mesh.h"

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
                                  const libMesh::Point & target)
{
  using namespace libMesh;

  libmesh_assert_equal_to(elem.dim(), 2);

  const auto ns = elem.n_sides();

  for (auto s : make_range(ns))
    {
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

      const Real sourcesdx = source(0)-v0(0),
                 sourcesdy = source(1)-v0(1);

      const Real t_num = sourcesdx * raydy -
                         sourcesdy * raydx;
      const Real t = t_num * one_over_denom;

      if (t < -TOLERANCE || t > 1 + TOLERANCE)
        continue;

      const Real u_num = sourcesdx * edgedy - sourcesdy * raydx;
      const Real u = u_num * one_over_denom;

      if (u < -TOLERANCE || u > 1 + TOLERANCE)
        continue;

      source(0) += raydx * u;
      source(1) += raydy * u;
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
    _n_boundary_nodes(n_boundary_nodes)
{
}




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

  // We assume node ids are contiguous here.
  libmesh_assert_equal_to(_mesh.n_nodes(), _mesh.max_node_id());

  // And if we have more nodes than boundary points, the rest will be
  // interior "Steiner points"
  std::vector<p2t::Point> steiner_points;

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
          p2t::Point * pt;

          // If we're out of boundary nodes, the rest are going to be
          // Steiner points
          if (node->id() < _n_boundary_nodes)
            {
              outer_boundary_points.emplace_back((*node)(0), (*node)(1));
              pt = &outer_boundary_points.back();
            }
          else
            {
              steiner_points.emplace_back((*node)(0), (*node)(1));
              pt = &steiner_points.back();
            }

          // We're not going to support overlapping nodes on the boundary
          if (point_node_map.count(*pt))
            libmesh_not_implemented();

          point_node_map.emplace(*pt, node);
        }
    }
  // If we have explicit segments defined, that's our outer polyline
  // ordering:
  else
    {
      // Let's make sure our segments are in order
      dof_id_type last_id = DofObject::invalid_id;

      // Add nodes from every segment, in order, to the outer polyline
      for (auto segment : this->segments)
        {
          const dof_id_type segment_start = segment.first,
                            segment_end = segment.second;

          libmesh_error_msg_if(segment_start == DofObject::invalid_id,
                               "Bad triangulator segment start");
          libmesh_error_msg_if(segment_end == DofObject::invalid_id,
                               "Bad triangulator segment end");

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
              steiner_points.push_back(pt);
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

  // Add any steiner points.
  for (auto & p : steiner_points)
    cdt.AddPoint(&p);

  // Add any holes
  for (auto h : make_range(n_holes))
    {
      const auto & our_hole = *((*_holes)[h]);
      auto & poly2tri_hole = inner_hole_points[h];
      for (auto i : make_range(our_hole.n_points()))
        {
          Point p = our_hole.point(i);
          poly2tri_hole.emplace_back(p(0), p(1));

          const auto & pt = poly2tri_hole.back();

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

  const Real area_target = this->desired_area();

  if (area_target == 0)
    return false;

  // We won't actually add these, lest we invalidate iterators on a
  // ReplicatedMesh.  They'll still be in the mesh neighbor topology
  // for the purpose of doing Delaunay cavity stuff, so we need to
  // manage memory here, but there's no point in adding them to the
  // Mesh just to remove them again afterward when we hit up poly2tri.
  std::vector<std::unique_ptr<Elem>> new_elems;

  for (auto & elem : mesh.element_ptr_range())
    {
      // element_ptr_range skips deleted elements ... right?
      libmesh_assert(elem);

      // We only handle triangles in our triangulation
      libmesh_assert_equal_to(elem->level(), 0u);
      libmesh_assert_equal_to(elem->type(), TRI3);

      const Real area = elem->volume();

      // If this triangle is as small as we desire, move along
      if (area <= area_target)
        continue;

      // Otherwise add a Steiner point.  We'd like to add the
      // circumcenter ...
      Point new_pt = elem->quasicircumcenter();

      // But that might be outside our triangle, or even outside the
      // boundary.  Let's find a triangle containing it, moving it to
      // the boundary if we have to.
      Elem * cavity_elem = elem; // Maybe?  Start looking there anyway
      unsigned int side = invalid_uint;
      Point ray_start = elem->vertex_average();
      while (!cavity_elem->contains_point(new_pt))
        {
          side = segment_intersection(*cavity_elem, ray_start, new_pt);

          libmesh_assert_not_equal_to (side, invalid_uint);

          Elem * neigh = cavity_elem->neighbor_ptr(side);
          if (!neigh)
            {
              new_pt = ray_start;
              break;
            }

          cavity_elem = neigh;
          side = invalid_uint;
        }

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

          cavity.insert(checking_cavity.begin(),
                        checking_cavity.end());
        }

      // Retriangulate the Delaunay cavity.
      // Each of our cavity triangle edges that are exterior to the
      // cavity will be a source of one new triangle.
      Node * new_node = mesh.add_point(new_pt);

      // Keep maps for doing neighbor pointer assignment
      std::map<Node *, Elem *> neighbors_CCW, neighbors_CW;

      for (Elem * old_elem : cavity)
        {
          for (auto s : make_range(3u))
            {
              Elem * neigh = old_elem->neighbor_ptr(s);
              if (neigh && !cavity.count(neigh))
                {
                  auto new_elem = Elem::build(TRI3);
                  new_elem->set_node(0) = new_node;
                  Node * node_CW = old_elem->node_ptr(s);
                  new_elem->set_node(1) = node_CW;
                  Node * node_CCW = old_elem->node_ptr((s+1)%3);
                  new_elem->set_node(2) = node_CCW;

                  // Set neighbor pointers
                  new_elem->set_neighbor(1, neigh);
                  const unsigned int neigh_s =
                    neigh->which_neighbor_am_i(old_elem);
                  neigh->set_neighbor(neigh_s, new_elem.get());

                  // Set clockwise neighbor and vice-versa if possible
                  {
                    auto it = neighbors_CW.find(node_CW);
                    if (it == neighbors_CW.end())
                      {
                        libmesh_assert(!neighbors_CCW.count(node_CW));
                        neighbors_CCW[node_CW] = new_elem.get();
                      }
                    else
                      {
                        Elem * neigh_CW = it->second;
                        new_elem->set_neighbor(0, neigh_CW);
                        neigh_CW->set_neighbor(2, new_elem.get());
                        neighbors_CW.erase(it);
                      }
                  }

                  // Set counter-CW neighbor and vice-versa if possible
                  {
                    auto it = neighbors_CCW.find(node_CCW);
                    if (it == neighbors_CCW.end())
                      {
                        libmesh_assert(!neighbors_CW.count(node_CCW));
                        neighbors_CW[node_CCW] = new_elem.get();
                      }
                    else
                      {
                        Elem * neigh_CCW = it->second;
                        new_elem->set_neighbor(2, neigh_CCW);
                        neigh_CCW->set_neighbor(0, new_elem.get());
                        neighbors_CCW.erase(it);
                      }

                    new_elems.push_back(std::move(new_elem));
                  }
                }
            }

          mesh.delete_elem(old_elem);
        }

      // Everybody found their match?
      libmesh_assert(neighbors_CW.empty());
      libmesh_assert(neighbors_CCW.empty());
    }

  // Did we add anything?
  return new_elems.empty();
}


} // namespace libMesh







#endif // LIBMESH_HAVE_POLY2TRI
