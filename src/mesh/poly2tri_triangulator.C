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
}

namespace libMesh
{

//
// Function definitions for the Poly2TriTriangulator class
//

// Constructor
Poly2TriTriangulator::Poly2TriTriangulator(UnstructuredMesh & mesh)
  : TriangulatorInterface(mesh),
    _serializer(_mesh)
{}




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

  // We currently require boundary nodes to be ordered by node id; we
  // don't yet support explicit segments
  if (this->segments.size())
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

  // Will the triangulation have holes?
  const std::size_t n_holes = _holes != nullptr ? _holes->size() : 0;

  // Insert additional new points in between existing boundary points,
  // if that is requested and reasonable
  this->insert_any_extra_boundary_points();

  // Mapping from Poly2Tri points to libMesh nodes, so we can get the
  // connectivity translated back later.
  std::map<const p2t::Point, Node *, P2TPointCompare> point_node_map;


  // Poly2Tri data structures
  // Poly2Tri takes vectors of pointers-to-Point for some reason, but
  // we'll just make those shims to vectors of Point rather than
  // individually/manually heap allocating everything.
  std::vector<p2t::Point> outer_boundary_points;
  std::vector<std::vector<p2t::Point>> inner_hole_points(n_holes);

  for (auto & node : _mesh.node_ptr_range())
    {
      outer_boundary_points.emplace_back((*node)(0), (*node)(1));

      const auto & pt = outer_boundary_points.back();

      // We're not going to support overlapping nodes yet.
      if (point_node_map.count(pt))
        libmesh_not_implemented();

      point_node_map.emplace(pt, node);
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
      const auto & our_hole = *((*_holes)[h]);
      auto & poly2tri_hole = inner_hole_points[h];
      for (auto i : make_range(our_hole.n_points()))
        {
          Point p = our_hole.point(i);
          poly2tri_hole.emplace_back(p(0), p(1));

          const auto & pt = poly2tri_hole.back();

          // We're not going to support overlapping nodes or
          // non-simple-manifold topologies yet
          if (point_node_map.count(pt))
            libmesh_not_implemented();

          Node * node = _mesh.add_point(p);
          point_node_map[pt] = node;
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

} // namespace libMesh







#endif // LIBMESH_HAVE_POLY2TRI
